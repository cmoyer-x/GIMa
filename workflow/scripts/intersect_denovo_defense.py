"""
intersect_denovo_defense.py

Intersects de novo genomic islands (found by MAbsIslandScanner but NOT
anchored to PPanGGOLiN RGPs) with DefenseFinder system coordinates to
identify defense systems in fixed islands — acquisitions present in too
many strains to appear as RGPs but still carrying detectable HGT evidence.

These are the defense systems your original 83 RGP-based analysis could
not see.

Usage:
    python intersect_denovo_defense.py \
        --islands   island_predictions/all_islands_combined.tsv \
        --defense   all_genomes_defense_finder_systems.tsv \
        --gff_dir   all_gffs \
        --outfile   denovo_defense_intersection.tsv \
        --report    denovo_defense_report.html \
        --min_conf  moderate
"""

import os, csv, re, json, argparse
from collections import defaultdict, Counter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--islands",   required=True,
                   help="all_islands_combined.tsv from MAbsIslandScanner")
    p.add_argument("--defense",   required=True,
                   help="DefenseFinder systems TSV (all_genomes_defense_finder_systems.tsv)")
    p.add_argument("--gff_dir",   required=True,
                   help="Directory of Prokka GFF files for coordinate lookup")
    p.add_argument("--outfile",   default="denovo_defense_intersection.tsv")
    p.add_argument("--report",    default="denovo_defense_report.html")
    p.add_argument("--min_conf",  default="moderate",
                   choices=["high","moderate","low"],
                   help="Minimum island confidence to consider (default: moderate)")
    return p.parse_args()


def canonical(s):
    return re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=re.IGNORECASE)


# ── Load de novo islands ──────────────────────────────────────────────────────
def load_denovo_islands(islands_file, min_conf):
    conf_rank = {"high": 3, "moderate": 2, "low": 1, "very_low": 0}
    min_rank  = conf_rank[min_conf]

    islands = []
    with open(islands_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # De novo = no RGP seed AND meets confidence threshold
            if row["rgp_seed"]:
                continue
            if conf_rank.get(row["confidence"], 0) < min_rank:
                continue
            try:
                islands.append({
                    "strain":       row["strain"],
                    "contig":       row["contig"],
                    "start":        int(row["start"]),
                    "end":          int(row["end"]),
                    "length":       int(row["length"]),
                    "confidence":   row["confidence"],
                    "n_evidence":   int(row["n_evidence"]),
                    "evidence":     row["evidence"],
                    "island_score": float(row["island_score"]),
                    "gc_z":         float(row["gc_z"]) if row["gc_z"] else 0,
                    "cai_ratio":    float(row["cai_ratio"]) if row["cai_ratio"] else 1.0,
                    "age_estimate": row["age_estimate"],
                    "trna_flanked": row["trna_flanked"],
                    "trna_product": row["trna_product"],
                    "has_dr":       row["has_dr"],
                    "dr_seq":       row["dr_seq"],
                    "mob_types":    row["mob_types"],
                })
            except: continue

    print(f"  De novo islands loaded: {len(islands)} (conf >= {min_conf})")
    return islands


# ── Load DefenseFinder results ────────────────────────────────────────────────
def load_defense_systems(defense_file):
    """Load DefenseFinder output. Returns list of system dicts.
    Extracts strain from sys_id since there is no genome column.
    sys_id format: GD01_Wadjet_I_38 -> strain GD01
    """
    systems = []
    with open(defense_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Extract strain from sys_id: GD01_Wadjet_I_38 -> GD01
            sys_id = row.get("sys_id","")
            # Strain is everything before the first defense system name
            # sys_beg looks like GD01_1451 - strain is prefix before _NNNN
            sys_beg = row.get("sys_beg","")
            if sys_beg:
                # Extract strain: GD01_1451 -> GD01 (everything before last _NNNN)
                parts = sys_beg.rsplit("_", 1)
                strain = parts[0] if len(parts) == 2 else sys_beg
            else:
                # Fall back to sys_id parsing
                strain = re.sub(r'_[A-Z].*$', '', sys_id)
            row["genome"] = strain
            systems.append(row)
    print(f"  Defense systems loaded: {len(systems)}")
    return systems


# ── Load GFF for coordinate lookup ───────────────────────────────────────────
def load_gff_coords(gff_dir):
    """
    Build locus_tag -> (contig, start, end) mapping from all GFF files.
    Used to convert DefenseFinder locus_tag gene names to genome coordinates.
    """
    locus_coords = {}
    n_files = 0
    for fname in os.listdir(gff_dir):
        if not fname.endswith(".gff"): continue
        strain = fname.replace("_prokka.gff","").replace(".gff","")
        try:
            with open(os.path.join(gff_dir, fname)) as f:
                for line in f:
                    if "##FASTA" in line: break
                    if line.startswith("#") or "\tCDS\t" not in line: continue
                    parts = line.strip().split("\t")
                    if len(parts) < 9: continue
                    contig = parts[0].replace("gnl|XXX|","")
                    start  = int(parts[3]) - 1
                    end    = int(parts[4])
                    attrs  = {}
                    for a in parts[8].split(";"):
                        if "=" in a:
                            k, v = a.split("=",1)
                            attrs[k.strip()] = v.strip()
                    lt = attrs.get("locus_tag","")
                    if lt:
                        entry = {"strain": strain, "contig": contig,
                                 "start": start, "end": end}
                        locus_coords[lt] = entry
                        # Index all DefenseFinder locus tag variants:
                        # GFF GD01_14510 -> DF GD01_1451 (strip trailing 0)
                        # GFF GD02_01330 -> DF GD02_133  (strip leading 0 from numeric part)
                        # GFF GD02_00133 -> DF GD02_133  (strip all leading zeros)
                        if "_" in lt:
                            prefix, num = lt.rsplit("_", 1)
                            if num.isdigit():
                                # Strip trailing zero
                                if num.endswith("0"):
                                    locus_coords[f"{prefix}_{num[:-1]}"] = entry
                                # Strip leading zeros
                                stripped = num.lstrip("0") or "0"
                                locus_coords[f"{prefix}_{stripped}"] = entry
                                # Strip leading zeros then trailing zero
                                if stripped.endswith("0") and len(stripped) > 1:
                                    locus_coords[f"{prefix}_{stripped[:-1]}"] = entry
            n_files += 1
        except: continue
    print(f"  GFF coordinate index: {len(locus_coords):,} locus tags from {n_files} files")
    return locus_coords


def contigs_match(a, b):
    if a == b: return True
    a_base = re.sub(r'_\d+$','',a); b_base = re.sub(r'_\d+$','',b)
    return (canonical(a)==canonical(b) or a_base==b_base
            or a_base==b or a==b_base)


# ── Intersection ──────────────────────────────────────────────────────────────
def intersect(islands, defense_systems, locus_coords):
    """
    For each defense system, look up the genomic coordinates of its genes
    and check if they overlap with any de novo island.
    """
    # Index de novo islands by strain
    islands_by_strain = defaultdict(list)
    for isl in islands:
        islands_by_strain[isl["strain"]].append(isl)
        islands_by_strain[canonical(isl["strain"])].append(isl)

    results = []
    n_checked = 0
    n_coord_missing = 0

    for sys in defense_systems:
        strain  = sys.get("genome") or sys.get("replicon","")
        subtype = sys.get("subtype","")
        sys_beg = sys.get("sys_beg","")
        sys_end = sys.get("sys_end","")

        # Get coordinates for the system's first and last gene
        beg_info = locus_coords.get(sys_beg)
        end_info = locus_coords.get(sys_end)

        if not beg_info and not end_info:
            n_coord_missing += 1
            continue

        info   = beg_info or end_info
        contig = info["contig"]
        sys_start = min(
            beg_info["start"] if beg_info else end_info["start"],
            end_info["start"] if end_info else beg_info["start"]
        )
        sys_stop = max(
            beg_info["end"] if beg_info else end_info["end"],
            end_info["end"] if end_info else beg_info["end"]
        )
        n_checked += 1

        # Check overlap with de novo islands for this strain
        strain_islands = (islands_by_strain.get(strain) or
                         islands_by_strain.get(canonical(strain)) or [])

        hit_island = None
        for isl in strain_islands:
            if not contigs_match(isl["contig"], contig):
                continue
            # Overlap check
            if sys_start <= isl["end"] and sys_stop >= isl["start"]:
                hit_island = isl
                break

        in_denovo = hit_island is not None

        results.append({
            "strain":          strain,
            "sys_id":          sys.get("sys_id",""),
            "type":            sys.get("type",""),
            "subtype":         subtype,
            "sys_beg":         sys_beg,
            "sys_end":         sys_end,
            "sys_start":       sys_start,
            "sys_stop":        sys_stop,
            "in_denovo_island":  "Yes" if in_denovo else "No",
            "island_start":    hit_island["start"]        if hit_island else "",
            "island_end":      hit_island["end"]          if hit_island else "",
            "island_length":   hit_island["length"]       if hit_island else "",
            "island_confidence":hit_island["confidence"]  if hit_island else "",
            "island_evidence": hit_island["evidence"]     if hit_island else "",
            "island_score":    hit_island["island_score"] if hit_island else "",
            "age_estimate":    hit_island["age_estimate"] if hit_island else "",
            "cai_ratio":       hit_island["cai_ratio"]    if hit_island else "",
            "trna_flanked":    hit_island["trna_flanked"] if hit_island else "",
            "has_dr":          hit_island["has_dr"]       if hit_island else "",
            "mob_types":       hit_island["mob_types"]    if hit_island else "",
        })

    print(f"  Systems checked:       {n_checked}")
    print(f"  Missing coordinates:   {n_coord_missing}")
    return results


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    print("\nLoading data...")
    islands      = load_denovo_islands(args.islands, args.min_conf)
    defense_sys  = load_defense_systems(args.defense)
    locus_coords = load_gff_coords(args.gff_dir)

    print("\nIntersecting de novo islands with defense systems...")
    results = intersect(islands, defense_sys, locus_coords)

    # Write TSV
    fieldnames = [
        "strain","sys_id","type","subtype","sys_beg","sys_end",
        "sys_start","sys_stop",
        "in_denovo_island","island_start","island_end","island_length",
        "island_confidence","island_evidence","island_score",
        "age_estimate","cai_ratio","trna_flanked","has_dr","mob_types"
    ]
    with open(args.outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"\nWritten: {args.outfile}")

    # ── Summary ───────────────────────────────────────────────────────────────
    in_denovo = [r for r in results if r["in_denovo_island"]=="Yes"]
    not_in    = [r for r in results if r["in_denovo_island"]=="No"]

    # De novo defense systems by type
    sys_types = Counter(r["subtype"] for r in in_denovo)
    age_dist  = Counter(r["age_estimate"] for r in in_denovo)
    conf_dist = Counter(r["island_confidence"] for r in in_denovo)

    print(f"\n{'='*60}")
    print(f"  Total defense systems:               {len(results)}")
    print(f"  In de novo island (fixed):            {len(in_denovo)}")
    print(f"  Not in de novo island:                {len(not_in)}")
    pct = 100*len(in_denovo)/len(results) if results else 0
    print(f"  % of defense systems in fixed islands:{pct:.1f}%")
    print(f"{'='*60}")

    print(f"\nTop defense system types in fixed islands:")
    for sys, n in sys_types.most_common(15):
        pct = 100*n/len(in_denovo) if in_denovo else 0
        print(f"  {sys:30}: {n:4} ({pct:.1f}%)")

    print(f"\nAge distribution of fixed defense islands:")
    for age in ["very_recent","recent","moderate","old"]:
        print(f"  {age:15}: {age_dist.get(age,0):4}")

    print(f"\nConfidence distribution:")
    for conf in ["high","moderate","low"]:
        print(f"  {conf:10}: {conf_dist.get(conf,0):4}")

    # Previously known vs newly found
    print(f"\nComparison to original RGP-based analysis:")
    print(f"  Original defense RGPs (PPanGGOLiN):  86")
    print(f"  New fixed defense islands (de novo):  {len(in_denovo)}")
    print(f"  Total defense islands:                {86 + len(in_denovo)}")
    print(f"  Increase:                             {100*len(in_denovo)/86:.0f}% more than RGP-only")

    # Unique strains with fixed defense islands
    strains_with_fixed = set(r["strain"] for r in in_denovo)
    print(f"\nStrains with fixed defense islands:    {len(strains_with_fixed)}")

    build_html(results, in_denovo, sys_types, age_dist, args.report, args.min_conf)
    print(f"Report: {args.report}")


def build_html(results, in_denovo, sys_types, age_dist, outpath, min_conf='moderate'):
    import json
    from collections import Counter

    conf_dist = Counter(r["island_confidence"] for r in in_denovo)
    top_def   = sorted(in_denovo,
                       key=lambda r: float(r["island_score"] or 0),
                       reverse=True)[:30]

    data = {
        "n_total":      len(results),
        "n_denovo":     len(in_denovo),
        "n_not":        len(results) - len(in_denovo),
        "sys_types":    dict(sys_types.most_common(12)),
        "age_dist":     dict(age_dist),
        "conf_dist":    dict(conf_dist),
        "top_def":      top_def,
        "n_strains":    len(set(r["strain"] for r in in_denovo)),
        "n_original":   86,
    }
    data_json = json.dumps(data)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>De novo defense island intersection</title>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@300;400;500&display=swap" rel="stylesheet">
<style>
  :root {{--bg:#f8f9fb;--surface:#fff;--border:#dde1ea;--text:#1a1f2e;--muted:#5a6380;
    --accent:#0e7a5a;--font-mono:'IBM Plex Mono',monospace;--font-sans:'IBM Plex Sans',sans-serif}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{background:var(--bg);color:var(--text);font-family:var(--font-sans);font-size:14px;padding:2rem}}
  h1{{font-size:16px;font-weight:500;color:var(--accent);font-family:var(--font-mono);
    text-transform:uppercase;letter-spacing:.08em;margin-bottom:.25rem}}
  h2{{font-size:13px;font-weight:500;margin:1.5rem 0 .75rem}}
  .sub{{font-size:12px;color:var(--muted);font-family:var(--font-mono);margin-bottom:1.5rem}}
  .stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-bottom:1.5rem}}
  .stat{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:.75rem 1rem}}
  .sv{{font-size:22px;font-weight:500;font-family:var(--font-mono)}}
  .sl{{font-size:10px;text-transform:uppercase;letter-spacing:.07em;color:var(--muted)}}
  .charts{{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-bottom:1.5rem}}
  .cw{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:1rem}}
  .ct{{font-size:11px;text-transform:uppercase;letter-spacing:.07em;color:var(--muted);margin-bottom:.75rem}}
  table{{width:100%;border-collapse:collapse;font-size:12px;font-family:var(--font-mono)}}
  th{{background:var(--surface);border-bottom:2px solid var(--border);padding:8px 10px;
    text-align:left;font-weight:500;font-size:11px;text-transform:uppercase;
    letter-spacing:.06em;color:var(--muted)}}
  td{{padding:7px 10px;border-bottom:1px solid var(--border)}}
  tr:hover td{{background:#f0faf6}}
  .badge{{display:inline-block;font-size:10px;padding:1px 7px;border-radius:3px;font-weight:600}}
  .bh{{background:#d4f2e8;color:#0e7a5a}}.bm{{background:#ddeeff;color:#1a5fa8}}
  .bl{{background:#eaecf0;color:#5a6380}}.bvr{{background:#fde8e4;color:#a82c1a}}
  .note{{background:#f0faf6;border:1px solid #9fe1cb;border-radius:8px;
    padding:1rem;font-size:12px;color:var(--muted);margin-top:1.5rem;line-height:1.7}}
</style>
</head>
<body>
<h1>Fixed defense islands — de novo intersection</h1>
<p class="sub">Defense systems in genomic islands invisible to PPanGGOLiN · MAbsIslandScanner de novo predictions · M. abscessus</p>

<div class="stats">
  <div class="stat"><div class="sv">{len(results):,}</div><div class="sl">Total defense systems</div></div>
  <div class="stat"><div class="sv" style="color:#0e7a5a">{len(in_denovo):,}</div>
    <div class="sl">In fixed (de novo) island</div></div>
  <div class="stat"><div class="sv">{86 + len(in_denovo):,}</div>
    <div class="sl">Total defense islands (RGP + fixed)</div></div>
  <div class="stat"><div class="sv">{len(set(r["strain"] for r in in_denovo))}</div>
    <div class="sl">Strains with fixed defense islands</div></div>
</div>

<div class="charts">
  <div class="cw">
    <p class="ct">Defense system types in fixed islands (top 12)</p>
    <div style="position:relative;height:280px">
      <canvas id="c1" role="img" aria-label="Horizontal bar chart of defense system types found in fixed de novo islands"></canvas>
    </div>
  </div>
  <div class="cw">
    <p class="ct">Age estimate of fixed defense islands</p>
    <div style="position:relative;height:280px">
      <canvas id="c2" role="img" aria-label="Bar chart of age estimates for fixed defense islands"></canvas>
    </div>
  </div>
</div>

<h2>Top fixed defense islands — ranked by island score</h2>
<table>
  <thead>
    <tr><th>Strain</th><th>Defense system</th><th>Island coords</th><th>Confidence</th>
        <th>Evidence</th><th>Age</th><th>CAI</th><th>tRNA</th><th>DR</th></tr>
  </thead>
  <tbody id="defTbl"></tbody>
</table>

<div class="note">
  <strong>What these represent:</strong>
  Fixed defense islands are genomic islands present in too many strains to appear as
  regions of genomic plasticity (RGPs) in PPanGGOLiN. They were acquired by HGT in an
  ancestral strain and have since spread through the clinical population — or have been
  present so long they are now fixed. These islands are invisible to any comparative
  genomics approach but are detectable by MAbsIslandScanner's compositional and structural
  evidence signals. The count here represents the lower bound — only islands with
  ≥{min_conf} confidence and ≥2 evidence lines are reported.
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.js"></script>
<script>
const D = {data_json};

new Chart(document.getElementById('c1'),{{
  type:'bar',
  data:{{
    labels:Object.keys(D.sys_types),
    datasets:[{{data:Object.values(D.sys_types),
      backgroundColor:'#0e7a5a',borderWidth:0,borderRadius:2}}]
  }},
  options:{{indexAxis:'y',responsive:true,maintainAspectRatio:false,
    plugins:{{legend:{{display:false}},
      tooltip:{{callbacks:{{label:ctx=>`${{ctx.parsed.x}} fixed defense islands`}}}}}},
    scales:{{x:{{grid:{{color:'rgba(128,128,128,0.1)'}},ticks:{{font:{{size:10}}}}}},
             y:{{grid:{{display:false}},ticks:{{font:{{size:10}}}}}}}}
  }}
}});

const AGE_COLORS={{very_recent:'#c0392b',recent:'#e67e22',moderate:'#3498db',old:'#bdc3c7'}};
const AGE_ORDER=['very_recent','recent','moderate','old'];
new Chart(document.getElementById('c2'),{{
  type:'bar',
  data:{{
    labels:AGE_ORDER.filter(k=>D.age_dist[k]).map(k=>k.replace('_',' ')),
    datasets:[{{data:AGE_ORDER.filter(k=>D.age_dist[k]).map(k=>D.age_dist[k]),
      backgroundColor:AGE_ORDER.filter(k=>D.age_dist[k]).map(k=>AGE_COLORS[k]),
      borderWidth:0,borderRadius:2}}]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    plugins:{{legend:{{display:false}},
      tooltip:{{callbacks:{{label:ctx=>`${{ctx.parsed.y}} defense systems`}}}}}},
    scales:{{x:{{grid:{{display:false}},ticks:{{font:{{size:11}}}}}},
             y:{{grid:{{color:'rgba(128,128,128,0.1)'}},ticks:{{font:{{size:10}}}}}}}}
  }}
}});

const confBadge=c=>{{
  const m={{high:'bh',moderate:'bm',low:'bl'}};
  return `<span class="badge ${{m[c]||'bl'}}">${{c}}</span>`;
}};
const ageBadge=a=>{{
  const m={{very_recent:'bvr',recent:'bm',moderate:'bl',old:'bl'}};
  return `<span class="badge ${{m[a]||'bl'}}">${{a?.replace('_',' ')||'—'}}</span>`;
}};

const tbl=document.getElementById('defTbl');
D.top_def.forEach(r=>{{
  const tr=document.createElement('tr');
  const len=r.island_length?Math.round(r.island_length/1000)+'kb':'—';
  tr.innerHTML=`
    <td style="font-size:11px;font-weight:500">${{r.strain}}</td>
    <td style="color:#0e7a5a">${{r.subtype}}</td>
    <td style="font-size:10px;color:#5a6380">${{r.island_start?r.island_start.toLocaleString():''}}–${{r.island_end?r.island_end.toLocaleString():''}} (${{len}})</td>
    <td>${{confBadge(r.island_confidence)}}</td>
    <td style="font-size:10px">${{r.island_evidence||'—'}}</td>
    <td>${{ageBadge(r.age_estimate)}}</td>
    <td style="font-family:monospace;font-size:11px">${{r.cai_ratio||'—'}}</td>
    <td style="font-size:11px">${{r.trna_flanked==='Yes'?r.trna_product||'Yes':'—'}}</td>
    <td style="font-size:11px">${{r.has_dr==='Yes'?r.dr_seq||'Yes':'—'}}</td>
  `;
  tbl.appendChild(tr);
}});
</script>
</body>
</html>"""

    with open(outpath,"w") as f:
        f.write(html)


if __name__ == "__main__":
    main()
