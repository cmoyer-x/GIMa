"""
build_patient_viewer.py

Patient longitudinal genome viewer with:
  - All patients (single and multiple isolates)
  - Pannable / zoomable genome track (scroll to zoom, drag to pan)
  - Minimap overview navigator
  - DefenseFinder, PADLOC, DefensePredictor tracks (toggleable)
  - Prophage regions track (from per-genome prediction CSVs)
  - RGP track
  - Defense system gained/lost badges vs first isolate

Usage:
    python build_patient_viewer.py \
        --json_dir genome_tracks \
        --intersection defense_rgp_intersection.tsv \
        --rgp_file abscessus_pangenome/rgp_output/regions_of_genomic_plasticity.tsv \
        --genome_stats abscessus_pangenome/genomes_statistics.tsv \
        --padloc padloc_all_systems.csv \
        --defense_predictor defense_gene_calls_tier1_2.csv \
        --prophage_dir prophage_csvs \
        --outfile patient_comparison_viewer.html
"""

import os, csv, json, re, argparse
from collections import defaultdict

SPOT_COLORS = {
    "spot_13": "#0e7a5a",
    "spot_7":  "#1a5fa8",
    "spot_48": "#a86000",
}
DEFENSE_COLORS = {
    "CBASS":"#0e7a5a","Hna":"#1a5fa8","Dnd":"#a86000",
    "RM":"#c0392b","RosmerTA":"#6c3db8","RloC":"#b5006e",
    "AbiAlpha":"#4a7c10","Wadjet":"#4a5568","Thoeris":"#0077aa",
    "BREX":"#8b4513","default":"#5a6380",
}

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json_dir",          required=True)
    p.add_argument("--intersection",      required=True)
    p.add_argument("--rgp_file",          required=True,
                   help="abscessus RGP TSV")
    p.add_argument("--rgp_file_mass",     default=None,
                   help="massiliense RGP TSV (optional)")
    p.add_argument("--spots_file_mass",   default=None,
                   help="massiliense spots TSV (optional)")
    p.add_argument("--intersection_mass", default=None,
                   help="massiliense defense intersection TSV (optional)")
    p.add_argument("--rgp_file_boll",     default=None,
                   help="bolletii RGP TSV (optional)")
    p.add_argument("--spots_file_boll",   default=None,
                   help="bolletii spots TSV (optional)")
    p.add_argument("--intersection_boll", default=None,
                   help="bolletii defense intersection TSV (optional)")
    p.add_argument("--island_catalog",      default=None,
                   help="island_viewer_data.json from build_island_catalog.py")
    p.add_argument("--scanner_islands",    default=None,
                   help="all_islands_combined.tsv from GIMa")
    p.add_argument("--denovo_defense",     default=None,
                   help="denovo_defense_intersection.tsv for fixed defense islands")
    p.add_argument("--trna_proximity_mass", default=None,
                   help="rgp_trna_proximity_massiliense.tsv (optional)")
    p.add_argument("--trna_proximity_boll", default=None,
                   help="rgp_trna_proximity_bolletii.tsv (optional)")
    p.add_argument("--min_evidence",       type=int, default=1,
                   help="Minimum HGT evidence lines to show RGP (0=all, 1=default, 2=stricter, 3=strongest)")
    p.add_argument("--trna_proximity",     default=None,
                   help="rgp_trna_proximity.tsv (optional)")
    p.add_argument("--genome_stats",      required=True)
    p.add_argument("--padloc",            default=None)
    p.add_argument("--defense_predictor", default=None)
    p.add_argument("--prophage_dir",      default=None,
                   help="Directory containing per-genome Depht HTML files (e.g. GD05.html)")
    p.add_argument("--outfile",           default="patient_comparison_viewer.html")
    return p.parse_args()


def canonical_strain(s):
    """Strip _WGS, _hybrid, _UNCUT suffixes for deduplication matching."""
    return re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=re.IGNORECASE)


def group_patients(strains):
    """
    Group ALL strains by patient ID (GD### prefix).
    Deduplicates strains that differ only by _WGS/_hybrid suffix
    (e.g. GD233A and GD233A_WGS are the same isolate).
    Prefers the suffixed version (has RGP/stats data) as the canonical key,
    keeping the shorter name only if no suffixed version exists.
    """
    # First pass: build canonical -> list of all name variants
    canonical_map = defaultdict(list)
    for s in strains:
        canonical_map[canonical_strain(s)].append(s)

    # Second pass: for each canonical group, pick one representative
    # Prefer _WGS > _hybrid > bare name
    deduped = {}
    for canon, variants in canonical_map.items():
        if len(variants) == 1:
            deduped[canon] = variants[0]
        else:
            # Prefer suffixed versions — they have more data attached
            wgs     = [v for v in variants if v.endswith('_WGS')]
            hybrid  = [v for v in variants if v.endswith('_hybrid')]
            if wgs:
                deduped[canon] = wgs[0]
            elif hybrid:
                deduped[canon] = hybrid[0]
            else:
                deduped[canon] = sorted(variants)[0]

    # Group deduplicated canonical strains by patient ID
    groups = defaultdict(list)
    for canon, representative in deduped.items():
        m = re.match(r'^(GD\d+)', canon)
        if m:
            groups[m.group(1)].append(representative)

    return {k: sorted(v) for k, v in groups.items()}


def load_genome_lengths(json_dir):
    lengths = {}
    for fname in os.listdir(json_dir):
        if not fname.endswith(".json"): continue
        strain = fname.replace(".json","")
        with open(os.path.join(json_dir, fname)) as f:
            data = json.load(f)
        lengths[strain] = data["genome_len"]
    print(f"  Loaded JSON tracks for {len(lengths)} strains")
    return lengths


def load_genome_stats(stats_file):
    stats = {}
    with open(stats_file) as f:
        header = None
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                col = {h: i for i, h in enumerate(header)}
                continue
            strain = parts[col["Genome_name"]]
            def g(key, default="NA"):
                try: return parts[col[key]]
                except: return default
            try:
                ctg = int(g("Contigs","999"))
                if ctg==1:      q="complete"
                elif ctg<=10:   q="high"
                elif ctg<=30:   q="moderate"
                elif ctg<=80:   q="fragmented"
                else:           q="highly_fragmented"
            except: q="unknown"
            stats[strain] = {
                "contigs":g("Contigs"),"genes":g("Genes"),
                "completeness":g("Completeness"),"contamination":g("Contamination"),
                "persistent":g("Persistent_families"),"shell":g("Shell_families"),
                "cloud":g("Cloud_families"),"frag_genes":g("Fragmented_genes"),
                "quality":q,
            }
    print(f"  Loaded assembly stats for {len(stats)} strains")
    return stats


def load_rgps(rgp_file):
    rgps = defaultdict(list)
    with open(rgp_file) as f:
        raw_header = f.readline().strip().split("\t")
    col = {h: i for i, h in enumerate(raw_header)}
    with open(rgp_file) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 8: continue
            genome = row[col["genome"]]
            rgps[genome].append({
                "rgp_id": row[col["region"]],
                "start": int(row[col["start"]]),
                "stop":  int(row[7]),
                "score": row[col["score"]] if "score" in col else "",
            })
    print(f"  Loaded RGPs for {len(rgps)} strains")
    return rgps


def load_defensefinder(intersection_file):
    df = defaultdict(list)
    with open(intersection_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["in_rgp"] not in ("Yes","No"): continue
            try:
                start = int(row["sys_start_coord"])
                stop  = int(row["sys_stop_coord"])
            except: continue
            df[row["strain"]].append({
                "sys_id":row["sys_id"],"type":row["type"],"subtype":row["subtype"],
                "start":start,"stop":stop,"in_rgp":row["in_rgp"],
                "rgp_id":row.get("rgp_id",""),"spot_id":row.get("spot_id",""),
                "tool":"DefenseFinder",
            })
    print(f"  Loaded DefenseFinder for {len(df)} strains")
    return df


def load_padloc(padloc_file):
    systems = defaultdict(lambda: {"starts":[],"stops":[],"proteins":[],"system":"","strain":""})
    with open(padloc_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            strain = row["strain"]
            key    = (strain, row["system"], row["system.number"])
            systems[key]["strain"]  = strain
            systems[key]["system"]  = row["system"]
            systems[key]["starts"].append(int(row["start"]))
            systems[key]["stops"].append(int(row["end"]))
            systems[key]["proteins"].append(row["protein.name"])
    pl = defaultdict(list)
    for (strain,system,snum), d in systems.items():
        pl[strain].append({
            "sys_id":f"{strain}_PADLOC_{system}_{snum}",
            "type":system,"subtype":system,
            "start":min(d["starts"]),"stop":max(d["stops"]),
            "proteins":", ".join(set(d["proteins"])),"tool":"PADLOC",
        })
    print(f"  Loaded PADLOC for {len(pl)} strains")
    return pl


def load_defense_predictor(dp_file):
    by_strain_cat = defaultdict(list)
    with open(dp_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                start = int(float(row["start"]))
                stop  = int(float(row["end"]))
            except: continue
            by_strain_cat[(row["strain"], row["defense_cat"])].append({
                "start":start,"stop":stop,"tier":row["confidence_tier"],
                "pfam":row["pfam_name"],"desc":row["defense_desc"],
                "prob":row["probability"],
            })
    dp = defaultdict(list)
    GAP = 10000
    for (strain, cat), genes in by_strain_cat.items():
        genes.sort(key=lambda g: g["start"])
        clusters, current = [], [genes[0]]
        for g in genes[1:]:
            if g["start"] - current[-1]["stop"] <= GAP: current.append(g)
            else: clusters.append(current); current = [g]
        clusters.append(current)
        for i, cluster in enumerate(clusters):
            probs = [float(g["prob"]) for g in cluster if g["prob"]]
            dp[strain].append({
                "sys_id":f"{strain}_DP_{cat}_{i+1}","type":cat,
                "subtype":cluster[0]["desc"],
                "start":min(g["start"] for g in cluster),
                "stop":max(g["stop"]  for g in cluster),
                "tiers":", ".join(set(g["tier"] for g in cluster)),
                "pfams":", ".join(set(g["pfam"] for g in cluster)),
                "mean_prob":round(sum(probs)/len(probs),1) if probs else 0,
                "n_genes":len(cluster),"tool":"DefensePredictor",
            })
    print(f"  Loaded DefensePredictor for {len(dp)} strains")
    return dp


def load_prophages(prophage_dir):
    """
    Load per-genome Depht prophage HTML files.
    Extracts genome-absolute coordinates from the summary table.
    HTML table format: Prophage Name | Left Coordinate | Right Coordinate | Length
    Falls back to CSV-based clustering if no HTML found.
    """
    import re

    prophages = defaultdict(list)
    n_files = 0

    for fname in os.listdir(prophage_dir):
        if not fname.endswith(".html"):
            continue
        strain = fname.replace(".html", "")
        fpath  = os.path.join(prophage_dir, fname)
        try:
            with open(fpath) as f:
                html = f.read()
            # Extract table rows: Prophage Name | Left | Right | Length
            rows = re.findall(r'<tr[^>]*>(.*?)</tr>', html, re.DOTALL)
            for row in rows:
                cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.DOTALL)
                clean = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]
                if len(clean) < 3:
                    continue
                # Skip header rows
                if clean[0].lower() in ('prophage name', 'name', 'prophage'):
                    continue
                try:
                    name  = clean[0]
                    start = int(clean[1])
                    stop  = int(clean[2])
                    length= int(clean[3]) if len(clean) > 3 else stop - start
                    # Confidence based on length — Depht high-confidence regions
                    # tend to be >20kb; use length as proxy since scores aren't in HTML
                    conf = "high" if length >= 20000 else "moderate"
                    prophages[strain].append({
                        "prophage_id": name,
                        "start":       start,
                        "stop":        stop,
                        "n_genes":     0,
                        "max_pred":    1.0,
                        "max_phage":   100.0,
                        "confidence":  conf,
                        "length":      length,
                    })
                except (ValueError, IndexError):
                    continue
            n_files += 1
        except Exception as e:
            print(f"  Warning: could not read {fname}: {e}")
            continue

    print(f"  Loaded Depht prophage HTML for {n_files} strains "
          f"({sum(len(v) for v in prophages.values())} regions total)")
    return prophages


def load_trna_proximity(trna_file):
    """
    Load tRNA proximity data.
    Returns dict: rgp_id -> {trna_flanked, min_trna_dist, trna_product,
                              hgt_evidence_count, hgt_evidence, gc_foreignness,
                              has_mobility, mobility_types}
    """
    trna_data = {}
    if not trna_file or not os.path.exists(trna_file):
        return trna_data
    with open(trna_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            prod = row.get("trna_left_product") or row.get("trna_right_product") or ""
            dist = row.get("min_trna_dist") or ""
            trna_data[row["rgp_id"]] = {
                "trna_flanked":      row["trna_flanked"],
                "min_trna_dist":     dist,
                "trna_product":      prod,
                "hgt_evidence_count":int(row.get("hgt_evidence_count") or 0),
                "hgt_evidence":      row.get("hgt_evidence",""),
                "gc_foreignness":    row.get("gc_foreignness",""),
                "has_mobility":      row.get("has_mobility","No"),
                "mobility_types":    row.get("mobility_types",""),
            }
    print(f"  Loaded tRNA proximity data for {len(trna_data)} RGPs")
    return trna_data


# ── Scanner island loaders ────────────────────────────────────────────────────
def load_island_catalog(path):
    """Load island_viewer_data.json from build_island_catalog.py.
    Returns dict: strain -> list of island dicts with nesting depth.
    """
    if not path or not os.path.exists(path):
        return {}
    import json as _json
    with open(path) as f:
        data = _json.load(f)
    total = sum(len(v) for v in data.values())
    print(f"  Island catalog: {total:,} islands across {len(data)} strains")
    return data


def load_scanner_islands(path):
    if not path or not os.path.exists(path):
        return {}
    from collections import defaultdict as dd
    islands_by_strain = dd(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                islands_by_strain[row["strain"]].append({
                    "start":       int(row["start"]),
                    "end":         int(row["end"]),
                    "length":      int(row["length"]),
                    "rgp_seed":    row.get("rgp_seed",""),
                    "confidence":  row.get("confidence",""),
                    "n_evidence":  int(row.get("n_evidence",0)),
                    "evidence":    row.get("evidence",""),
                    "age_estimate":row.get("age_estimate",""),
                    "cai_ratio":   row.get("cai_ratio",""),
                    "trna_flanked":row.get("trna_flanked","No"),
                    "has_dr":      row.get("has_dr","No"),
                    "mob_types":   row.get("mob_types",""),
                    "is_denovo":   not bool(row.get("rgp_seed","")),
                })
            except: continue
    total = sum(len(v) for v in islands_by_strain.values())
    print(f"  Scanner islands: {total:,} across {len(islands_by_strain)} strains")
    return dict(islands_by_strain)


def load_denovo_defense(path):
    if not path or not os.path.exists(path):
        return {}
    from collections import defaultdict as dd
    defense_by_strain = dd(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("in_denovo_island") != "Yes": continue
            try:
                defense_by_strain[row["strain"]].append({
                    "subtype":     row["subtype"],
                    "sys_start":   int(row["sys_start"])    if row["sys_start"]    else 0,
                    "sys_stop":    int(row["sys_stop"])     if row["sys_stop"]     else 0,
                    "island_start":int(row["island_start"]) if row["island_start"] else 0,
                    "island_end":  int(row["island_end"])   if row["island_end"]   else 0,
                    "confidence":  row.get("island_confidence",""),
                    "age_estimate":row.get("age_estimate",""),
                })
            except: continue
    total = sum(len(v) for v in defense_by_strain.values())
    print(f"  Fixed defense: {total:,} across {len(defense_by_strain)} strains")
    return dict(defense_by_strain)


def build_html(genome_lengths, rgps, df, pl, dp, prophages, stats, patient_groups, trna_data=None, min_evidence=1, scanner_islands=None, denovo_defense=None, island_catalog=None):
    # Build strain_data from UNION of all data sources
    # so strains with prophage/defense data but no JSON track are still included
    all_strains = set(genome_lengths.keys())
    all_strains |= set(rgps.keys())
    all_strains |= set(df.keys())
    all_strains |= set(pl.keys())
    all_strains |= set(dp.keys())
    all_strains |= set(prophages.keys())
    all_strains |= set(stats.keys())

    # Build canonical -> all variants map for data merging
    canon_variants = defaultdict(list)
    for s in all_strains:
        canon_variants[canonical_strain(s)].append(s)

    strain_data = {}
    for strain in all_strains:
        asm = stats.get(strain, {})
        fallback_len = 5000000
        try:
            genes = int(asm.get("genes", 0) or 0)
            if genes > 0:
                fallback_len = genes * 1000
        except: pass

        # Merge data from all name variants of this strain
        # e.g. GD233A_WGS gets defense data from GD233A too
        canon = canonical_strain(strain)
        variants = canon_variants.get(canon, [strain])

        merged_df      = []
        merged_padloc  = []
        merged_dp      = []
        merged_rgps    = rgps.get(strain, [])
        merged_prophage= prophages.get(strain, [])
        merged_asm     = asm

        for v in variants:
            merged_df     += df.get(v, [])
            merged_padloc += pl.get(v, [])
            merged_dp     += dp.get(v, [])
            if not merged_rgps:
                merged_rgps = rgps.get(v, [])
            if not merged_prophage:
                merged_prophage = prophages.get(v, [])
            if not merged_asm and stats.get(v):
                merged_asm = stats.get(v, {})

        # Deduplicate by sys_id
        seen = set()
        dedup_df = []
        for d in merged_df:
            if d["sys_id"] not in seen:
                seen.add(d["sys_id"]); dedup_df.append(d)

        seen = set()
        dedup_pl = []
        for d in merged_padloc:
            if d["sys_id"] not in seen:
                seen.add(d["sys_id"]); dedup_pl.append(d)

        seen = set()
        dedup_dp = []
        for d in merged_dp:
            if d["sys_id"] not in seen:
                seen.add(d["sys_id"]); dedup_dp.append(d)

        # Add tRNA proximity data to each RGP
        rgps_with_trna = []
        for rgp in merged_rgps:
            rgp_copy = dict(rgp)
            td = (trna_data or {}).get(rgp["rgp_id"], {})
            rgp_copy["trna_flanked"]       = td.get("trna_flanked", "No")
            rgp_copy["min_trna_dist"]      = td.get("min_trna_dist", "")
            rgp_copy["trna_product"]       = td.get("trna_product", "")
            rgp_copy["hgt_evidence_count"] = td.get("hgt_evidence_count", 0)
            rgp_copy["hgt_evidence"]       = td.get("hgt_evidence", "")
            rgp_copy["gc_foreignness"]     = td.get("gc_foreignness", "")
            rgp_copy["has_mobility"]       = td.get("has_mobility", "No")
            rgp_copy["mobility_types"]     = td.get("mobility_types", "")
            rgps_with_trna.append(rgp_copy)

        # Filter RGPs by minimum evidence level
        rgps_filtered = [r for r in rgps_with_trna
                         if r.get("hgt_evidence_count", 0) >= min_evidence
                         or r.get("has_defense_in_rgp", False)]

        # Scanner islands for this strain
        def _canon(s):
            import re as _re
            return _re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=_re.IGNORECASE)
        scan_isls = None
        if scanner_islands:
            scan_isls = scanner_islands.get(strain) or scanner_islands.get(_canon(strain))
            if not scan_isls:
                sc = _canon(strain)
                for k in scanner_islands:
                    if _canon(k) == sc:
                        scan_isls = scanner_islands[k]; break
        scan_isls = scan_isls or []
        fix_def = None
        if denovo_defense:
            fix_def = denovo_defense.get(strain) or denovo_defense.get(_canon(strain))
            if not fix_def:
                sc = _canon(strain)
                for k in denovo_defense:
                    if _canon(k) == sc:
                        fix_def = denovo_defense[k]; break
        fix_def = fix_def or []

        # Island catalog data (nested-resolved)
        cat_isls = None
        if island_catalog:
            cat_isls = island_catalog.get(strain) or island_catalog.get(_canon(strain))
            if not cat_isls:
                sc = _canon(strain)
                for k in island_catalog:
                    if _canon(k) == sc:
                        cat_isls = island_catalog[k]; break
        cat_isls = cat_isls or []

        # Build defense coordinate list for annotation
        # (used to flag which islands carry defense systems)
        all_defense = []
        for d in dedup_df:
            s = d.get("start")
            e = d.get("stop") or d.get("end")
            if s and e:
                all_defense.append({"start": s, "end": e})
        for d in fix_def:
            s = d.get("sys_start") or d.get("island_start", 0)
            e = d.get("sys_stop")  or d.get("island_end",   0)
            if s and e:
                all_defense.append({"start": s, "end": e})

        def has_defense_overlap(isl, defense_list):
            for d in defense_list:
                ds = d.get("start", 0)
                de = d.get("stop") or d.get("end") or 0
                if ds and de:
                    if isl["start"] <= de and isl["end"] >= ds:
                        return True
            return False

        # Show ALL catalog islands — flag which ones carry defense systems
        # Genomic islands are defined by HGT origin, not cargo type
        defense_cat_isls = []
        for isl in cat_isls:
            isl_copy = dict(isl)
            isl_copy["has_defense"] = has_defense_overlap(isl, all_defense)
            defense_cat_isls.append(isl_copy)

        sdata = {
            "length":         genome_lengths.get(strain) or genome_lengths.get(_canon(strain)) or fallback_len,
            "rgps":           rgps_filtered,
            "scanner_islands":scan_isls,
            "fixed_defense":  fix_def,
            "catalog_islands":defense_cat_isls,
            "df":       dedup_df,
            "padloc":   dedup_pl,
            "dp":       dedup_dp,
            "prophage": merged_prophage,
            "asm":      merged_asm,
        }
        strain_data[strain] = sdata

    # All patients — single and multi isolate
    patient_index = {}
    for pid, isos in sorted(patient_groups.items(),
                             key=lambda x: (x[0].replace("GD","").zfill(6))):
        valid = [s for s in isos if s in strain_data]
        if valid:
            patient_index[pid] = valid

    patient_json        = json.dumps(patient_index)
    # strain_data is written to per-strain JSON files — not embedded in HTML
    spot_colors_json    = json.dumps(SPOT_COLORS)
    defense_colors_json = json.dumps(DEFENSE_COLORS)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Patient longitudinal defense viewer</title>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap" rel="stylesheet">
<style>
  :root {{
    --bg:#f4f6fa;--surface:#fff;--border:#dde1ea;
    --text:#1a1f2e;--muted:#5a6380;--accent:#0e7a5a;
    --font-mono:'IBM Plex Mono',monospace;--font-sans:'IBM Plex Sans',sans-serif;
  }}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{background:var(--bg);color:var(--text);font-family:var(--font-sans);font-size:14px}}
  header{{border-bottom:1px solid var(--border);padding:1.25rem 2rem;
    display:flex;align-items:center;gap:1.5rem;background:var(--surface)}}
  header h1{{font-size:14px;font-weight:500;letter-spacing:.08em;text-transform:uppercase;
    color:var(--accent);font-family:var(--font-mono)}}
  header span{{font-size:12px;color:var(--muted);font-family:var(--font-mono)}}
  .controls{{padding:1rem 2rem;display:flex;align-items:center;gap:1rem;
    border-bottom:1px solid var(--border);flex-wrap:wrap;background:var(--surface)}}
  .controls label{{font-size:11px;text-transform:uppercase;letter-spacing:.08em;
    color:var(--muted);font-family:var(--font-mono)}}
  select{{background:var(--bg);border:1px solid var(--border);color:var(--text);
    padding:6px 12px;border-radius:4px;font-family:var(--font-mono);
    font-size:13px;cursor:pointer;min-width:200px}}
  select:focus{{outline:none;border-color:var(--accent)}}
  .toggles{{display:flex;gap:6px;align-items:center;flex-wrap:wrap}}
  .toggle-btn{{display:flex;align-items:center;gap:5px;padding:4px 10px;border-radius:20px;
    border:1.5px solid;font-family:var(--font-mono);font-size:11px;cursor:pointer;
    font-weight:500;transition:opacity .15s;user-select:none}}
  .toggle-btn.df    {{border-color:#0e7a5a;color:#0e7a5a;background:#f0faf6}}
  .toggle-btn.padloc{{border-color:#1a5fa8;color:#1a5fa8;background:#f0f5ff}}
  .toggle-btn.dp    {{border-color:#8b4513;color:#8b4513;background:#fff8f0}}
  .toggle-btn.rgp   {{border-color:#888780;color:#5a6380;background:#f5f5f5}}
  .toggle-btn.pp    {{border-color:#7c2d92;color:#7c2d92;background:#faf0ff}}
  .toggle-btn.off   {{opacity:.3}}
  .toggle-dot{{width:7px;height:7px;border-radius:50%;flex-shrink:0}}
  .df .toggle-dot{{background:#0e7a5a}}.padloc .toggle-dot{{background:#1a5fa8}}
  .dp .toggle-dot{{background:#8b4513}}.rgp .toggle-dot{{background:#888780}}
  .pp .toggle-dot{{background:#7c2d92}}
  .patient-summary{{padding:.6rem 2rem;border-bottom:1px solid var(--border);
    background:#f8fafd;display:flex;gap:2rem;align-items:center;flex-wrap:wrap}}
  .ps-item{{font-family:var(--font-mono);font-size:12px;color:var(--muted)}}
  .ps-item span{{color:var(--text);font-weight:500}}
  .change-badge{{display:inline-block;font-size:10px;font-family:var(--font-mono);
    padding:1px 7px;border-radius:3px;font-weight:600;margin:1px}}
  .gained{{background:#d4f2e8;color:#0e7a5a}}
  .lost{{background:#fde8e4;color:#a82c1a}}
  .stable{{background:#eaecf0;color:#5a6380}}
  .isolate-list{{padding:1.5rem 2rem;display:flex;flex-direction:column;gap:1.5rem}}
  .isolate-block{{background:var(--surface);border:1px solid var(--border);border-radius:10px;overflow:hidden}}
  .isolate-header{{display:flex;align-items:center;gap:.75rem;padding:.6rem 1.25rem;
    border-bottom:1px solid var(--border);background:#f8fafd;flex-wrap:wrap}}
  .isolate-name{{font-family:var(--font-mono);font-size:14px;font-weight:500;
    color:var(--text);min-width:110px}}
  .iso-stats{{display:flex;gap:.75rem;flex-wrap:wrap}}
  .iso-stat{{font-family:var(--font-mono);font-size:11px;color:var(--muted)}}
  .iso-stat span{{color:var(--text);font-weight:500}}
  .asm-pill{{font-size:10px;font-family:var(--font-mono);padding:2px 8px;border-radius:10px;font-weight:500}}
  .pill-complete,.pill-high{{background:#d4f2e8;color:#0e7a5a}}
  .pill-moderate{{background:#fff3cc;color:#8a5c00}}
  .pill-fragmented{{background:#fde8e4;color:#a82c1a}}
  .pill-highly_fragmented{{background:#fbd1cb;color:#a82c1a}}
  .pill-unknown{{background:#eaecf0;color:#5a6380}}
  /* Zoom controls */
  .zoom-controls{{display:flex;gap:4px;align-items:center;margin-left:auto}}
  .zoom-btn{{background:var(--surface);border:1px solid var(--border);color:var(--text);
    width:26px;height:26px;border-radius:4px;cursor:pointer;font-size:14px;
    display:flex;align-items:center;justify-content:center;font-family:var(--font-mono);
    user-select:none}}
  .zoom-btn:hover{{background:var(--bg);border-color:var(--accent)}}
  .zoom-label{{font-family:var(--font-mono);font-size:11px;color:var(--muted);min-width:38px;text-align:center}}
  /* Track canvas wrapper */
  .track-wrap{{padding:.75rem 1.25rem 1rem;position:relative}}
  canvas{{display:block;cursor:grab}}
  canvas.panning{{cursor:grabbing}}
  /* Minimap */
  .minimap-wrap{{padding:0 1.25rem .75rem;}}
  canvas.minimap{{cursor:pointer;border-radius:3px;border:1px solid var(--border)}}
  .legend{{display:flex;flex-wrap:wrap;gap:.75rem;padding:.75rem 2rem;
    border-top:1px solid var(--border);background:var(--surface)}}
  .legend-item{{display:flex;align-items:center;gap:5px;font-size:11px;
    color:var(--muted);font-family:var(--font-mono)}}
  .legend-swatch{{width:11px;height:11px;border-radius:2px;flex-shrink:0}}
  .legend-divider{{width:1px;height:16px;background:var(--border);margin:0 2px}}
  .tooltip{{position:fixed;background:#fff;border:1px solid #c5cad6;border-radius:6px;
    padding:10px 14px;font-family:var(--font-mono);font-size:12px;
    pointer-events:none;display:none;z-index:1000;max-width:300px;
    line-height:1.8;box-shadow:0 4px 16px rgba(0,0,0,0.12)}}
  .tt-title{{font-weight:500;color:#1a1f2e;margin-bottom:4px;font-size:13px}}
  .tt-badge{{display:inline-block;font-size:10px;padding:1px 6px;border-radius:3px;
    font-weight:600;margin-left:6px}}
  .tt-row{{color:#5a6380}}.tt-row span{{color:#1a1f2e}}
  .tt-df{{background:#d4f2e8;color:#0e7a5a}}
  .tt-padloc{{background:#ddeeff;color:#1a5fa8}}
  .tt-dp{{background:#fff0e0;color:#8b4513}}
  .tt-pp{{background:#f3e8ff;color:#7c2d92}}
  .no-data{{padding:2rem;text-align:center;color:var(--muted);font-family:var(--font-mono)}}
</style>
</head>
<body>

<header>
  <h1>Patient longitudinal viewer</h1>
  <span>M. abscessus — serial isolate defense &amp; prophage comparison</span>
</header>

<div class="controls">
  <label>Patient</label>
  <select id="patientSelect"></select>
  <div class="toggles">
    <span style="font-size:11px;color:var(--muted);font-family:var(--font-mono);
      text-transform:uppercase;letter-spacing:.07em;white-space:nowrap">Tracks:</span>
    <div class="toggle-btn df"     id="tog-df"     onclick="toggleTrack('df')"><div class="toggle-dot"></div>DefenseFinder</div>
    <div class="toggle-btn padloc" id="tog-padloc"  onclick="toggleTrack('padloc')"><div class="toggle-dot"></div>PADLOC</div>
    <div class="toggle-btn dp"     id="tog-dp"     onclick="toggleTrack('dp')"><div class="toggle-dot"></div>DefensePredictor</div>
    <div class="toggle-btn pp"     id="tog-pp"     onclick="toggleTrack('pp')"><div class="toggle-dot"></div>Prophage</div>
    <div class="toggle-btn rgp"    id="tog-rgp"    onclick="toggleTrack('rgp')"><div class="toggle-dot"></div>RGPs</div>
    <div class="toggle-btn"        id="tog-fix"    onclick="toggleTrack('fix')"
         style="background:#c0392b22;border-color:#c0392b"><div class="toggle-dot" style="background:#c0392b"></div>Fixed defense</div>
    <div class="toggle-btn"        id="tog-isl"    onclick="toggleTrack('isl')"
         style="background:#1D9E7522;border-color:#1D9E75"><div class="toggle-dot" style="background:#1D9E75"></div>Genomic islands</div>
  </div>
  <div style="display:flex;align-items:center;gap:8px;margin-left:8px">
    <span style="font-size:11px;color:var(--muted);font-family:var(--font-mono);text-transform:uppercase;letter-spacing:.07em;white-space:nowrap">RGP evidence:</span>
    <select id="evidenceFilter" onchange="MIN_EVIDENCE=parseInt(this.value);const pid=document.getElementById('patientSelect').value;if(pid)redrawAll(pid);"
      style="font-size:11px;padding:4px 8px;min-width:0">
      <option value="0">All RGPs</option>
      <option value="1" selected>1+ lines</option>
      <option value="2">2+ lines</option>
      <option value="3">3+ lines (strong)</option>
      <option value="4">4 lines (triple+defense)</option>
    </select>
  </div>
</div>

<div class="patient-summary" id="patientSummary"></div>
<div class="isolate-list"    id="isolateList"></div>

<div class="legend">
  <div class="legend-item"><div class="legend-swatch" style="background:#0e7a5a"></div>spot_13 CBASS_II</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#1a5fa8"></div>spot_7 Hna</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#a86000"></div>spot_48 Dnd</div>
  <div class="legend-divider"></div>
  <div class="legend-item"><div class="legend-swatch" style="background:#0e7a5a;opacity:.5"></div>DF</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#1a5fa8"></div>PADLOC</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#8b4513"></div>DefPredictor</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#7c2d92"></div>Prophage (high)</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#b57cc4"></div>Prophage (mod)</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#95a5a6"></div>RGP — 1 evidence line</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#2980b9"></div>RGP — tRNA flanked</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#e67e22"></div>RGP — 2 evidence lines</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#c0392b"></div>RGP — 3–4 lines (strongest)</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#0e7a5a"></div>GI — defense cargo</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#7c2d92"></div>GI — mobility cargo</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#c0392b"></div>GI — TA system cargo</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#378ADD"></div>GI — other cargo</div>
  <div class="legend-item" style="font-size:10px;color:var(--muted)">Dashed border = nested child island · GI+1 label = nesting depth</div>
  <div class="legend-divider"></div>
  <div class="legend-item"><div class="legend-swatch" style="background:#c0392b"></div>Fixed defense — very recent</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#e67e22"></div>Fixed defense — recent</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#3498db"></div>Fixed defense — moderate/old</div>
  <div class="legend-divider"></div>
  <div class="legend-item" style="font-size:11px;color:var(--muted);font-family:var(--font-mono)">HGT lines: GC foreign · tRNA flanked · mobility gene · defense system</div>
  <div class="legend-divider"></div>
  <div class="legend-item"><span class="change-badge gained">+gained</span></div>
  <div class="legend-item"><span class="change-badge lost">-lost</span></div>
  <div class="legend-item"><span class="change-badge stable">stable</span></div>
</div>

<div class="tooltip" id="tooltip"></div>

<script>
const STRAIN_DATA      = {{}};  // populated on demand via fetch
const PATIENT_GROUPS   = {patient_json};
const SPOT_COLORS      = {spot_colors_json};
const DEFENSE_COLORS   = {defense_colors_json};
const DATA_DIR         = 'viewer_data';  // relative path to per-strain JSON files
let   _fetchCache      = {{}};           // cache fetched strain data
const KEY_SPOTS        = new Set(['spot_13','spot_7','spot_48']);
const visible          = {{df:true,padloc:true,dp:true,pp:true,rgp:true,fix:true,isl:true}};
let   MIN_EVIDENCE     = {min_evidence};  // filter RGPs by HGT evidence count

// Per-isolate view state (zoom/pan per canvas)
const viewState = {{}};
const canvasHits = {{}};

function toggleTrack(track) {{
  visible[track] = !visible[track];
  document.getElementById('tog-'+track).classList.toggle('off', !visible[track]);
  const pid = document.getElementById('patientSelect').value;
  if (pid) redrawAll(pid);
}}

// ── Colors ───────────────────────────────────────────────────────────────────
function defColor(d) {{
  if (d.tool==='DefenseFinder') {{
    if (d.in_rgp==='Yes' && KEY_SPOTS.has(d.spot_id)) return SPOT_COLORS[d.spot_id];
    return DEFENSE_COLORS[d.type] || DEFENSE_COLORS.default;
  }}
  if (d.tool==='PADLOC')           return '#1a5fa8';
  if (d.tool==='DefensePredictor') return '#8b4513';
  return '#5a6380';
}}
function ppColor(pp) {{
  if (pp.confidence==='high')     return '#7c2d92';
  if (pp.confidence==='moderate') return '#b57cc4';
  return '#d4b8e0';
}}
function rgpColor(rgp, defInRGP) {{
  // Returns [fillColor, strokeColor] as solid colors like prophage track
  const s = defInRGP.filter(d=>d.rgp_id===rgp.rgp_id&&KEY_SPOTS.has(d.spot_id)).map(d=>d.spot_id);
  if (s.length>0) {{
    const c = SPOT_COLORS[s[0]];
    return [c, c];
  }}
  const ev = rgp.hgt_evidence_count || 0;
  if (ev >= 3) return ['#c0392b', '#922b21'];   // red — 3+ evidence
  if (ev >= 2) return ['#e67e22', '#b7600d'];   // amber — 2 evidence
  if (rgp.trna_flanked === 'Yes') return ['#2980b9','#1a5fa8']; // blue — tRNA flanked
  return ['#95a5a6', '#7f8c8d'];                // gray — 1 evidence
}}

function rgpFill(rgp, defInRGP) {{
  return rgpColor(rgp, defInRGP)[0];
}}

// ── Draw a single track canvas ────────────────────────────────────────────────
function drawTrack(canvasEl, minimapEl, strain, maxLen) {{
  const data = STRAIN_DATA[strain];
  if (!data) return;

  const gLen  = data.length;
  const rgps         = (data.rgps || []).filter(r =>
    (r.hgt_evidence_count || 0) >= MIN_EVIDENCE ||
    (data.df||[]).some(d => d.rgp_id === r.rgp_id)
  );
  const scannerIslands = data.scanner_islands || [];
  const fixedDefense   = data.fixed_defense   || [];
  const catalogIslands = data.catalog_islands || [];
  const df    = data.df     || [];
  const pl    = data.padloc || [];
  const dp    = data.dp     || [];
  const pp    = data.prophage|| [];
  const defInRGP = df.filter(d=>d.in_rgp==='Yes');

  const id    = canvasEl.id;
  if (!viewState[id]) viewState[id] = {{zoom:1, pan:0}};
  const vs    = viewState[id];

  const DPR   = window.devicePixelRatio || 1;
  const W     = canvasEl.parentElement.clientWidth - 4;
  const PAD   = 56;
  const ARROW_H=18, GENOME_H=12, RGP_H=14, PP_H=16, TRACK_GAP=8;

  // Estimate height
  const dfR  = visible.df     ? Math.min(df.length,5)  : 0;
  const plR  = visible.padloc ? Math.min(pl.length,5)  : 0;
  const dpR  = visible.dp     ? Math.min(dp.length,5)  : 0;
  const ppR  = visible.pp     ? 1 : 0;
  const aboveH = (dfR+plR+dpR)*( ARROW_H+3) + (dfR>0?TRACK_GAP:0)+(plR>0?TRACK_GAP:0)+(dpR>0?TRACK_GAP:0) + 24;
  // Dynamic height: estimate lanes needed from island count and genome length
  const catIslands  = data.catalog_islands || [];
  const ISL_ROW_H   = 12;
  const ISL_GAP     = 3;
  // Estimate max lanes: islands / (genome_length / avg_island_width_px)
  // Conservative: assume up to 4 lanes for dense genomes, 1 for sparse
  const islDensity  = catIslands.length / Math.max(1, (data.length||5000000) / 50000);
  const estLanes    = visible.isl && catIslands.length > 0
    ? Math.min(5, Math.max(1, Math.ceil(islDensity)))
    : 0;
  // Reserve minimum height for island track even when lane count is unknown
  const islTrackH   = visible.isl
    ? Math.max(4 * (ISL_ROW_H + ISL_GAP) + 8, estLanes * (ISL_ROW_H + ISL_GAP) + (estLanes > 0 ? 8 : 0))
    : 0;
  const belowH = (visible.rgp?RGP_H+6:0) + (visible.pp?PP_H+6:0) +
                 (visible.fix?8+4:0) + islTrackH + 24;
  const H = Math.max(120, aboveH + GENOME_H + belowH);

  canvasEl.width  = W*DPR; canvasEl.height = H*DPR;
  canvasEl.style.width=W+'px'; canvasEl.style.height=H+'px';
  // Ensure container expands to fit canvas on retina displays
  canvasEl.style.minHeight=H+'px';

  const ctx = canvasEl.getContext('2d');
  ctx.scale(DPR, DPR);
  ctx.clearRect(0,0,W,H);

  // Viewport: pan is in bp, zoom scales bp→px
  const trackW = W - PAD - 8;
  const visibleBp = maxLen / vs.zoom;
  const bpPerPx   = visibleBp / trackW;
  const panClamped= Math.max(0, Math.min(vs.pan, maxLen - visibleBp));
  vs.pan = panClamped;

  function bpToX(bp) {{ return PAD + (bp - panClamped) / bpPerPx; }}
  function bpWidth(bp) {{ return bp / bpPerPx; }}
  function inView(start, stop) {{ return stop >= panClamped && start <= panClamped + visibleBp; }}

  const hits = [];
  const GENOME_Y = aboveH + 8;

  // Reference line
  ctx.fillStyle='#e8eaf0';
  ctx.fillRect(PAD, GENOME_Y-1, trackW, 1);

  // Axis ticks
  ctx.font='9px IBM Plex Mono,monospace';
  const tickCount = 7;
  for (let i=0;i<=tickCount;i++) {{
    const bp = Math.round((panClamped + visibleBp*i/tickCount));
    const x  = PAD + trackW*i/tickCount;
    ctx.fillStyle='#c5cad6'; ctx.fillRect(x, GENOME_Y+GENOME_H, 1, 4);
    ctx.fillStyle='#5a6380';
    ctx.textAlign = i===0?'left':i===tickCount?'right':'center';
    const lbl=bp>=1e6?(bp/1e6).toFixed(2)+'M':bp>=1e3?(bp/1e3).toFixed(0)+'k':bp;
    ctx.fillText(lbl, x, GENOME_Y+GENOME_H+13);
  }}

  // Genome backbone — only show actual genome length
  const genomeEndX = Math.min(bpToX(gLen), PAD+trackW);
  const grad=ctx.createLinearGradient(PAD,0,genomeEndX,0);
  grad.addColorStop(0,'#d0d5e0');grad.addColorStop(0.5,'#bcc2d0');grad.addColorStop(1,'#d0d5e0');
  ctx.fillStyle=grad; ctx.strokeStyle='#a0a8bc'; ctx.lineWidth=1;
  ctx.beginPath();
  ctx.roundRect(PAD, GENOME_Y, Math.max(0,genomeEndX-PAD), GENOME_H, 3);
  ctx.fill(); ctx.stroke();
  hits.push({{type:'genome',data:{{length:gLen,strain}},
    x:PAD,y:GENOME_Y,w:genomeEndX-PAD,h:GENOME_H,bpPerPx,panClamped}});

  // RGP track
  if (visible.rgp) {{
    const RGP_Y = GENOME_Y+GENOME_H+4;
    rgps.forEach(rgp => {{
      if (!inView(rgp.start,rgp.stop)) return;
      const x=bpToX(rgp.start), w=Math.max(2,bpWidth(rgp.stop-rgp.start));
      const [fill,stroke] = rgpColor(rgp,defInRGP);
      ctx.fillStyle   = fill + 'cc';
      ctx.strokeStyle = stroke;
      ctx.lineWidth=1;
      ctx.beginPath();ctx.roundRect(x,RGP_Y,w,RGP_H,3);ctx.fill();ctx.stroke();
      hits.push({{type:'rgp',data:rgp,x,y:RGP_Y,w,h:RGP_H}});
    }});
    ctx.fillStyle='#5a6380';ctx.font='9px IBM Plex Mono,monospace';
    ctx.textAlign='right';ctx.fillText('RGPs',PAD-4,RGP_Y+RGP_H/2+3);
  }}

  // Prophage track
  if (visible.pp && pp.length>0) {{
    const PP_Y = GENOME_Y+GENOME_H+(visible.rgp?RGP_H+8:4);
    pp.forEach(p => {{
      if (!inView(p.start,p.stop)) return;
      const x=bpToX(p.start), w=Math.max(3,bpWidth(p.stop-p.start));
      const c=ppColor(p);
      ctx.fillStyle=c+'bb'; ctx.strokeStyle=c; ctx.lineWidth=1;
      ctx.beginPath();ctx.roundRect(x,PP_Y,w,PP_H,3);ctx.fill();ctx.stroke();
      if (w>40) {{
        ctx.fillStyle='#fff';ctx.font='8px IBM Plex Mono,monospace';ctx.textAlign='left';
        ctx.fillText(p.confidence+' conf',x+4,PP_Y+PP_H/2+3);
      }}
      hits.push({{type:'prophage',data:p,x,y:PP_Y,w,h:PP_H}});
    }});

    // ── Fixed defense islands track ────────────────────────────────────────
    const FIX_Y = PP_Y + PP_H + 3;
    const FIX_H = 7;
    if (visible.fix) {{
      fixedDefense.forEach(fd => {{
        const s = fd.island_start || fd.sys_start || 0;
        const e = fd.island_end   || fd.sys_stop  || 0;
        if (!s || !e || !inView(s,e)) return;
        const x = bpToX(s), w = Math.max(3, bpWidth(e-s));
        const age = fd.age_estimate || '';
        const col = age==='very_recent'?'#c0392b':age==='recent'?'#e67e22':
                    age==='moderate'?'#3498db':'#888780';
        ctx.fillStyle=col+'dd'; ctx.strokeStyle=col; ctx.lineWidth=1;
        ctx.beginPath(); ctx.roundRect(x,FIX_Y,w,FIX_H,2); ctx.fill(); ctx.stroke();
        hits.push({{type:'fixed_defense',data:fd,x,y:FIX_Y,w,h:FIX_H}});
      }});
    }}

    // ── Nested genomic island catalog tracks ─────────────────────────────────
    // Each nesting depth gets its own row. Top-level (depth 0) are full height.
    // Nested children are shorter and offset downward, visually inside the parent.
    const ISL_ROW_H  = 12;
    const ISL_GAP    = 3;
    const ISL_BASE_Y = FIX_Y + (visible.fix ? FIX_H + 6 : 0);

    function islColor(isl) {{
      // Color by dominant cargo — all islands shown
      const cargo = isl.dominant_cargo || '';
      const hasDef = isl.has_defense;
      if (hasDef)              return ['#0e7a5a','#085041'];
      if (cargo==='mobility')  return ['#7c2d92','#5b1f6e'];
      if (cargo==='ta_system') return ['#c0392b','#922b21'];
      if (cargo==='metal')     return ['#e67e22','#b7600d'];
      if (cargo==='efflux')    return ['#2980b9','#1a5fa8'];
      if (cargo==='phage')     return ['#8e44ad','#6c3483'];
      if (cargo==='regulatory')return ['#16a085','#0d7566'];
      if (cargo==='defense')   return ['#0e7a5a','#085041'];
      // Unknown/hypothetical
      const ev = isl.n_evidence || 0;
      if (ev>=4) return ['#1D9E75','#0F6E56'];
      if (ev>=3) return ['#378ADD','#185FA5'];
      return ['#888780','#5F5E5A'];
    }}

    if (visible.isl) {{
      const catIslands = (data.catalog_islands || [])
        .filter(i => inView(i.start, i.end))
        .sort((a,b) => a.start - b.start);

      if (catIslands.length > 0) {{

        // ── Lane assignment ─────────────────────────────────────────────────
        // Assign each island to the lowest lane where it doesn't overlap
        // any previously placed island. This prevents visual overlap.
        const laneEnds = [];   // laneEnds[lane] = pixel x where last island ends
        const LANE_PAD = 2;    // minimum pixel gap between islands in same lane
        const islLanes = catIslands.map(isl => {{
          const x = bpToX(isl.start);
          const w = Math.max(4, bpWidth(isl.end - isl.start));
          // Find lowest lane with room
          let lane = 0;
          while (laneEnds[lane] !== undefined && laneEnds[lane] > x - LANE_PAD) {{
            lane++;
          }}
          laneEnds[lane] = x + w;
          return lane;
        }});

        const maxLane = Math.max(...islLanes);

        // Draw islands in their assigned lanes
        catIslands.forEach((isl, idx) => {{
          const lane  = islLanes[idx];
          const x     = bpToX(isl.start);
          const w     = Math.max(4, bpWidth(isl.end - isl.start));
          const h     = Math.max(6, ISL_ROW_H - lane);   // slightly shorter in deeper lanes
          const y     = ISL_BASE_Y + lane * (ISL_ROW_H + ISL_GAP);
          const [fill, stroke] = islColor(isl);

          // Nested children get dashed border; top-level solid
          const isNested = (isl.status === 'nested_child') || lane > 0;
          ctx.fillStyle   = fill + (isNested ? '99' : 'cc');
          ctx.strokeStyle = stroke;
          ctx.lineWidth   = 1;
          ctx.setLineDash(isNested ? [2,2] : []);
          ctx.beginPath();
          const r2 = Math.min(2, w/2, h/2);
          ctx.moveTo(x+r2, y);
          ctx.lineTo(x+w-r2, y);
          ctx.arcTo(x+w, y, x+w, y+r2, r2);
          ctx.lineTo(x+w, y+h-r2);
          ctx.arcTo(x+w, y+h, x+w-r2, y+h, r2);
          ctx.lineTo(x+r2, y+h);
          ctx.arcTo(x, y+h, x, y+h-r2, r2);
          ctx.lineTo(x, y+r2);
          ctx.arcTo(x, y, x+r2, y, r2);
          ctx.closePath();
          ctx.fill(); ctx.stroke();
          ctx.setLineDash([]);

          // BLAST-validated badge — small star on top-right corner
          if (isl.blast_validated) {{
            ctx.fillStyle = '#f59e0b';
            ctx.font = 'bold 8px sans-serif';
            ctx.textAlign = 'right';
            ctx.textBaseline = 'top';
            ctx.fillText('★', x+w-1, y+1);
            ctx.textBaseline = 'middle';
          }}

          // Label wide islands
          if (w > 80) {{
            ctx.fillStyle = '#fff';
            ctx.font = `bold 8px IBM Plex Mono,monospace`;
            ctx.textBaseline = 'middle';
            ctx.textAlign = 'left';
            // Show island_id if wide enough, otherwise cargo type
            const label = w > 150
              ? (isl.island_id || isl.dominant_cargo || '').substring(0,12)
              : (isl.dominant_cargo || '').substring(0,8);
            ctx.fillText(label, x+3, y+h/2);
          }}

          hits.push({{type:'catalog_island',data:isl,x,y,w,h}});
        }});

        // Lane labels on left margin
        for (let lane=0; lane<=maxLane; lane++) {{
          const labelY = ISL_BASE_Y + lane*(ISL_ROW_H+ISL_GAP) + ISL_ROW_H/2;
          ctx.fillStyle = 'rgba(128,128,128,0.5)';
          ctx.font = '7px IBM Plex Mono,monospace';
          ctx.textBaseline = 'middle';
          ctx.textAlign = 'right';
          ctx.fillText(lane===0 ? 'GI' : 'GI+'+lane, PAD-2, labelY);
          ctx.textAlign = 'left';
        }}
      }}
    }}
    ctx.fillStyle='#7c2d92';ctx.font='9px IBM Plex Mono,monospace';
    ctx.textAlign='right';
    const PP_Y2 = GENOME_Y+GENOME_H+(visible.rgp?RGP_H+8:4);
    ctx.fillText('Phage',PAD-4,PP_Y2+PP_H/2+3);
  }}

  // Arrow helper — draws above genome
  function drawArrows(systems, baseY, colorFn, trackLabel) {{
    if (!systems.length) return 0;
    const levels={{}};
    const sorted=[...systems].filter(d=>inView(d.start,d.stop)).sort((a,b)=>a.start-b.start);
    if (!sorted.length) return 0;
    sorted.forEach(d => {{
      const x=bpToX(d.start), w=Math.max(4,bpWidth(d.stop-d.start));
      let lv=0;
      while(true){{const occ=levels[lv]||0;if(x>=occ){{levels[lv]=x+w+2;break;}}lv++;}}
      const y=baseY-lv*(ARROW_H+3);
      const c=colorFn(d);
      ctx.fillStyle=c+'cc';ctx.strokeStyle=c;ctx.lineWidth=1;
      const aw=Math.min(7,w*0.3);
      ctx.beginPath();
      if(w>10){{ctx.moveTo(x,y);ctx.lineTo(x+w-aw,y);ctx.lineTo(x+w,y+ARROW_H/2);
        ctx.lineTo(x+w-aw,y+ARROW_H);ctx.lineTo(x,y+ARROW_H);}}
      else ctx.rect(x,y,w,ARROW_H);
      ctx.closePath();ctx.fill();ctx.stroke();
      if(w>28){{ctx.fillStyle='#fff';ctx.font='8px IBM Plex Mono,monospace';ctx.textAlign='left';
        ctx.fillText((d.subtype||d.type||'').replace(/_/g,' ').substring(0,14),x+3,y+ARROW_H/2+3);}}
      ctx.strokeStyle=c+'44';ctx.lineWidth=1;ctx.setLineDash([2,3]);
      ctx.beginPath();ctx.moveTo(x+w/2,y+ARROW_H);ctx.lineTo(x+w/2,GENOME_Y);
      ctx.stroke();ctx.setLineDash([]);
      hits.push({{type:'defense',data:d,x,y,w,h:ARROW_H}});
    }});
    const rows=Object.keys(levels).length;
    ctx.fillStyle='#5a6380';ctx.font='9px IBM Plex Mono,monospace';ctx.textAlign='right';
    ctx.fillText(trackLabel,PAD-4,baseY-2);
    return rows;
  }}

  let curY=GENOME_Y-TRACK_GAP;
  if(visible.df && df.length>0) {{
    const rows=drawArrows(df,curY-ARROW_H,defColor,'DF');
    curY-=Math.min(rows||1,5)*(ARROW_H+3)+TRACK_GAP;
  }}
  if(visible.padloc && pl.length>0) {{
    const rows=drawArrows(pl,curY-ARROW_H,()=>'#1a5fa8','PL');
    curY-=Math.min(rows||1,5)*(ARROW_H+3)+TRACK_GAP;
  }}
  if(visible.dp && dp.length>0) {{
    drawArrows(dp,curY-ARROW_H,()=>'#8b4513','DP');
  }}

  canvasHits[id] = hits;

  // Zoom label
  const zl = document.getElementById('zl_'+id);
  if (zl) zl.textContent = vs.zoom===1?'1×':vs.zoom.toFixed(1)+'×';

  // Draw minimap
  if (minimapEl) drawMinimap(minimapEl, strain, gLen, maxLen, panClamped, visibleBp);
}}

// ── Minimap ──────────────────────────────────────────────────────────────────
function drawMinimap(mmEl, strain, gLen, maxLen, panStart, visBp) {{
  const data=STRAIN_DATA[strain];
  const W=mmEl.parentElement.clientWidth-4, H=28;
  const DPR=window.devicePixelRatio||1;
  mmEl.width=W*DPR; mmEl.height=H*DPR;
  mmEl.style.width=W+'px'; mmEl.style.height=H+'px';
  const ctx=mmEl.getContext('2d');
  ctx.scale(DPR,DPR);
  ctx.clearRect(0,0,W,H);

  const scale=W/maxLen;

  // Genome bar
  ctx.fillStyle='#d0d5e0';
  ctx.fillRect(0,8,gLen*scale,12);

  // Prophages
  if(visible.pp)(data.prophage||[]).forEach(p=>{{
    ctx.fillStyle=ppColor(p)+'99';
    ctx.fillRect(p.start*scale,8,Math.max(1,(p.stop-p.start)*scale),12);
  }});

  // RGPs
  if(visible.rgp)(data.rgps||[]).filter(r=>(r.hgt_evidence_count||0)>=MIN_EVIDENCE||(data.df||[]).some(d=>d.rgp_id===r.rgp_id)).forEach(r=>{{
    const ev=r.hgt_evidence_count||0;
    const hasSpotDef=(data.df||[]).some(d=>d.rgp_id===r.rgp_id&&['spot_13','spot_7','spot_48'].includes(d.spot_id));
    ctx.fillStyle=hasSpotDef?'#0e7a5a':ev>=3?'#c0392b':ev>=2?'#e67e22':ev>=1?'#2980b9':'#95a5a6';
    ctx.fillRect(r.start*scale,12,Math.max(1,(r.stop-r.start)*scale),4);
  }});

  // Defense systems
  if(visible.df)(data.df||[]).forEach(d=>{{
    const c=defColor(d);
    ctx.fillStyle=c;
    ctx.fillRect(d.start*scale,8,Math.max(2,(d.stop-d.start)*scale),12);
  }});

  // Viewport highlight
  const vpX=panStart*scale, vpW=visBp*scale;
  ctx.fillStyle='rgba(30,95,168,0.15)';
  ctx.fillRect(vpX,0,vpW,H);
  ctx.strokeStyle='#1a5fa8';ctx.lineWidth=1.5;
  ctx.strokeRect(vpX,1,vpW,H-2);
}}

// ── Zoom / Pan interaction ────────────────────────────────────────────────────
function setupInteraction(canvasEl, minimapEl, strain, maxLen) {{
  const id=canvasEl.id;
  if (!viewState[id]) viewState[id]={{zoom:1,pan:0}};

  // Mouse wheel zoom
  canvasEl.addEventListener('wheel', e=>{{
    e.preventDefault();
    const vs=viewState[id];
    const rect=canvasEl.getBoundingClientRect();
    const mx=e.clientX-rect.left;
    const PAD=56, trackW=rect.width-PAD-8;
    const bpAtCursor=vs.pan + (mx-PAD)/( (trackW)/(maxLen/vs.zoom) );
    const factor=e.deltaY<0?1.25:0.8;
    vs.zoom=Math.max(1,Math.min(50,vs.zoom*factor));
    const newVisBp=maxLen/vs.zoom;
    vs.pan=Math.max(0,Math.min(maxLen-newVisBp, bpAtCursor - newVisBp*(mx-PAD)/trackW));
    drawTrack(canvasEl,minimapEl,strain,maxLen);
  }},{{passive:false}});

  // Click-drag pan
  let dragging=false, lastX=0;
  canvasEl.addEventListener('mousedown', e=>{{
    dragging=true; lastX=e.clientX;
    canvasEl.classList.add('panning');
  }});
  window.addEventListener('mousemove', e=>{{
    if(!dragging) return;
    const vs=viewState[id];
    const rect=canvasEl.getBoundingClientRect();
    const trackW=rect.width-56-8;
    const bpPerPx=(maxLen/vs.zoom)/trackW;
    vs.pan=Math.max(0,Math.min(maxLen-maxLen/vs.zoom, vs.pan-(e.clientX-lastX)*bpPerPx));
    lastX=e.clientX;
    drawTrack(canvasEl,minimapEl,strain,maxLen);
  }});
  window.addEventListener('mouseup',()=>{{
    dragging=false;
    canvasEl.classList.remove('panning');
  }});

  // Minimap click to jump
  if (minimapEl) {{
    minimapEl.addEventListener('click', e=>{{
      const vs=viewState[id];
      const rect=minimapEl.getBoundingClientRect();
      const mx=e.clientX-rect.left;
      const frac=mx/rect.width;
      const newPan=frac*maxLen - (maxLen/vs.zoom)/2;
      vs.pan=Math.max(0,Math.min(maxLen-maxLen/vs.zoom,newPan));
      drawTrack(canvasEl,minimapEl,strain,maxLen);
    }});
  }}

  // Zoom buttons
  const zIn=document.getElementById('zi_'+id);
  const zOut=document.getElementById('zo_'+id);
  const zReset=document.getElementById('zr_'+id);
  if(zIn)    zIn.addEventListener('click',()=>{{const vs=viewState[id];const mid=vs.pan+(maxLen/vs.zoom)/2;vs.zoom=Math.min(50,vs.zoom*1.5);const nv=maxLen/vs.zoom;vs.pan=Math.max(0,Math.min(maxLen-nv,mid-nv/2));drawTrack(canvasEl,minimapEl,strain,maxLen);}});
  if(zOut)   zOut.addEventListener('click',()=>{{const vs=viewState[id];const mid=vs.pan+(maxLen/vs.zoom)/2;vs.zoom=Math.max(1,vs.zoom/1.5);const nv=maxLen/vs.zoom;vs.pan=Math.max(0,Math.min(maxLen-nv,mid-nv/2));drawTrack(canvasEl,minimapEl,strain,maxLen);}});
  if(zReset) zReset.addEventListener('click',()=>{{viewState[id]={{zoom:1,pan:0}};drawTrack(canvasEl,minimapEl,strain,maxLen);}});
}}

// ── Change analysis ──────────────────────────────────────────────────────────
function analyzeChanges(isolates) {{
  const perIsolate=isolates.map(s=>{{
    const d=STRAIN_DATA[s];
    return d ? new Set((d.df||[]).map(f=>f.subtype)) : new Set();
  }});
  const baseline=perIsolate[0];
  return isolates.map((s,i)=>{{
    if(i===0) return {{strain:s,gained:[],lost:[],stable:[...baseline]}};
    const cur=perIsolate[i];
    return {{strain:s,
      gained:[...cur].filter(x=>!baseline.has(x)),
      lost:[...baseline].filter(x=>!cur.has(x)),
      stable:[...cur].filter(x=>baseline.has(x))}};
  }});
}}

// ── Render patient ────────────────────────────────────────────────────────────
// Fetch strain data from per-strain JSON file (with cache)
async function fetchStrainData(strain) {{
  if (STRAIN_DATA[strain]) return STRAIN_DATA[strain];
  if (_fetchCache[strain]) return _fetchCache[strain];
  try {{
    const resp = await fetch(`${{DATA_DIR}}/${{strain}}.json`);
    if (!resp.ok) throw new Error(`HTTP ${{resp.status}}`);
    const data = await resp.json();
    STRAIN_DATA[strain] = data;
    _fetchCache[strain] = data;
    return data;
  }} catch(e) {{
    console.warn(`Could not load data for ${{strain}}:`, e);
    return {{}};
  }}
}}

async function redrawAll(pid) {{
  const isolates = PATIENT_GROUPS[pid] || [];

  // Show loading indicator
  isolates.forEach(strain => {{
    const cid = 'canvas_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const c = document.getElementById(cid);
    if (c) {{
      const ctx = c.getContext('2d');
      ctx.clearRect(0,0,c.width,c.height);
      ctx.fillStyle = 'rgba(128,128,128,0.3)';
      ctx.font = '11px IBM Plex Mono,monospace';
      ctx.fillText('Loading...', 10, 20);
    }}
  }});

  // Fetch all strains in parallel
  await Promise.all(isolates.map(s => fetchStrainData(s)));

  // Now draw with loaded data
  const maxLen = Math.max(...isolates.map(s => STRAIN_DATA[s]?.length || 5e6));
  isolates.forEach(strain => {{
    const cid = 'canvas_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const mid = 'mini_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const c = document.getElementById(cid);
    const m = document.getElementById(mid);
    if(c) drawTrack(c, m, strain, maxLen);
  }});
}}

function renderPatient(patientId) {{
  const isolates=PATIENT_GROUPS[patientId]||[];
  const maxLen=Math.max(...isolates.map(s=>STRAIN_DATA[s]?.length||5e6));
  const changes=analyzeChanges(isolates);

  // Patient summary
  const totalDF=new Set(isolates.flatMap(s=>(STRAIN_DATA[s]?.df||[]).map(d=>d.subtype))).size;
  const totalPP=isolates.reduce((n,s)=>n+(STRAIN_DATA[s]?.prophage||[]).length,0);
  const hasChanges=changes.some(c=>c.gained.length>0||c.lost.length>0);
  document.getElementById('patientSummary').innerHTML=`
    <div class="ps-item">Patient <span>${{patientId}}</span></div>
    <div class="ps-item">Isolates <span>${{isolates.length}}</span></div>
    <div class="ps-item">Unique DF subtypes <span>${{totalDF}}</span></div>
    <div class="ps-item">Prophage regions <span>${{totalPP}}</span></div>
    <div class="ps-item">Defense changes <span>${{hasChanges?'Yes':'None detected'}}</span></div>
  `;

  const list=document.getElementById('isolateList');
  list.innerHTML='';

  isolates.forEach((strain,idx)=>{{
    const data=STRAIN_DATA[strain];
    if(!data) return;
    const asm=data.asm||{{}};
    const q=asm.quality||'unknown';
    const ch=changes[idx];
    const compPct=asm.completeness?parseFloat(asm.completeness).toFixed(1)+'%':'—';
    const cid='canvas_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const mid='mini_'+strain.replace(/[^a-zA-Z0-9]/g,'_');

    let badges='';
    if(idx===0) badges='<span class="change-badge stable">baseline</span>';
    else {{
      ch.gained.forEach(t=>badges+=`<span class="change-badge gained">+${{t}}</span>`);
      ch.lost.forEach(t=>badges+=`<span class="change-badge lost">-${{t}}</span>`);
      if(!ch.gained.length&&!ch.lost.length) badges='<span class="change-badge stable">no change</span>';
    }}

    const block=document.createElement('div');
    block.className='isolate-block';
    block.innerHTML=`
      <div class="isolate-header">
        <div class="isolate-name">${{strain}}</div>
        <span class="asm-pill pill-${{q}}">${{q.replace(/_/g,' ')}}</span>
        <div class="iso-stats">
          <div class="iso-stat">Contigs <span>${{asm.contigs||'—'}}</span></div>
          <div class="iso-stat">Complete <span>${{compPct}}</span></div>
          <div class="iso-stat">RGPs <span>${{(data.rgps||[]).length}}</span></div>
          <div class="iso-stat">DF <span>${{(data.df||[]).length}}</span></div>
          <div class="iso-stat">PADLOC <span>${{(data.padloc||[]).length}}</span></div>
          <div class="iso-stat">Prophage <span>${{(data.prophage||[]).length}}</span></div>
        </div>
        <div class="zoom-controls">
          <div class="zoom-btn" id="zo_${{cid}}" title="Zoom out">−</div>
          <div class="zoom-label" id="zl_${{cid}}">1×</div>
          <div class="zoom-btn" id="zi_${{cid}}" title="Zoom in">+</div>
          <div class="zoom-btn" id="zr_${{cid}}" title="Reset zoom" style="font-size:11px;width:36px">Reset</div>
        </div>
        <div style="margin-left:8px;display:flex;gap:4px;flex-wrap:wrap">${{badges}}</div>
      </div>
      <div class="minimap-wrap"><canvas id="${{mid}}" class="minimap" style="height:28px"></canvas></div>
      <div class="track-wrap"><canvas id="${{cid}}"></canvas></div>
    `;
    list.appendChild(block);

    requestAnimationFrame(()=>{{
      const c=document.getElementById(cid);
      const m=document.getElementById(mid);
      if(c){{
        drawTrack(c,m,strain,maxLen);
        setupInteraction(c,m,strain,maxLen);
      }}
    }});
  }});
}}

// ── Tooltip ──────────────────────────────────────────────────────────────────
const tooltip=document.getElementById('tooltip');
document.addEventListener('mousemove',e=>{{
  const canvas=e.target.closest('canvas:not(.minimap)');
  if(!canvas){{tooltip.style.display='none';return;}}
  const hits=canvasHits[canvas.id];
  if(!hits){{tooltip.style.display='none';return;}}
  const rect=canvas.getBoundingClientRect();
  const mx=e.clientX-rect.left,my=e.clientY-rect.top;
  let hit=null;
  for(const r of hits){{if(r.type==='genome')continue;if(mx>=r.x&&mx<=r.x+r.w&&my>=r.y&&my<=r.y+r.h){{hit=r;break;}}}}
  if(!hit)for(const r of hits){{if(r.type==='genome'&&mx>=r.x&&mx<=r.x+r.w&&my>=r.y&&my<=r.y+r.h){{hit=r;break;}}}}
  if(!hit){{tooltip.style.display='none';return;}}
  const d=hit.data;
  let html='';
  if(hit.type==='genome'){{
    const bp=Math.round(d.bpPerPx*(mx-56)+d.panClamped)||0;
    const bpFmt=bp>=1e6?(bp/1e6).toFixed(3)+' Mb':bp>=1e3?(bp/1e3).toFixed(1)+' kb':bp+' bp';
    html=`<div class="tt-title">Chromosome — ${{d.strain}}</div>
          <div class="tt-row">Position: <span>${{bpFmt}}</span></div>
          <div class="tt-row">Genome %: <span>${{(100*bp/d.length).toFixed(1)}}%</span></div>
          <div class="tt-row">Total length: <span>${{(d.length/1e6).toFixed(3)}} Mb</span></div>`;
    const nearby=[...(STRAIN_DATA[d.strain]?.df||[]),...(STRAIN_DATA[d.strain]?.prophage||[])]
      .filter(f=>bp>=f.start-5000&&bp<=f.stop+5000).slice(0,4);
    if(nearby.length){{html+=`<div class="tt-row" style="margin-top:5px;border-top:1px solid #eee;padding-top:4px">Nearby:</div>`;
      nearby.forEach(f=>html+=`<div class="tt-row" style="padding-left:8px">• <span>${{f.type||f.prophage_id}} (${{f.tool||'Prophage'}})</span></div>`);}}
  }} else if(hit.type==='catalog_island'){{
    const d=hit.data;
    html=`<div class="tt-title">${{d.island_id||''}}<span class="tt-badge" style="background:#1D9E7522;color:#0e7a5a">${{d.status==='nested_child'?'Nested GI':'Genomic Island'}}</span></div>
          <div class="tt-row">Group: <span class="tt-val">${{d.group_id||'—'}}</span></div>
          <div class="tt-row">Cargo: <span class="tt-val">${{d.dominant_cargo||'unknown'}}</span></div>
          <div class="tt-row">Depth: <span class="tt-val">${{d.depth===0?'Top-level':'Nested (depth '+d.depth+')'}}</span></div>
          <div class="tt-row">Evidence: <span class="tt-val">${{d.n_evidence}} lines</span></div>
          <div class="tt-row">Age: <span class="tt-val">${{d.age||'—'}}</span></div>
          <div class="tt-row">Coords: <span class="tt-val">${{d.start.toLocaleString()}}–${{d.end.toLocaleString()}}</span></div>
          <div class="tt-row">Genes: <span class="tt-val">${{d.n_genes}} (${{d.dominant_cargo}})</span></div>
          ${{d.status==='nested_child'&&d.parent?'<div class="tt-row">Parent: <span class="tt-val">'+d.parent+'</span></div>':''}}`;
  }} else if(hit.type==='fixed_defense'){{
    html=`<div class="tt-title">${{d.subtype}}<span class="tt-badge" style="background:#c0392b22;color:#c0392b">Fixed defense</span></div>
          <div class="tt-row">Age: <span class="tt-val">${{d.age_estimate||'—'}}</span></div>
          <div class="tt-row">Island: <span class="tt-val">${{(d.island_start||0).toLocaleString()}}–${{(d.island_end||0).toLocaleString()}}</span></div>
          <div class="tt-row">Confidence: <span class="tt-val">${{d.confidence||'—'}}</span></div>`;
  }} else if(hit.type==='prophage'){{
    const sz=d.stop-d.start;
    html=`<div class="tt-title">${{d.prophage_id}}<span class="tt-badge tt-pp">Prophage</span></div>
          <div class="tt-row">Confidence: <span>${{d.confidence}}</span></div>
          <div class="tt-row">Max pred: <span>${{d.max_pred}}</span></div>
          <div class="tt-row">Phage homology: <span>${{d.max_phage}}%</span></div>
          <div class="tt-row">Genes: <span>${{d.n_genes}}</span></div>
          <div class="tt-row">Size: <span>${{sz>=1000?(sz/1000).toFixed(1)+' kb':sz+' bp'}}</span></div>
          <div class="tt-row">Coords: <span>${{d.start.toLocaleString()}}–${{d.stop.toLocaleString()}}</span></div>`;
  }} else if(hit.type==='defense'){{
    const tc=d.tool==='DefenseFinder'?'tt-df':d.tool==='PADLOC'?'tt-padloc':'tt-dp';
    html=`<div class="tt-title">${{d.sys_id}}<span class="tt-badge ${{tc}}">${{d.tool}}</span></div>
          <div class="tt-row">Type: <span>${{d.type}}</span></div>
          <div class="tt-row">Subtype: <span>${{(d.subtype||'').substring(0,40)}}</span></div>`;
    if(d.tool==='DefenseFinder'){{
      const sb=d.spot_id&&KEY_SPOTS.has(d.spot_id)?`<span style="color:${{SPOT_COLORS[d.spot_id]}};font-weight:600">${{d.spot_id}}</span>`:(d.spot_id||'—');
      html+=`<div class="tt-row">In RGP: <span>${{d.in_rgp}}</span></div><div class="tt-row">Spot: ${{sb}}</div>`;
    }}
    if(d.tool==='PADLOC') html+=`<div class="tt-row">Proteins: <span>${{(d.proteins||'').substring(0,50)}}</span></div>`;
    if(d.tool==='DefensePredictor'){{html+=`<div class="tt-row">Tier: <span>${{d.tiers}}</span></div><div class="tt-row">Pfam: <span>${{d.pfams}}</span></div><div class="tt-row">Prob: <span>${{d.mean_prob}}</span></div>`;}}
    const sz=d.stop-d.start;
    html+=`<div class="tt-row">Size: <span>${{sz>=1000?(sz/1000).toFixed(1)+' kb':sz+' bp'}}</span></div>`;
  }} else if(hit.type==='rgp'){{
    const sz=d.stop-d.start;
    const ev=d.hgt_evidence_count||0;
    const ec=ev>=3?'#c0392b':ev>=2?'#e67e22':'#5a6380';
    html=`<div class="tt-title">${{d.rgp_id}}</div>
          <div class="tt-row">Size: <span>${{sz>=1000?(sz/1000).toFixed(1)+' kb':sz+' bp'}}</span></div>
          <div class="tt-row">GC foreignness: <span>${{d.gc_foreignness||'—'}}</span></div>
          <div class="tt-row">tRNA flanked: <span style="color:${{d.trna_flanked==='Yes'?'#0e7a5a':'inherit'}}">${{d.trna_flanked||'No'}}</span></div>
          ${{d.trna_product?`<div class="tt-row">tRNA: <span>${{d.trna_product}}</span></div>`:''}}
          ${{d.min_trna_dist?`<div class="tt-row">tRNA dist: <span>${{parseInt(d.min_trna_dist).toLocaleString()}} bp</span></div>`:''}}
          <div class="tt-row">Mobility gene: <span>${{d.has_mobility||'No'}}</span></div>
          <div class="tt-row" style="margin-top:4px;border-top:1px solid #eee;padding-top:4px">HGT evidence: <span style="color:${{ec}};font-weight:600">${{ev}} line${{ev!==1?'s':''}}</span></div>
          ${{d.hgt_evidence?`<div class="tt-row" style="font-size:10px;color:#888">${{d.hgt_evidence.replace(/;/g,' ·')}}</div>`:''}}`;
  }}ld_patient_viewer.py

Patient longitudinal genome viewer with:
  - All patients (single and multiple isolates)
  - Pannable / zoomable genome track (scroll to zoom, drag to pan)
  - Minimap overview navigator
  - DefenseFinder, PADLOC, DefensePredictor tracks (toggleable)
  - Prophage regions track (from per-genome prediction CSVs)
  - RGP track
  - Defense system gained/lost badges vs first isolate

Usage:
    python build_patient_viewer.py \
        --json_dir genome_tracks \
        --intersection defense_rgp_intersection.tsv \
        --rgp_file abscessus_pangenome/rgp_output/regions_of_genomic_plasticity.tsv \
        --genome_stats abscessus_pangenome/genomes_statistics.tsv \
        --padloc padloc_all_systems.csv \
        --defense_predictor defense_gene_calls_tier1_2.csv \
        --prophage_dir prophage_csvs \
        --outfile patient_comparison_viewer.html
"""

import os, csv, json, re, argparse
from collections import defaultdict

SPOT_COLORS = {
    "spot_13": "#0e7a5a",
    "spot_7":  "#1a5fa8",
    "spot_48": "#a86000",
}
DEFENSE_COLORS = {
    "CBASS":"#0e7a5a","Hna":"#1a5fa8","Dnd":"#a86000",
    "RM":"#c0392b","RosmerTA":"#6c3db8","RloC":"#b5006e",
    "AbiAlpha":"#4a7c10","Wadjet":"#4a5568","Thoeris":"#0077aa",
    "BREX":"#8b4513","default":"#5a6380",
}

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json_dir",          required=True)
    p.add_argument("--intersection",      required=True)
    p.add_argument("--rgp_file",          required=True,
                   help="abscessus RGP TSV")
    p.add_argument("--rgp_file_mass",     default=None,
                   help="massiliense RGP TSV (optional)")
    p.add_argument("--spots_file_mass",   default=None,
                   help="massiliense spots TSV (optional)")
    p.add_argument("--intersection_mass", default=None,
                   help="massiliense defense intersection TSV (optional)")
    p.add_argument("--rgp_file_boll",     default=None,
                   help="bolletii RGP TSV (optional)")
    p.add_argument("--spots_file_boll",   default=None,
                   help="bolletii spots TSV (optional)")
    p.add_argument("--intersection_boll", default=None,
                   help="bolletii defense intersection TSV (optional)")
    p.add_argument("--island_catalog",      default=None,
                   help="island_viewer_data.json from build_island_catalog.py")
    p.add_argument("--scanner_islands",    default=None,
                   help="all_islands_combined.tsv from GIMa")
    p.add_argument("--denovo_defense",     default=None,
                   help="denovo_defense_intersection.tsv for fixed defense islands")
    p.add_argument("--trna_proximity_mass", default=None,
                   help="rgp_trna_proximity_massiliense.tsv (optional)")
    p.add_argument("--trna_proximity_boll", default=None,
                   help="rgp_trna_proximity_bolletii.tsv (optional)")
    p.add_argument("--min_evidence",       type=int, default=1,
                   help="Minimum HGT evidence lines to show RGP (0=all, 1=default, 2=stricter, 3=strongest)")
    p.add_argument("--trna_proximity",     default=None,
                   help="rgp_trna_proximity.tsv (optional)")
    p.add_argument("--genome_stats",      required=True)
    p.add_argument("--padloc",            default=None)
    p.add_argument("--defense_predictor", default=None)
    p.add_argument("--prophage_dir",      default=None,
                   help="Directory containing per-genome Depht HTML files (e.g. GD05.html)")
    p.add_argument("--outfile",           default="patient_comparison_viewer.html")
    return p.parse_args()


def canonical_strain(s):
    """Strip _WGS, _hybrid, _UNCUT suffixes for deduplication matching."""
    return re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=re.IGNORECASE)


def group_patients(strains):
    """
    Group ALL strains by patient ID (GD### prefix).
    Deduplicates strains that differ only by _WGS/_hybrid suffix
    (e.g. GD233A and GD233A_WGS are the same isolate).
    Prefers the suffixed version (has RGP/stats data) as the canonical key,
    keeping the shorter name only if no suffixed version exists.
    """
    # First pass: build canonical -> list of all name variants
    canonical_map = defaultdict(list)
    for s in strains:
        canonical_map[canonical_strain(s)].append(s)

    # Second pass: for each canonical group, pick one representative
    # Prefer _WGS > _hybrid > bare name
    deduped = {}
    for canon, variants in canonical_map.items():
        if len(variants) == 1:
            deduped[canon] = variants[0]
        else:
            # Prefer suffixed versions — they have more data attached
            wgs     = [v for v in variants if v.endswith('_WGS')]
            hybrid  = [v for v in variants if v.endswith('_hybrid')]
            if wgs:
                deduped[canon] = wgs[0]
            elif hybrid:
                deduped[canon] = hybrid[0]
            else:
                deduped[canon] = sorted(variants)[0]

    # Group deduplicated canonical strains by patient ID
    groups = defaultdict(list)
    for canon, representative in deduped.items():
        m = re.match(r'^(GD\d+)', canon)
        if m:
            groups[m.group(1)].append(representative)

    return {k: sorted(v) for k, v in groups.items()}


def load_genome_lengths(json_dir):
    lengths = {}
    for fname in os.listdir(json_dir):
        if not fname.endswith(".json"): continue
        strain = fname.replace(".json","")
        with open(os.path.join(json_dir, fname)) as f:
            data = json.load(f)
        lengths[strain] = data["genome_len"]
    print(f"  Loaded JSON tracks for {len(lengths)} strains")
    return lengths


def load_genome_stats(stats_file):
    stats = {}
    with open(stats_file) as f:
        header = None
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                col = {h: i for i, h in enumerate(header)}
                continue
            strain = parts[col["Genome_name"]]
            def g(key, default="NA"):
                try: return parts[col[key]]
                except: return default
            try:
                ctg = int(g("Contigs","999"))
                if ctg==1:      q="complete"
                elif ctg<=10:   q="high"
                elif ctg<=30:   q="moderate"
                elif ctg<=80:   q="fragmented"
                else:           q="highly_fragmented"
            except: q="unknown"
            stats[strain] = {
                "contigs":g("Contigs"),"genes":g("Genes"),
                "completeness":g("Completeness"),"contamination":g("Contamination"),
                "persistent":g("Persistent_families"),"shell":g("Shell_families"),
                "cloud":g("Cloud_families"),"frag_genes":g("Fragmented_genes"),
                "quality":q,
            }
    print(f"  Loaded assembly stats for {len(stats)} strains")
    return stats


def load_rgps(rgp_file):
    rgps = defaultdict(list)
    with open(rgp_file) as f:
        raw_header = f.readline().strip().split("\t")
    col = {h: i for i, h in enumerate(raw_header)}
    with open(rgp_file) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 8: continue
            genome = row[col["genome"]]
            rgps[genome].append({
                "rgp_id": row[col["region"]],
                "start": int(row[col["start"]]),
                "stop":  int(row[7]),
                "score": row[col["score"]] if "score" in col else "",
            })
    print(f"  Loaded RGPs for {len(rgps)} strains")
    return rgps


def load_defensefinder(intersection_file):
    df = defaultdict(list)
    with open(intersection_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["in_rgp"] not in ("Yes","No"): continue
            try:
                start = int(row["sys_start_coord"])
                stop  = int(row["sys_stop_coord"])
            except: continue
            df[row["strain"]].append({
                "sys_id":row["sys_id"],"type":row["type"],"subtype":row["subtype"],
                "start":start,"stop":stop,"in_rgp":row["in_rgp"],
                "rgp_id":row.get("rgp_id",""),"spot_id":row.get("spot_id",""),
                "tool":"DefenseFinder",
            })
    print(f"  Loaded DefenseFinder for {len(df)} strains")
    return df


def load_padloc(padloc_file):
    systems = defaultdict(lambda: {"starts":[],"stops":[],"proteins":[],"system":"","strain":""})
    with open(padloc_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            strain = row["strain"]
            key    = (strain, row["system"], row["system.number"])
            systems[key]["strain"]  = strain
            systems[key]["system"]  = row["system"]
            systems[key]["starts"].append(int(row["start"]))
            systems[key]["stops"].append(int(row["end"]))
            systems[key]["proteins"].append(row["protein.name"])
    pl = defaultdict(list)
    for (strain,system,snum), d in systems.items():
        pl[strain].append({
            "sys_id":f"{strain}_PADLOC_{system}_{snum}",
            "type":system,"subtype":system,
            "start":min(d["starts"]),"stop":max(d["stops"]),
            "proteins":", ".join(set(d["proteins"])),"tool":"PADLOC",
        })
    print(f"  Loaded PADLOC for {len(pl)} strains")
    return pl


def load_defense_predictor(dp_file):
    by_strain_cat = defaultdict(list)
    with open(dp_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                start = int(float(row["start"]))
                stop  = int(float(row["end"]))
            except: continue
            by_strain_cat[(row["strain"], row["defense_cat"])].append({
                "start":start,"stop":stop,"tier":row["confidence_tier"],
                "pfam":row["pfam_name"],"desc":row["defense_desc"],
                "prob":row["probability"],
            })
    dp = defaultdict(list)
    GAP = 10000
    for (strain, cat), genes in by_strain_cat.items():
        genes.sort(key=lambda g: g["start"])
        clusters, current = [], [genes[0]]
        for g in genes[1:]:
            if g["start"] - current[-1]["stop"] <= GAP: current.append(g)
            else: clusters.append(current); current = [g]
        clusters.append(current)
        for i, cluster in enumerate(clusters):
            probs = [float(g["prob"]) for g in cluster if g["prob"]]
            dp[strain].append({
                "sys_id":f"{strain}_DP_{cat}_{i+1}","type":cat,
                "subtype":cluster[0]["desc"],
                "start":min(g["start"] for g in cluster),
                "stop":max(g["stop"]  for g in cluster),
                "tiers":", ".join(set(g["tier"] for g in cluster)),
                "pfams":", ".join(set(g["pfam"] for g in cluster)),
                "mean_prob":round(sum(probs)/len(probs),1) if probs else 0,
                "n_genes":len(cluster),"tool":"DefensePredictor",
            })
    print(f"  Loaded DefensePredictor for {len(dp)} strains")
    return dp


def load_prophages(prophage_dir):
    """
    Load per-genome Depht prophage HTML files.
    Extracts genome-absolute coordinates from the summary table.
    HTML table format: Prophage Name | Left Coordinate | Right Coordinate | Length
    Falls back to CSV-based clustering if no HTML found.
    """
    import re

    prophages = defaultdict(list)
    n_files = 0

    for fname in os.listdir(prophage_dir):
        if not fname.endswith(".html"):
            continue
        strain = fname.replace(".html", "")
        fpath  = os.path.join(prophage_dir, fname)
        try:
            with open(fpath) as f:
                html = f.read()
            # Extract table rows: Prophage Name | Left | Right | Length
            rows = re.findall(r'<tr[^>]*>(.*?)</tr>', html, re.DOTALL)
            for row in rows:
                cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.DOTALL)
                clean = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]
                if len(clean) < 3:
                    continue
                # Skip header rows
                if clean[0].lower() in ('prophage name', 'name', 'prophage'):
                    continue
                try:
                    name  = clean[0]
                    start = int(clean[1])
                    stop  = int(clean[2])
                    length= int(clean[3]) if len(clean) > 3 else stop - start
                    # Confidence based on length — Depht high-confidence regions
                    # tend to be >20kb; use length as proxy since scores aren't in HTML
                    conf = "high" if length >= 20000 else "moderate"
                    prophages[strain].append({
                        "prophage_id": name,
                        "start":       start,
                        "stop":        stop,
                        "n_genes":     0,
                        "max_pred":    1.0,
                        "max_phage":   100.0,
                        "confidence":  conf,
                        "length":      length,
                    })
                except (ValueError, IndexError):
                    continue
            n_files += 1
        except Exception as e:
            print(f"  Warning: could not read {fname}: {e}")
            continue

    print(f"  Loaded Depht prophage HTML for {n_files} strains "
          f"({sum(len(v) for v in prophages.values())} regions total)")
    return prophages


def load_trna_proximity(trna_file):
    """
    Load tRNA proximity data.
    Returns dict: rgp_id -> {trna_flanked, min_trna_dist, trna_product,
                              hgt_evidence_count, hgt_evidence, gc_foreignness,
                              has_mobility, mobility_types}
    """
    trna_data = {}
    if not trna_file or not os.path.exists(trna_file):
        return trna_data
    with open(trna_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            prod = row.get("trna_left_product") or row.get("trna_right_product") or ""
            dist = row.get("min_trna_dist") or ""
            trna_data[row["rgp_id"]] = {
                "trna_flanked":      row["trna_flanked"],
                "min_trna_dist":     dist,
                "trna_product":      prod,
                "hgt_evidence_count":int(row.get("hgt_evidence_count") or 0),
                "hgt_evidence":      row.get("hgt_evidence",""),
                "gc_foreignness":    row.get("gc_foreignness",""),
                "has_mobility":      row.get("has_mobility","No"),
                "mobility_types":    row.get("mobility_types",""),
            }
    print(f"  Loaded tRNA proximity data for {len(trna_data)} RGPs")
    return trna_data


# ── Scanner island loaders ────────────────────────────────────────────────────
def load_island_catalog(path):
    """Load island_viewer_data.json from build_island_catalog.py.
    Returns dict: strain -> list of island dicts with nesting depth.
    """
    if not path or not os.path.exists(path):
        return {}
    import json as _json
    with open(path) as f:
        data = _json.load(f)
    total = sum(len(v) for v in data.values())
    print(f"  Island catalog: {total:,} islands across {len(data)} strains")
    return data


def load_scanner_islands(path):
    if not path or not os.path.exists(path):
        return {}
    from collections import defaultdict as dd
    islands_by_strain = dd(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                islands_by_strain[row["strain"]].append({
                    "start":       int(row["start"]),
                    "end":         int(row["end"]),
                    "length":      int(row["length"]),
                    "rgp_seed":    row.get("rgp_seed",""),
                    "confidence":  row.get("confidence",""),
                    "n_evidence":  int(row.get("n_evidence",0)),
                    "evidence":    row.get("evidence",""),
                    "age_estimate":row.get("age_estimate",""),
                    "cai_ratio":   row.get("cai_ratio",""),
                    "trna_flanked":row.get("trna_flanked","No"),
                    "has_dr":      row.get("has_dr","No"),
                    "mob_types":   row.get("mob_types",""),
                    "is_denovo":   not bool(row.get("rgp_seed","")),
                })
            except: continue
    total = sum(len(v) for v in islands_by_strain.values())
    print(f"  Scanner islands: {total:,} across {len(islands_by_strain)} strains")
    return dict(islands_by_strain)


def load_denovo_defense(path):
    if not path or not os.path.exists(path):
        return {}
    from collections import defaultdict as dd
    defense_by_strain = dd(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("in_denovo_island") != "Yes": continue
            try:
                defense_by_strain[row["strain"]].append({
                    "subtype":     row["subtype"],
                    "sys_start":   int(row["sys_start"])    if row["sys_start"]    else 0,
                    "sys_stop":    int(row["sys_stop"])     if row["sys_stop"]     else 0,
                    "island_start":int(row["island_start"]) if row["island_start"] else 0,
                    "island_end":  int(row["island_end"])   if row["island_end"]   else 0,
                    "confidence":  row.get("island_confidence",""),
                    "age_estimate":row.get("age_estimate",""),
                })
            except: continue
    total = sum(len(v) for v in defense_by_strain.values())
    print(f"  Fixed defense: {total:,} across {len(defense_by_strain)} strains")
    return dict(defense_by_strain)


def build_html(genome_lengths, rgps, df, pl, dp, prophages, stats, patient_groups, trna_data=None, min_evidence=1, scanner_islands=None, denovo_defense=None, island_catalog=None):
    # Build strain_data from UNION of all data sources
    # so strains with prophage/defense data but no JSON track are still included
    all_strains = set(genome_lengths.keys())
    all_strains |= set(rgps.keys())
    all_strains |= set(df.keys())
    all_strains |= set(pl.keys())
    all_strains |= set(dp.keys())
    all_strains |= set(prophages.keys())
    all_strains |= set(stats.keys())

    # Build canonical -> all variants map for data merging
    canon_variants = defaultdict(list)
    for s in all_strains:
        canon_variants[canonical_strain(s)].append(s)

    strain_data = {}
    for strain in all_strains:
        asm = stats.get(strain, {})
        fallback_len = 5000000
        try:
            genes = int(asm.get("genes", 0) or 0)
            if genes > 0:
                fallback_len = genes * 1000
        except: pass

        # Merge data from all name variants of this strain
        # e.g. GD233A_WGS gets defense data from GD233A too
        canon = canonical_strain(strain)
        variants = canon_variants.get(canon, [strain])

        merged_df      = []
        merged_padloc  = []
        merged_dp      = []
        merged_rgps    = rgps.get(strain, [])
        merged_prophage= prophages.get(strain, [])
        merged_asm     = asm

        for v in variants:
            merged_df     += df.get(v, [])
            merged_padloc += pl.get(v, [])
            merged_dp     += dp.get(v, [])
            if not merged_rgps:
                merged_rgps = rgps.get(v, [])
            if not merged_prophage:
                merged_prophage = prophages.get(v, [])
            if not merged_asm and stats.get(v):
                merged_asm = stats.get(v, {})

        # Deduplicate by sys_id
        seen = set()
        dedup_df = []
        for d in merged_df:
            if d["sys_id"] not in seen:
                seen.add(d["sys_id"]); dedup_df.append(d)

        seen = set()
        dedup_pl = []
        for d in merged_padloc:
            if d["sys_id"] not in seen:
                seen.add(d["sys_id"]); dedup_pl.append(d)

        seen = set()
        dedup_dp = []
        for d in merged_dp:
            if d["sys_id"] not in seen:
                seen.add(d["sys_id"]); dedup_dp.append(d)

        # Add tRNA proximity data to each RGP
        rgps_with_trna = []
        for rgp in merged_rgps:
            rgp_copy = dict(rgp)
            td = (trna_data or {}).get(rgp["rgp_id"], {})
            rgp_copy["trna_flanked"]       = td.get("trna_flanked", "No")
            rgp_copy["min_trna_dist"]      = td.get("min_trna_dist", "")
            rgp_copy["trna_product"]       = td.get("trna_product", "")
            rgp_copy["hgt_evidence_count"] = td.get("hgt_evidence_count", 0)
            rgp_copy["hgt_evidence"]       = td.get("hgt_evidence", "")
            rgp_copy["gc_foreignness"]     = td.get("gc_foreignness", "")
            rgp_copy["has_mobility"]       = td.get("has_mobility", "No")
            rgp_copy["mobility_types"]     = td.get("mobility_types", "")
            rgps_with_trna.append(rgp_copy)

        # Filter RGPs by minimum evidence level
        rgps_filtered = [r for r in rgps_with_trna
                         if r.get("hgt_evidence_count", 0) >= min_evidence
                         or r.get("has_defense_in_rgp", False)]

        # Scanner islands for this strain
        def _canon(s):
            import re as _re
            return _re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=_re.IGNORECASE)
        scan_isls = None
        if scanner_islands:
            scan_isls = scanner_islands.get(strain) or scanner_islands.get(_canon(strain))
            if not scan_isls:
                sc = _canon(strain)
                for k in scanner_islands:
                    if _canon(k) == sc:
                        scan_isls = scanner_islands[k]; break
        scan_isls = scan_isls or []
        fix_def = None
        if denovo_defense:
            fix_def = denovo_defense.get(strain) or denovo_defense.get(_canon(strain))
            if not fix_def:
                sc = _canon(strain)
                for k in denovo_defense:
                    if _canon(k) == sc:
                        fix_def = denovo_defense[k]; break
        fix_def = fix_def or []

        # Island catalog data (nested-resolved)
        cat_isls = None
        if island_catalog:
            cat_isls = island_catalog.get(strain) or island_catalog.get(_canon(strain))
            if not cat_isls:
                sc = _canon(strain)
                for k in island_catalog:
                    if _canon(k) == sc:
                        cat_isls = island_catalog[k]; break
        cat_isls = cat_isls or []

        # Build defense coordinate list for annotation
        # (used to flag which islands carry defense systems)
        all_defense = []
        for d in dedup_df:
            s = d.get("start")
            e = d.get("stop") or d.get("end")
            if s and e:
                all_defense.append({"start": s, "end": e})
        for d in fix_def:
            s = d.get("sys_start") or d.get("island_start", 0)
            e = d.get("sys_stop")  or d.get("island_end",   0)
            if s and e:
                all_defense.append({"start": s, "end": e})

        def has_defense_overlap(isl, defense_list):
            for d in defense_list:
                ds = d.get("start", 0)
                de = d.get("stop") or d.get("end") or 0
                if ds and de:
                    if isl["start"] <= de and isl["end"] >= ds:
                        return True
            return False

        # Show ALL catalog islands — flag which ones carry defense systems
        # Genomic islands are defined by HGT origin, not cargo type
        defense_cat_isls = []
        for isl in cat_isls:
            isl_copy = dict(isl)
            isl_copy["has_defense"] = has_defense_overlap(isl, all_defense)
            defense_cat_isls.append(isl_copy)

        sdata = {
            "length":         genome_lengths.get(strain) or genome_lengths.get(_canon(strain)) or fallback_len,
            "rgps":           rgps_filtered,
            "scanner_islands":scan_isls,
            "fixed_defense":  fix_def,
            "catalog_islands":defense_cat_isls,
            "df":       dedup_df,
            "padloc":   dedup_pl,
            "dp":       dedup_dp,
            "prophage": merged_prophage,
            "asm":      merged_asm,
        }
        strain_data[strain] = sdata

    # All patients — single and multi isolate
    patient_index = {}
    for pid, isos in sorted(patient_groups.items(),
                             key=lambda x: (x[0].replace("GD","").zfill(6))):
        valid = [s for s in isos if s in strain_data]
        if valid:
            patient_index[pid] = valid

    patient_json        = json.dumps(patient_index)
    # strain_data is written to per-strain JSON files — not embedded in HTML
    spot_colors_json    = json.dumps(SPOT_COLORS)
    defense_colors_json = json.dumps(DEFENSE_COLORS)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Patient longitudinal defense viewer</title>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap" rel="stylesheet">
<style>
  :root {{
    --bg:#f4f6fa;--surface:#fff;--border:#dde1ea;
    --text:#1a1f2e;--muted:#5a6380;--accent:#0e7a5a;
    --font-mono:'IBM Plex Mono',monospace;--font-sans:'IBM Plex Sans',sans-serif;
  }}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{background:var(--bg);color:var(--text);font-family:var(--font-sans);font-size:14px}}
  header{{border-bottom:1px solid var(--border);padding:1.25rem 2rem;
    display:flex;align-items:center;gap:1.5rem;background:var(--surface)}}
  header h1{{font-size:14px;font-weight:500;letter-spacing:.08em;text-transform:uppercase;
    color:var(--accent);font-family:var(--font-mono)}}
  header span{{font-size:12px;color:var(--muted);font-family:var(--font-mono)}}
  .controls{{padding:1rem 2rem;display:flex;align-items:center;gap:1rem;
    border-bottom:1px solid var(--border);flex-wrap:wrap;background:var(--surface)}}
  .controls label{{font-size:11px;text-transform:uppercase;letter-spacing:.08em;
    color:var(--muted);font-family:var(--font-mono)}}
  select{{background:var(--bg);border:1px solid var(--border);color:var(--text);
    padding:6px 12px;border-radius:4px;font-family:var(--font-mono);
    font-size:13px;cursor:pointer;min-width:200px}}
  select:focus{{outline:none;border-color:var(--accent)}}
  .toggles{{display:flex;gap:6px;align-items:center;flex-wrap:wrap}}
  .toggle-btn{{display:flex;align-items:center;gap:5px;padding:4px 10px;border-radius:20px;
    border:1.5px solid;font-family:var(--font-mono);font-size:11px;cursor:pointer;
    font-weight:500;transition:opacity .15s;user-select:none}}
  .toggle-btn.df    {{border-color:#0e7a5a;color:#0e7a5a;background:#f0faf6}}
  .toggle-btn.padloc{{border-color:#1a5fa8;color:#1a5fa8;background:#f0f5ff}}
  .toggle-btn.dp    {{border-color:#8b4513;color:#8b4513;background:#fff8f0}}
  .toggle-btn.rgp   {{border-color:#888780;color:#5a6380;background:#f5f5f5}}
  .toggle-btn.pp    {{border-color:#7c2d92;color:#7c2d92;background:#faf0ff}}
  .toggle-btn.off   {{opacity:.3}}
  .toggle-dot{{width:7px;height:7px;border-radius:50%;flex-shrink:0}}
  .df .toggle-dot{{background:#0e7a5a}}.padloc .toggle-dot{{background:#1a5fa8}}
  .dp .toggle-dot{{background:#8b4513}}.rgp .toggle-dot{{background:#888780}}
  .pp .toggle-dot{{background:#7c2d92}}
  .patient-summary{{padding:.6rem 2rem;border-bottom:1px solid var(--border);
    background:#f8fafd;display:flex;gap:2rem;align-items:center;flex-wrap:wrap}}
  .ps-item{{font-family:var(--font-mono);font-size:12px;color:var(--muted)}}
  .ps-item span{{color:var(--text);font-weight:500}}
  .change-badge{{display:inline-block;font-size:10px;font-family:var(--font-mono);
    padding:1px 7px;border-radius:3px;font-weight:600;margin:1px}}
  .gained{{background:#d4f2e8;color:#0e7a5a}}
  .lost{{background:#fde8e4;color:#a82c1a}}
  .stable{{background:#eaecf0;color:#5a6380}}
  .isolate-list{{padding:1.5rem 2rem;display:flex;flex-direction:column;gap:1.5rem}}
  .isolate-block{{background:var(--surface);border:1px solid var(--border);border-radius:10px;overflow:hidden}}
  .isolate-header{{display:flex;align-items:center;gap:.75rem;padding:.6rem 1.25rem;
    border-bottom:1px solid var(--border);background:#f8fafd;flex-wrap:wrap}}
  .isolate-name{{font-family:var(--font-mono);font-size:14px;font-weight:500;
    color:var(--text);min-width:110px}}
  .iso-stats{{display:flex;gap:.75rem;flex-wrap:wrap}}
  .iso-stat{{font-family:var(--font-mono);font-size:11px;color:var(--muted)}}
  .iso-stat span{{color:var(--text);font-weight:500}}
  .asm-pill{{font-size:10px;font-family:var(--font-mono);padding:2px 8px;border-radius:10px;font-weight:500}}
  .pill-complete,.pill-high{{background:#d4f2e8;color:#0e7a5a}}
  .pill-moderate{{background:#fff3cc;color:#8a5c00}}
  .pill-fragmented{{background:#fde8e4;color:#a82c1a}}
  .pill-highly_fragmented{{background:#fbd1cb;color:#a82c1a}}
  .pill-unknown{{background:#eaecf0;color:#5a6380}}
  /* Zoom controls */
  .zoom-controls{{display:flex;gap:4px;align-items:center;margin-left:auto}}
  .zoom-btn{{background:var(--surface);border:1px solid var(--border);color:var(--text);
    width:26px;height:26px;border-radius:4px;cursor:pointer;font-size:14px;
    display:flex;align-items:center;justify-content:center;font-family:var(--font-mono);
    user-select:none}}
  .zoom-btn:hover{{background:var(--bg);border-color:var(--accent)}}
  .zoom-label{{font-family:var(--font-mono);font-size:11px;color:var(--muted);min-width:38px;text-align:center}}
  /* Track canvas wrapper */
  .track-wrap{{padding:.75rem 1.25rem 1rem;position:relative}}
  canvas{{display:block;cursor:grab}}
  canvas.panning{{cursor:grabbing}}
  /* Minimap */
  .minimap-wrap{{padding:0 1.25rem .75rem;}}
  canvas.minimap{{cursor:pointer;border-radius:3px;border:1px solid var(--border)}}
  .legend{{display:flex;flex-wrap:wrap;gap:.75rem;padding:.75rem 2rem;
    border-top:1px solid var(--border);background:var(--surface)}}
  .legend-item{{display:flex;align-items:center;gap:5px;font-size:11px;
    color:var(--muted);font-family:var(--font-mono)}}
  .legend-swatch{{width:11px;height:11px;border-radius:2px;flex-shrink:0}}
  .legend-divider{{width:1px;height:16px;background:var(--border);margin:0 2px}}
  .tooltip{{position:fixed;background:#fff;border:1px solid #c5cad6;border-radius:6px;
    padding:10px 14px;font-family:var(--font-mono);font-size:12px;
    pointer-events:none;display:none;z-index:1000;max-width:300px;
    line-height:1.8;box-shadow:0 4px 16px rgba(0,0,0,0.12)}}
  .tt-title{{font-weight:500;color:#1a1f2e;margin-bottom:4px;font-size:13px}}
  .tt-badge{{display:inline-block;font-size:10px;padding:1px 6px;border-radius:3px;
    font-weight:600;margin-left:6px}}
  .tt-row{{color:#5a6380}}.tt-row span{{color:#1a1f2e}}
  .tt-df{{background:#d4f2e8;color:#0e7a5a}}
  .tt-padloc{{background:#ddeeff;color:#1a5fa8}}
  .tt-dp{{background:#fff0e0;color:#8b4513}}
  .tt-pp{{background:#f3e8ff;color:#7c2d92}}
  .no-data{{padding:2rem;text-align:center;color:var(--muted);font-family:var(--font-mono)}}
</style>
</head>
<body>

<header>
  <h1>Patient longitudinal viewer</h1>
  <span>M. abscessus — serial isolate defense &amp; prophage comparison</span>
</header>

<div class="controls">
  <label>Patient</label>
  <select id="patientSelect"></select>
  <div class="toggles">
    <span style="font-size:11px;color:var(--muted);font-family:var(--font-mono);
      text-transform:uppercase;letter-spacing:.07em;white-space:nowrap">Tracks:</span>
    <div class="toggle-btn df"     id="tog-df"     onclick="toggleTrack('df')"><div class="toggle-dot"></div>DefenseFinder</div>
    <div class="toggle-btn padloc" id="tog-padloc"  onclick="toggleTrack('padloc')"><div class="toggle-dot"></div>PADLOC</div>
    <div class="toggle-btn dp"     id="tog-dp"     onclick="toggleTrack('dp')"><div class="toggle-dot"></div>DefensePredictor</div>
    <div class="toggle-btn pp"     id="tog-pp"     onclick="toggleTrack('pp')"><div class="toggle-dot"></div>Prophage</div>
    <div class="toggle-btn rgp"    id="tog-rgp"    onclick="toggleTrack('rgp')"><div class="toggle-dot"></div>RGPs</div>
    <div class="toggle-btn"        id="tog-fix"    onclick="toggleTrack('fix')"
         style="background:#c0392b22;border-color:#c0392b"><div class="toggle-dot" style="background:#c0392b"></div>Fixed defense</div>
    <div class="toggle-btn"        id="tog-isl"    onclick="toggleTrack('isl')"
         style="background:#1D9E7522;border-color:#1D9E75"><div class="toggle-dot" style="background:#1D9E75"></div>Genomic islands</div>
  </div>
  <div style="display:flex;align-items:center;gap:8px;margin-left:8px">
    <span style="font-size:11px;color:var(--muted);font-family:var(--font-mono);text-transform:uppercase;letter-spacing:.07em;white-space:nowrap">RGP evidence:</span>
    <select id="evidenceFilter" onchange="MIN_EVIDENCE=parseInt(this.value);const pid=document.getElementById('patientSelect').value;if(pid)redrawAll(pid);"
      style="font-size:11px;padding:4px 8px;min-width:0">
      <option value="0">All RGPs</option>
      <option value="1" selected>1+ lines</option>
      <option value="2">2+ lines</option>
      <option value="3">3+ lines (strong)</option>
      <option value="4">4 lines (triple+defense)</option>
    </select>
  </div>
</div>

<div class="patient-summary" id="patientSummary"></div>
<div class="isolate-list"    id="isolateList"></div>

<div class="legend">
  <div class="legend-item"><div class="legend-swatch" style="background:#0e7a5a"></div>spot_13 CBASS_II</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#1a5fa8"></div>spot_7 Hna</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#a86000"></div>spot_48 Dnd</div>
  <div class="legend-divider"></div>
  <div class="legend-item"><div class="legend-swatch" style="background:#0e7a5a;opacity:.5"></div>DF</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#1a5fa8"></div>PADLOC</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#8b4513"></div>DefPredictor</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#7c2d92"></div>Prophage (high)</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#b57cc4"></div>Prophage (mod)</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#95a5a6"></div>RGP — 1 evidence line</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#2980b9"></div>RGP — tRNA flanked</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#e67e22"></div>RGP — 2 evidence lines</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#c0392b"></div>RGP — 3–4 lines (strongest)</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#0e7a5a"></div>GI — defense cargo</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#7c2d92"></div>GI — mobility cargo</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#c0392b"></div>GI — TA system cargo</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#378ADD"></div>GI — other cargo</div>
  <div class="legend-item" style="font-size:10px;color:var(--muted)">Dashed border = nested child island · GI+1 label = nesting depth</div>
  <div class="legend-divider"></div>
  <div class="legend-item"><div class="legend-swatch" style="background:#c0392b"></div>Fixed defense — very recent</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#e67e22"></div>Fixed defense — recent</div>
  <div class="legend-item"><div class="legend-swatch" style="background:#3498db"></div>Fixed defense — moderate/old</div>
  <div class="legend-divider"></div>
  <div class="legend-item" style="font-size:11px;color:var(--muted);font-family:var(--font-mono)">HGT lines: GC foreign · tRNA flanked · mobility gene · defense system</div>
  <div class="legend-divider"></div>
  <div class="legend-item"><span class="change-badge gained">+gained</span></div>
  <div class="legend-item"><span class="change-badge lost">-lost</span></div>
  <div class="legend-item"><span class="change-badge stable">stable</span></div>
</div>

<div class="tooltip" id="tooltip"></div>

<script>
const STRAIN_DATA      = {{}};  // populated on demand via fetch
const PATIENT_GROUPS   = {patient_json};
const SPOT_COLORS      = {spot_colors_json};
const DEFENSE_COLORS   = {defense_colors_json};
const DATA_DIR         = 'viewer_data';  // relative path to per-strain JSON files
let   _fetchCache      = {{}};           // cache fetched strain data
const KEY_SPOTS        = new Set(['spot_13','spot_7','spot_48']);
const visible          = {{df:true,padloc:true,dp:true,pp:true,rgp:true,fix:true,isl:true}};
let   MIN_EVIDENCE     = {min_evidence};  // filter RGPs by HGT evidence count

// Per-isolate view state (zoom/pan per canvas)
const viewState = {{}};
const canvasHits = {{}};

function toggleTrack(track) {{
  visible[track] = !visible[track];
  document.getElementById('tog-'+track).classList.toggle('off', !visible[track]);
  const pid = document.getElementById('patientSelect').value;
  if (pid) redrawAll(pid);
}}

// ── Colors ───────────────────────────────────────────────────────────────────
function defColor(d) {{
  if (d.tool==='DefenseFinder') {{
    if (d.in_rgp==='Yes' && KEY_SPOTS.has(d.spot_id)) return SPOT_COLORS[d.spot_id];
    return DEFENSE_COLORS[d.type] || DEFENSE_COLORS.default;
  }}
  if (d.tool==='PADLOC')           return '#1a5fa8';
  if (d.tool==='DefensePredictor') return '#8b4513';
  return '#5a6380';
}}
function ppColor(pp) {{
  if (pp.confidence==='high')     return '#7c2d92';
  if (pp.confidence==='moderate') return '#b57cc4';
  return '#d4b8e0';
}}
function rgpColor(rgp, defInRGP) {{
  // Returns [fillColor, strokeColor] as solid colors like prophage track
  const s = defInRGP.filter(d=>d.rgp_id===rgp.rgp_id&&KEY_SPOTS.has(d.spot_id)).map(d=>d.spot_id);
  if (s.length>0) {{
    const c = SPOT_COLORS[s[0]];
    return [c, c];
  }}
  const ev = rgp.hgt_evidence_count || 0;
  if (ev >= 3) return ['#c0392b', '#922b21'];   // red — 3+ evidence
  if (ev >= 2) return ['#e67e22', '#b7600d'];   // amber — 2 evidence
  if (rgp.trna_flanked === 'Yes') return ['#2980b9','#1a5fa8']; // blue — tRNA flanked
  return ['#95a5a6', '#7f8c8d'];                // gray — 1 evidence
}}

function rgpFill(rgp, defInRGP) {{
  return rgpColor(rgp, defInRGP)[0];
}}

// ── Draw a single track canvas ────────────────────────────────────────────────
function drawTrack(canvasEl, minimapEl, strain, maxLen) {{
  const data = STRAIN_DATA[strain];
  if (!data) return;

  const gLen  = data.length;
  const rgps         = (data.rgps || []).filter(r =>
    (r.hgt_evidence_count || 0) >= MIN_EVIDENCE ||
    (data.df||[]).some(d => d.rgp_id === r.rgp_id)
  );
  const scannerIslands = data.scanner_islands || [];
  const fixedDefense   = data.fixed_defense   || [];
  const catalogIslands = data.catalog_islands || [];
  const df    = data.df     || [];
  const pl    = data.padloc || [];
  const dp    = data.dp     || [];
  const pp    = data.prophage|| [];
  const defInRGP = df.filter(d=>d.in_rgp==='Yes');

  const id    = canvasEl.id;
  if (!viewState[id]) viewState[id] = {{zoom:1, pan:0}};
  const vs    = viewState[id];

  const DPR   = window.devicePixelRatio || 1;
  const W     = canvasEl.parentElement.clientWidth - 4;
  const PAD   = 56;
  const ARROW_H=18, GENOME_H=12, RGP_H=14, PP_H=16, TRACK_GAP=8;

  // Estimate height
  const dfR  = visible.df     ? Math.min(df.length,5)  : 0;
  const plR  = visible.padloc ? Math.min(pl.length,5)  : 0;
  const dpR  = visible.dp     ? Math.min(dp.length,5)  : 0;
  const ppR  = visible.pp     ? 1 : 0;
  const aboveH = (dfR+plR+dpR)*( ARROW_H+3) + (dfR>0?TRACK_GAP:0)+(plR>0?TRACK_GAP:0)+(dpR>0?TRACK_GAP:0) + 24;
  // Dynamic height: estimate lanes needed from island count and genome length
  const catIslands  = data.catalog_islands || [];
  const ISL_ROW_H   = 12;
  const ISL_GAP     = 3;
  // Estimate max lanes: islands / (genome_length / avg_island_width_px)
  // Conservative: assume up to 4 lanes for dense genomes, 1 for sparse
  const islDensity  = catIslands.length / Math.max(1, (data.length||5000000) / 50000);
  const estLanes    = visible.isl && catIslands.length > 0
    ? Math.min(5, Math.max(1, Math.ceil(islDensity)))
    : 0;
  // Reserve minimum height for island track even when lane count is unknown
  const islTrackH   = visible.isl
    ? Math.max(4 * (ISL_ROW_H + ISL_GAP) + 8, estLanes * (ISL_ROW_H + ISL_GAP) + (estLanes > 0 ? 8 : 0))
    : 0;
  const belowH = (visible.rgp?RGP_H+6:0) + (visible.pp?PP_H+6:0) +
                 (visible.fix?8+4:0) + islTrackH + 24;
  const H = Math.max(120, aboveH + GENOME_H + belowH);

  canvasEl.width  = W*DPR; canvasEl.height = H*DPR;
  canvasEl.style.width=W+'px'; canvasEl.style.height=H+'px';
  // Ensure container expands to fit canvas on retina displays
  canvasEl.style.minHeight=H+'px';

  const ctx = canvasEl.getContext('2d');
  ctx.scale(DPR, DPR);
  ctx.clearRect(0,0,W,H);

  // Viewport: pan is in bp, zoom scales bp→px
  const trackW = W - PAD - 8;
  const visibleBp = maxLen / vs.zoom;
  const bpPerPx   = visibleBp / trackW;
  const panClamped= Math.max(0, Math.min(vs.pan, maxLen - visibleBp));
  vs.pan = panClamped;

  function bpToX(bp) {{ return PAD + (bp - panClamped) / bpPerPx; }}
  function bpWidth(bp) {{ return bp / bpPerPx; }}
  function inView(start, stop) {{ return stop >= panClamped && start <= panClamped + visibleBp; }}

  const hits = [];
  const GENOME_Y = aboveH + 8;

  // Reference line
  ctx.fillStyle='#e8eaf0';
  ctx.fillRect(PAD, GENOME_Y-1, trackW, 1);

  // Axis ticks
  ctx.font='9px IBM Plex Mono,monospace';
  const tickCount = 7;
  for (let i=0;i<=tickCount;i++) {{
    const bp = Math.round((panClamped + visibleBp*i/tickCount));
    const x  = PAD + trackW*i/tickCount;
    ctx.fillStyle='#c5cad6'; ctx.fillRect(x, GENOME_Y+GENOME_H, 1, 4);
    ctx.fillStyle='#5a6380';
    ctx.textAlign = i===0?'left':i===tickCount?'right':'center';
    const lbl=bp>=1e6?(bp/1e6).toFixed(2)+'M':bp>=1e3?(bp/1e3).toFixed(0)+'k':bp;
    ctx.fillText(lbl, x, GENOME_Y+GENOME_H+13);
  }}

  // Genome backbone — only show actual genome length
  const genomeEndX = Math.min(bpToX(gLen), PAD+trackW);
  const grad=ctx.createLinearGradient(PAD,0,genomeEndX,0);
  grad.addColorStop(0,'#d0d5e0');grad.addColorStop(0.5,'#bcc2d0');grad.addColorStop(1,'#d0d5e0');
  ctx.fillStyle=grad; ctx.strokeStyle='#a0a8bc'; ctx.lineWidth=1;
  ctx.beginPath();
  ctx.roundRect(PAD, GENOME_Y, Math.max(0,genomeEndX-PAD), GENOME_H, 3);
  ctx.fill(); ctx.stroke();
  hits.push({{type:'genome',data:{{length:gLen,strain}},
    x:PAD,y:GENOME_Y,w:genomeEndX-PAD,h:GENOME_H,bpPerPx,panClamped}});

  // RGP track
  if (visible.rgp) {{
    const RGP_Y = GENOME_Y+GENOME_H+4;
    rgps.forEach(rgp => {{
      if (!inView(rgp.start,rgp.stop)) return;
      const x=bpToX(rgp.start), w=Math.max(2,bpWidth(rgp.stop-rgp.start));
      const [fill,stroke] = rgpColor(rgp,defInRGP);
      ctx.fillStyle   = fill + 'cc';
      ctx.strokeStyle = stroke;
      ctx.lineWidth=1;
      ctx.beginPath();ctx.roundRect(x,RGP_Y,w,RGP_H,3);ctx.fill();ctx.stroke();
      hits.push({{type:'rgp',data:rgp,x,y:RGP_Y,w,h:RGP_H}});
    }});
    ctx.fillStyle='#5a6380';ctx.font='9px IBM Plex Mono,monospace';
    ctx.textAlign='right';ctx.fillText('RGPs',PAD-4,RGP_Y+RGP_H/2+3);
  }}

  // Prophage track
  if (visible.pp && pp.length>0) {{
    const PP_Y = GENOME_Y+GENOME_H+(visible.rgp?RGP_H+8:4);
    pp.forEach(p => {{
      if (!inView(p.start,p.stop)) return;
      const x=bpToX(p.start), w=Math.max(3,bpWidth(p.stop-p.start));
      const c=ppColor(p);
      ctx.fillStyle=c+'bb'; ctx.strokeStyle=c; ctx.lineWidth=1;
      ctx.beginPath();ctx.roundRect(x,PP_Y,w,PP_H,3);ctx.fill();ctx.stroke();
      if (w>40) {{
        ctx.fillStyle='#fff';ctx.font='8px IBM Plex Mono,monospace';ctx.textAlign='left';
        ctx.fillText(p.confidence+' conf',x+4,PP_Y+PP_H/2+3);
      }}
      hits.push({{type:'prophage',data:p,x,y:PP_Y,w,h:PP_H}});
    }});

    // ── Fixed defense islands track ────────────────────────────────────────
    const FIX_Y = PP_Y + PP_H + 3;
    const FIX_H = 7;
    if (visible.fix) {{
      fixedDefense.forEach(fd => {{
        const s = fd.island_start || fd.sys_start || 0;
        const e = fd.island_end   || fd.sys_stop  || 0;
        if (!s || !e || !inView(s,e)) return;
        const x = bpToX(s), w = Math.max(3, bpWidth(e-s));
        const age = fd.age_estimate || '';
        const col = age==='very_recent'?'#c0392b':age==='recent'?'#e67e22':
                    age==='moderate'?'#3498db':'#888780';
        ctx.fillStyle=col+'dd'; ctx.strokeStyle=col; ctx.lineWidth=1;
        ctx.beginPath(); ctx.roundRect(x,FIX_Y,w,FIX_H,2); ctx.fill(); ctx.stroke();
        hits.push({{type:'fixed_defense',data:fd,x,y:FIX_Y,w,h:FIX_H}});
      }});
    }}

    // ── Nested genomic island catalog tracks ─────────────────────────────────
    // Each nesting depth gets its own row. Top-level (depth 0) are full height.
    // Nested children are shorter and offset downward, visually inside the parent.
    const ISL_ROW_H  = 12;
    const ISL_GAP    = 3;
    const ISL_BASE_Y = FIX_Y + (visible.fix ? FIX_H + 6 : 0);

    function islColor(isl) {{
      // Color by dominant cargo — all islands shown
      const cargo = isl.dominant_cargo || '';
      const hasDef = isl.has_defense;
      if (hasDef)              return ['#0e7a5a','#085041'];
      if (cargo==='mobility')  return ['#7c2d92','#5b1f6e'];
      if (cargo==='ta_system') return ['#c0392b','#922b21'];
      if (cargo==='metal')     return ['#e67e22','#b7600d'];
      if (cargo==='efflux')    return ['#2980b9','#1a5fa8'];
      if (cargo==='phage')     return ['#8e44ad','#6c3483'];
      if (cargo==='regulatory')return ['#16a085','#0d7566'];
      if (cargo==='defense')   return ['#0e7a5a','#085041'];
      // Unknown/hypothetical
      const ev = isl.n_evidence || 0;
      if (ev>=4) return ['#1D9E75','#0F6E56'];
      if (ev>=3) return ['#378ADD','#185FA5'];
      return ['#888780','#5F5E5A'];
    }}

    if (visible.isl) {{
      const catIslands = (data.catalog_islands || [])
        .filter(i => inView(i.start, i.end))
        .sort((a,b) => a.start - b.start);

      if (catIslands.length > 0) {{

        // ── Lane assignment ─────────────────────────────────────────────────
        // Assign each island to the lowest lane where it doesn't overlap
        // any previously placed island. This prevents visual overlap.
        const laneEnds = [];   // laneEnds[lane] = pixel x where last island ends
        const LANE_PAD = 2;    // minimum pixel gap between islands in same lane
        const islLanes = catIslands.map(isl => {{
          const x = bpToX(isl.start);
          const w = Math.max(4, bpWidth(isl.end - isl.start));
          // Find lowest lane with room
          let lane = 0;
          while (laneEnds[lane] !== undefined && laneEnds[lane] > x - LANE_PAD) {{
            lane++;
          }}
          laneEnds[lane] = x + w;
          return lane;
        }});

        const maxLane = Math.max(...islLanes);

        // Draw islands in their assigned lanes
        catIslands.forEach((isl, idx) => {{
          const lane  = islLanes[idx];
          const x     = bpToX(isl.start);
          const w     = Math.max(4, bpWidth(isl.end - isl.start));
          const h     = Math.max(6, ISL_ROW_H - lane);   // slightly shorter in deeper lanes
          const y     = ISL_BASE_Y + lane * (ISL_ROW_H + ISL_GAP);
          const [fill, stroke] = islColor(isl);

          // Nested children get dashed border; top-level solid
          const isNested = (isl.status === 'nested_child') || lane > 0;
          ctx.fillStyle   = fill + (isNested ? '99' : 'cc');
          ctx.strokeStyle = stroke;
          ctx.lineWidth   = 1;
          ctx.setLineDash(isNested ? [2,2] : []);
          ctx.beginPath();
          const r2 = Math.min(2, w/2, h/2);
          ctx.moveTo(x+r2, y);
          ctx.lineTo(x+w-r2, y);
          ctx.arcTo(x+w, y, x+w, y+r2, r2);
          ctx.lineTo(x+w, y+h-r2);
          ctx.arcTo(x+w, y+h, x+w-r2, y+h, r2);
          ctx.lineTo(x+r2, y+h);
          ctx.arcTo(x, y+h, x, y+h-r2, r2);
          ctx.lineTo(x, y+r2);
          ctx.arcTo(x, y, x+r2, y, r2);
          ctx.closePath();
          ctx.fill(); ctx.stroke();
          ctx.setLineDash([]);

          // BLAST-validated badge — small star on top-right corner
          if (isl.blast_validated) {{
            ctx.fillStyle = '#f59e0b';
            ctx.font = 'bold 8px sans-serif';
            ctx.textAlign = 'right';
            ctx.textBaseline = 'top';
            ctx.fillText('★', x+w-1, y+1);
            ctx.textBaseline = 'middle';
          }}

          // Label wide islands
          if (w > 80) {{
            ctx.fillStyle = '#fff';
            ctx.font = `bold 8px IBM Plex Mono,monospace`;
            ctx.textBaseline = 'middle';
            ctx.textAlign = 'left';
            // Show island_id if wide enough, otherwise cargo type
            const label = w > 150
              ? (isl.island_id || isl.dominant_cargo || '').substring(0,12)
              : (isl.dominant_cargo || '').substring(0,8);
            ctx.fillText(label, x+3, y+h/2);
          }}

          hits.push({{type:'catalog_island',data:isl,x,y,w,h}});
        }});

        // Lane labels on left margin
        for (let lane=0; lane<=maxLane; lane++) {{
          const labelY = ISL_BASE_Y + lane*(ISL_ROW_H+ISL_GAP) + ISL_ROW_H/2;
          ctx.fillStyle = 'rgba(128,128,128,0.5)';
          ctx.font = '7px IBM Plex Mono,monospace';
          ctx.textBaseline = 'middle';
          ctx.textAlign = 'right';
          ctx.fillText(lane===0 ? 'GI' : 'GI+'+lane, PAD-2, labelY);
          ctx.textAlign = 'left';
        }}
      }}
    }}
    ctx.fillStyle='#7c2d92';ctx.font='9px IBM Plex Mono,monospace';
    ctx.textAlign='right';
    const PP_Y2 = GENOME_Y+GENOME_H+(visible.rgp?RGP_H+8:4);
    ctx.fillText('Phage',PAD-4,PP_Y2+PP_H/2+3);
  }}

  // Arrow helper — draws above genome
  function drawArrows(systems, baseY, colorFn, trackLabel) {{
    if (!systems.length) return 0;
    const levels={{}};
    const sorted=[...systems].filter(d=>inView(d.start,d.stop)).sort((a,b)=>a.start-b.start);
    if (!sorted.length) return 0;
    sorted.forEach(d => {{
      const x=bpToX(d.start), w=Math.max(4,bpWidth(d.stop-d.start));
      let lv=0;
      while(true){{const occ=levels[lv]||0;if(x>=occ){{levels[lv]=x+w+2;break;}}lv++;}}
      const y=baseY-lv*(ARROW_H+3);
      const c=colorFn(d);
      ctx.fillStyle=c+'cc';ctx.strokeStyle=c;ctx.lineWidth=1;
      const aw=Math.min(7,w*0.3);
      ctx.beginPath();
      if(w>10){{ctx.moveTo(x,y);ctx.lineTo(x+w-aw,y);ctx.lineTo(x+w,y+ARROW_H/2);
        ctx.lineTo(x+w-aw,y+ARROW_H);ctx.lineTo(x,y+ARROW_H);}}
      else ctx.rect(x,y,w,ARROW_H);
      ctx.closePath();ctx.fill();ctx.stroke();
      if(w>28){{ctx.fillStyle='#fff';ctx.font='8px IBM Plex Mono,monospace';ctx.textAlign='left';
        ctx.fillText((d.subtype||d.type||'').replace(/_/g,' ').substring(0,14),x+3,y+ARROW_H/2+3);}}
      ctx.strokeStyle=c+'44';ctx.lineWidth=1;ctx.setLineDash([2,3]);
      ctx.beginPath();ctx.moveTo(x+w/2,y+ARROW_H);ctx.lineTo(x+w/2,GENOME_Y);
      ctx.stroke();ctx.setLineDash([]);
      hits.push({{type:'defense',data:d,x,y,w,h:ARROW_H}});
    }});
    const rows=Object.keys(levels).length;
    ctx.fillStyle='#5a6380';ctx.font='9px IBM Plex Mono,monospace';ctx.textAlign='right';
    ctx.fillText(trackLabel,PAD-4,baseY-2);
    return rows;
  }}

  let curY=GENOME_Y-TRACK_GAP;
  if(visible.df && df.length>0) {{
    const rows=drawArrows(df,curY-ARROW_H,defColor,'DF');
    curY-=Math.min(rows||1,5)*(ARROW_H+3)+TRACK_GAP;
  }}
  if(visible.padloc && pl.length>0) {{
    const rows=drawArrows(pl,curY-ARROW_H,()=>'#1a5fa8','PL');
    curY-=Math.min(rows||1,5)*(ARROW_H+3)+TRACK_GAP;
  }}
  if(visible.dp && dp.length>0) {{
    drawArrows(dp,curY-ARROW_H,()=>'#8b4513','DP');
  }}

  canvasHits[id] = hits;

  // Zoom label
  const zl = document.getElementById('zl_'+id);
  if (zl) zl.textContent = vs.zoom===1?'1×':vs.zoom.toFixed(1)+'×';

  // Draw minimap
  if (minimapEl) drawMinimap(minimapEl, strain, gLen, maxLen, panClamped, visibleBp);
}}

// ── Minimap ──────────────────────────────────────────────────────────────────
function drawMinimap(mmEl, strain, gLen, maxLen, panStart, visBp) {{
  const data=STRAIN_DATA[strain];
  const W=mmEl.parentElement.clientWidth-4, H=28;
  const DPR=window.devicePixelRatio||1;
  mmEl.width=W*DPR; mmEl.height=H*DPR;
  mmEl.style.width=W+'px'; mmEl.style.height=H+'px';
  const ctx=mmEl.getContext('2d');
  ctx.scale(DPR,DPR);
  ctx.clearRect(0,0,W,H);

  const scale=W/maxLen;

  // Genome bar
  ctx.fillStyle='#d0d5e0';
  ctx.fillRect(0,8,gLen*scale,12);

  // Prophages
  if(visible.pp)(data.prophage||[]).forEach(p=>{{
    ctx.fillStyle=ppColor(p)+'99';
    ctx.fillRect(p.start*scale,8,Math.max(1,(p.stop-p.start)*scale),12);
  }});

  // RGPs
  if(visible.rgp)(data.rgps||[]).filter(r=>(r.hgt_evidence_count||0)>=MIN_EVIDENCE||(data.df||[]).some(d=>d.rgp_id===r.rgp_id)).forEach(r=>{{
    const ev=r.hgt_evidence_count||0;
    const hasSpotDef=(data.df||[]).some(d=>d.rgp_id===r.rgp_id&&['spot_13','spot_7','spot_48'].includes(d.spot_id));
    ctx.fillStyle=hasSpotDef?'#0e7a5a':ev>=3?'#c0392b':ev>=2?'#e67e22':ev>=1?'#2980b9':'#95a5a6';
    ctx.fillRect(r.start*scale,12,Math.max(1,(r.stop-r.start)*scale),4);
  }});

  // Defense systems
  if(visible.df)(data.df||[]).forEach(d=>{{
    const c=defColor(d);
    ctx.fillStyle=c;
    ctx.fillRect(d.start*scale,8,Math.max(2,(d.stop-d.start)*scale),12);
  }});

  // Viewport highlight
  const vpX=panStart*scale, vpW=visBp*scale;
  ctx.fillStyle='rgba(30,95,168,0.15)';
  ctx.fillRect(vpX,0,vpW,H);
  ctx.strokeStyle='#1a5fa8';ctx.lineWidth=1.5;
  ctx.strokeRect(vpX,1,vpW,H-2);
}}

// ── Zoom / Pan interaction ────────────────────────────────────────────────────
function setupInteraction(canvasEl, minimapEl, strain, maxLen) {{
  const id=canvasEl.id;
  if (!viewState[id]) viewState[id]={{zoom:1,pan:0}};

  // Mouse wheel zoom
  canvasEl.addEventListener('wheel', e=>{{
    e.preventDefault();
    const vs=viewState[id];
    const rect=canvasEl.getBoundingClientRect();
    const mx=e.clientX-rect.left;
    const PAD=56, trackW=rect.width-PAD-8;
    const bpAtCursor=vs.pan + (mx-PAD)/( (trackW)/(maxLen/vs.zoom) );
    const factor=e.deltaY<0?1.25:0.8;
    vs.zoom=Math.max(1,Math.min(50,vs.zoom*factor));
    const newVisBp=maxLen/vs.zoom;
    vs.pan=Math.max(0,Math.min(maxLen-newVisBp, bpAtCursor - newVisBp*(mx-PAD)/trackW));
    drawTrack(canvasEl,minimapEl,strain,maxLen);
  }},{{passive:false}});

  // Click-drag pan
  let dragging=false, lastX=0;
  canvasEl.addEventListener('mousedown', e=>{{
    dragging=true; lastX=e.clientX;
    canvasEl.classList.add('panning');
  }});
  window.addEventListener('mousemove', e=>{{
    if(!dragging) return;
    const vs=viewState[id];
    const rect=canvasEl.getBoundingClientRect();
    const trackW=rect.width-56-8;
    const bpPerPx=(maxLen/vs.zoom)/trackW;
    vs.pan=Math.max(0,Math.min(maxLen-maxLen/vs.zoom, vs.pan-(e.clientX-lastX)*bpPerPx));
    lastX=e.clientX;
    drawTrack(canvasEl,minimapEl,strain,maxLen);
  }});
  window.addEventListener('mouseup',()=>{{
    dragging=false;
    canvasEl.classList.remove('panning');
  }});

  // Minimap click to jump
  if (minimapEl) {{
    minimapEl.addEventListener('click', e=>{{
      const vs=viewState[id];
      const rect=minimapEl.getBoundingClientRect();
      const mx=e.clientX-rect.left;
      const frac=mx/rect.width;
      const newPan=frac*maxLen - (maxLen/vs.zoom)/2;
      vs.pan=Math.max(0,Math.min(maxLen-maxLen/vs.zoom,newPan));
      drawTrack(canvasEl,minimapEl,strain,maxLen);
    }});
  }}

  // Zoom buttons
  const zIn=document.getElementById('zi_'+id);
  const zOut=document.getElementById('zo_'+id);
  const zReset=document.getElementById('zr_'+id);
  if(zIn)    zIn.addEventListener('click',()=>{{const vs=viewState[id];const mid=vs.pan+(maxLen/vs.zoom)/2;vs.zoom=Math.min(50,vs.zoom*1.5);const nv=maxLen/vs.zoom;vs.pan=Math.max(0,Math.min(maxLen-nv,mid-nv/2));drawTrack(canvasEl,minimapEl,strain,maxLen);}});
  if(zOut)   zOut.addEventListener('click',()=>{{const vs=viewState[id];const mid=vs.pan+(maxLen/vs.zoom)/2;vs.zoom=Math.max(1,vs.zoom/1.5);const nv=maxLen/vs.zoom;vs.pan=Math.max(0,Math.min(maxLen-nv,mid-nv/2));drawTrack(canvasEl,minimapEl,strain,maxLen);}});
  if(zReset) zReset.addEventListener('click',()=>{{viewState[id]={{zoom:1,pan:0}};drawTrack(canvasEl,minimapEl,strain,maxLen);}});
}}

// ── Change analysis ──────────────────────────────────────────────────────────
function analyzeChanges(isolates) {{
  const perIsolate=isolates.map(s=>{{
    const d=STRAIN_DATA[s];
    return d ? new Set((d.df||[]).map(f=>f.subtype)) : new Set();
  }});
  const baseline=perIsolate[0];
  return isolates.map((s,i)=>{{
    if(i===0) return {{strain:s,gained:[],lost:[],stable:[...baseline]}};
    const cur=perIsolate[i];
    return {{strain:s,
      gained:[...cur].filter(x=>!baseline.has(x)),
      lost:[...baseline].filter(x=>!cur.has(x)),
      stable:[...cur].filter(x=>baseline.has(x))}};
  }});
}}

// ── Render patient ────────────────────────────────────────────────────────────
// Fetch strain data from per-strain JSON file (with cache)
async function fetchStrainData(strain) {{
  if (STRAIN_DATA[strain]) return STRAIN_DATA[strain];
  if (_fetchCache[strain]) return _fetchCache[strain];
  try {{
    const resp = await fetch(`${{DATA_DIR}}/${{strain}}.json`);
    if (!resp.ok) throw new Error(`HTTP ${{resp.status}}`);
    const data = await resp.json();
    STRAIN_DATA[strain] = data;
    _fetchCache[strain] = data;
    return data;
  }} catch(e) {{
    console.warn(`Could not load data for ${{strain}}:`, e);
    return {{}};
  }}
}}

async function redrawAll(pid) {{
  const isolates = PATIENT_GROUPS[pid] || [];

  // Show loading indicator
  isolates.forEach(strain => {{
    const cid = 'canvas_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const c = document.getElementById(cid);
    if (c) {{
      const ctx = c.getContext('2d');
      ctx.clearRect(0,0,c.width,c.height);
      ctx.fillStyle = 'rgba(128,128,128,0.3)';
      ctx.font = '11px IBM Plex Mono,monospace';
      ctx.fillText('Loading...', 10, 20);
    }}
  }});

  // Fetch all strains in parallel
  await Promise.all(isolates.map(s => fetchStrainData(s)));

  // Now draw with loaded data
  const maxLen = Math.max(...isolates.map(s => STRAIN_DATA[s]?.length || 5e6));
  isolates.forEach(strain => {{
    const cid = 'canvas_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const mid = 'mini_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const c = document.getElementById(cid);
    const m = document.getElementById(mid);
    if(c) drawTrack(c, m, strain, maxLen);
  }});
}}

function renderPatient(patientId) {{
  const isolates=PATIENT_GROUPS[patientId]||[];
  const maxLen=Math.max(...isolates.map(s=>STRAIN_DATA[s]?.length||5e6));
  const changes=analyzeChanges(isolates);

  // Patient summary
  const totalDF=new Set(isolates.flatMap(s=>(STRAIN_DATA[s]?.df||[]).map(d=>d.subtype))).size;
  const totalPP=isolates.reduce((n,s)=>n+(STRAIN_DATA[s]?.prophage||[]).length,0);
  const hasChanges=changes.some(c=>c.gained.length>0||c.lost.length>0);
  document.getElementById('patientSummary').innerHTML=`
    <div class="ps-item">Patient <span>${{patientId}}</span></div>
    <div class="ps-item">Isolates <span>${{isolates.length}}</span></div>
    <div class="ps-item">Unique DF subtypes <span>${{totalDF}}</span></div>
    <div class="ps-item">Prophage regions <span>${{totalPP}}</span></div>
    <div class="ps-item">Defense changes <span>${{hasChanges?'Yes':'None detected'}}</span></div>
  `;

  const list=document.getElementById('isolateList');
  list.innerHTML='';

  isolates.forEach((strain,idx)=>{{
    const data=STRAIN_DATA[strain];
    if(!data) return;
    const asm=data.asm||{{}};
    const q=asm.quality||'unknown';
    const ch=changes[idx];
    const compPct=asm.completeness?parseFloat(asm.completeness).toFixed(1)+'%':'—';
    const cid='canvas_'+strain.replace(/[^a-zA-Z0-9]/g,'_');
    const mid='mini_'+strain.replace(/[^a-zA-Z0-9]/g,'_');

    let badges='';
    if(idx===0) badges='<span class="change-badge stable">baseline</span>';
    else {{
      ch.gained.forEach(t=>badges+=`<span class="change-badge gained">+${{t}}</span>`);
      ch.lost.forEach(t=>badges+=`<span class="change-badge lost">-${{t}}</span>`);
      if(!ch.gained.length&&!ch.lost.length) badges='<span class="change-badge stable">no change</span>';
    }}

    const block=document.createElement('div');
    block.className='isolate-block';
    block.innerHTML=`
      <div class="isolate-header">
        <div class="isolate-name">${{strain}}</div>
        <span class="asm-pill pill-${{q}}">${{q.replace(/_/g,' ')}}</span>
        <div class="iso-stats">
          <div class="iso-stat">Contigs <span>${{asm.contigs||'—'}}</span></div>
          <div class="iso-stat">Complete <span>${{compPct}}</span></div>
          <div class="iso-stat">RGPs <span>${{(data.rgps||[]).length}}</span></div>
          <div class="iso-stat">DF <span>${{(data.df||[]).length}}</span></div>
          <div class="iso-stat">PADLOC <span>${{(data.padloc||[]).length}}</span></div>
          <div class="iso-stat">Prophage <span>${{(data.prophage||[]).length}}</span></div>
        </div>
        <div class="zoom-controls">
          <div class="zoom-btn" id="zo_${{cid}}" title="Zoom out">−</div>
          <div class="zoom-label" id="zl_${{cid}}">1×</div>
          <div class="zoom-btn" id="zi_${{cid}}" title="Zoom in">+</div>
          <div class="zoom-btn" id="zr_${{cid}}" title="Reset zoom" style="font-size:11px;width:36px">Reset</div>
        </div>
        <div style="margin-left:8px;display:flex;gap:4px;flex-wrap:wrap">${{badges}}</div>
      </div>
      <div class="minimap-wrap"><canvas id="${{mid}}" class="minimap" style="height:28px"></canvas></div>
      <div class="track-wrap"><canvas id="${{cid}}"></canvas></div>
    `;
    list.appendChild(block);

    requestAnimationFrame(()=>{{
      const c=document.getElementById(cid);
      const m=document.getElementById(mid);
      if(c){{
        drawTrack(c,m,strain,maxLen);
        setupInteraction(c,m,strain,maxLen);
      }}
    }});
  }});
}}

// ── Tooltip ──────────────────────────────────────────────────────────────────
const tooltip=document.getElementById('tooltip');
document.addEventListener('mousemove',e=>{{
  const canvas=e.target.closest('canvas:not(.minimap)');
  if(!canvas){{tooltip.style.display='none';return;}}
  const hits=canvasHits[canvas.id];
  if(!hits){{tooltip.style.display='none';return;}}
  const rect=canvas.getBoundingClientRect();
  const mx=e.clientX-rect.left,my=e.clientY-rect.top;
  let hit=null;
  for(const r of hits){{if(r.type==='genome')continue;if(mx>=r.x&&mx<=r.x+r.w&&my>=r.y&&my<=r.y+r.h){{hit=r;break;}}}}
  if(!hit)for(const r of hits){{if(r.type==='genome'&&mx>=r.x&&mx<=r.x+r.w&&my>=r.y&&my<=r.y+r.h){{hit=r;break;}}}}
  if(!hit){{tooltip.style.display='none';return;}}
  const d=hit.data;
  let html='';
  if(hit.type==='genome'){{
    const bp=Math.round(d.bpPerPx*(mx-56)+d.panClamped)||0;
    const bpFmt=bp>=1e6?(bp/1e6).toFixed(3)+' Mb':bp>=1e3?(bp/1e3).toFixed(1)+' kb':bp+' bp';
    html=`<div class="tt-title">Chromosome — ${{d.strain}}</div>
          <div class="tt-row">Position: <span>${{bpFmt}}</span></div>
          <div class="tt-row">Genome %: <span>${{(100*bp/d.length).toFixed(1)}}%</span></div>
          <div class="tt-row">Total length: <span>${{(d.length/1e6).toFixed(3)}} Mb</span></div>`;
    const nearby=[...(STRAIN_DATA[d.strain]?.df||[]),...(STRAIN_DATA[d.strain]?.prophage||[])]
      .filter(f=>bp>=f.start-5000&&bp<=f.stop+5000).slice(0,4);
    if(nearby.length){{html+=`<div class="tt-row" style="margin-top:5px;border-top:1px solid #eee;padding-top:4px">Nearby:</div>`;
      nearby.forEach(f=>html+=`<div class="tt-row" style="padding-left:8px">• <span>${{f.type||f.prophage_id}} (${{f.tool||'Prophage'}})</span></div>`);}}
  }} else if(hit.type==='catalog_island'){{
    const d=hit.data;
    html=`<div class="tt-title">${{d.island_id||''}}<span class="tt-badge" style="background:#1D9E7522;color:#0e7a5a">${{d.status==='nested_child'?'Nested GI':'Genomic Island'}}</span></div>
          <div class="tt-row">Group: <span class="tt-val">${{d.group_id||'—'}}</span></div>
          <div class="tt-row">Cargo: <span class="tt-val">${{d.dominant_cargo||'unknown'}}</span></div>
          <div class="tt-row">Depth: <span class="tt-val">${{d.depth===0?'Top-level':'Nested (depth '+d.depth+')'}}</span></div>
          <div class="tt-row">Evidence: <span class="tt-val">${{d.n_evidence}} lines</span></div>
          <div class="tt-row">Age: <span class="tt-val">${{d.age||'—'}}</span></div>
          <div class="tt-row">Coords: <span class="tt-val">${{d.start.toLocaleString()}}–${{d.end.toLocaleString()}}</span></div>
          <div class="tt-row">Genes: <span class="tt-val">${{d.n_genes}} (${{d.dominant_cargo}})</span></div>
          ${{d.status==='nested_child'&&d.parent?'<div class="tt-row">Parent: <span class="tt-val">'+d.parent+'</span></div>':''}}`;
  }} else if(hit.type==='fixed_defense'){{
    html=`<div class="tt-title">${{d.subtype}}<span class="tt-badge" style="background:#c0392b22;color:#c0392b">Fixed defense</span></div>
          <div class="tt-row">Age: <span class="tt-val">${{d.age_estimate||'—'}}</span></div>
          <div class="tt-row">Island: <span class="tt-val">${{(d.island_start||0).toLocaleString()}}–${{(d.island_end||0).toLocaleString()}}</span></div>
          <div class="tt-row">Confidence: <span class="tt-val">${{d.confidence||'—'}}</span></div>`;
  }} else if(hit.type==='prophage'){{
    const sz=d.stop-d.start;
    html=`<div class="tt-title">${{d.prophage_id}}<span class="tt-badge tt-pp">Prophage</span></div>
          <div class="tt-row">Confidence: <span>${{d.confidence}}</span></div>
          <div class="tt-row">Max pred: <span>${{d.max_pred}}</span></div>
          <div class="tt-row">Phage homology: <span>${{d.max_phage}}%</span></div>
          <div class="tt-row">Genes: <span>${{d.n_genes}}</span></div>
          <div class="tt-row">Size: <span>${{sz>=1000?(sz/1000).toFixed(1)+' kb':sz+' bp'}}</span></div>
          <div class="tt-row">Coords: <span>${{d.start.toLocaleString()}}–${{d.stop.toLocaleString()}}</span></div>`;
  }} else if(hit.type==='defense'){{
    const tc=d.tool==='DefenseFinder'?'tt-df':d.tool==='PADLOC'?'tt-padloc':'tt-dp';
    html=`<div class="tt-title">${{d.sys_id}}<span class="tt-badge ${{tc}}">${{d.tool}}</span></div>
          <div class="tt-row">Type: <span>${{d.type}}</span></div>
          <div class="tt-row">Subtype: <span>${{(d.subtype||'').substring(0,40)}}</span></div>`;
    if(d.tool==='DefenseFinder'){{
      const sb=d.spot_id&&KEY_SPOTS.has(d.spot_id)?`<span style="color:${{SPOT_COLORS[d.spot_id]}};font-weight:600">${{d.spot_id}}</span>`:(d.spot_id||'—');
      html+=`<div class="tt-row">In RGP: <span>${{d.in_rgp}}</span></div><div class="tt-row">Spot: ${{sb}}</div>`;
    }}
    if(d.tool==='PADLOC') html+=`<div class="tt-row">Proteins: <span>${{(d.proteins||'').substring(0,50)}}</span></div>`;
    if(d.tool==='DefensePredictor'){{html+=`<div class="tt-row">Tier: <span>${{d.tiers}}</span></div><div class="tt-row">Pfam: <span>${{d.pfams}}</span></div><div class="tt-row">Prob: <span>${{d.mean_prob}}</span></div>`;}}
    const sz=d.stop-d.start;
    html+=`<div class="tt-row">Size: <span>${{sz>=1000?(sz/1000).toFixed(1)+' kb':sz+' bp'}}</span></div>`;
  }} else if(hit.type==='rgp'){{
    const sz=d.stop-d.start;
    html=`<div class="tt-title">${{d.rgp_id}}</div>
          <div class="tt-row">Size: <span>${{sz>=1000?(sz/1000).toFixed(1)+' kb':sz+' bp'}}</span></div>
          <div class="tt-row">Score: <span>${{d.score}}</span></div>
          <div class="tt-row">Coords: <span>${{d.start.toLocaleString()}}–${{d.stop.toLocaleString()}}</span></div>`;
  }}
  tooltip.innerHTML=html;
  tooltip.style.display='block';
  tooltip.style.left=(e.clientX+16)+'px';
  tooltip.style.top=(e.clientY-10)+'px';
  const tr=tooltip.getBoundingClientRect();
  if(tr.right>window.innerWidth-10)  tooltip.style.left=(e.clientX-tr.width-10)+'px';
  if(tr.bottom>window.innerHeight-10) tooltip.style.top=(e.clientY-tr.height-10)+'px';
}});
document.addEventListener('mouseleave',()=>tooltip.style.display='none');

// ── Populate dropdown — all patients, alphabetical ───────────────────────────
const sel=document.getElementById('patientSelect');
Object.keys(PATIENT_GROUPS).sort((a,b)=>{{
  const na=parseInt(a.replace('GD',''))||0, nb=parseInt(b.replace('GD',''))||0;
  return na-nb;
}}).forEach(pid=>{{
  const opt=document.createElement('option');
  opt.value=pid;
  const n=PATIENT_GROUPS[pid].length;
  opt.textContent=n>1?`${{pid}} (${{n}} isolates)`:pid;
  sel.appendChild(opt);
}});
sel.addEventListener('change',()=>renderPatient(sel.value));
window.addEventListener('resize',()=>{{const pid=sel.value;if(pid)redrawAll(pid);}});

// Default to GD233
const def=PATIENT_GROUPS['GD233']?'GD233':Object.keys(PATIENT_GROUPS).sort((a,b)=>parseInt(a.replace('GD',''))-parseInt(b.replace('GD','')))[0];
sel.value=def;
renderPatient(def);
</script>
</body>
</html>"""
    return html, strain_data


def main():
    args = parse_args()
    print("\nLoading data...")
    genome_lengths = load_genome_lengths(args.json_dir)
    genome_stats   = load_genome_stats(args.genome_stats)
    rgps           = load_rgps(args.rgp_file)
    df             = load_defensefinder(args.intersection)
    pl             = load_padloc(args.padloc)
    dp             = load_defense_predictor(args.defense_predictor)
    prophages      = load_prophages(args.prophage_dir)

    # Load optional massiliense data
    if args.rgp_file_mass and args.spots_file_mass:
        print("  Loading massiliense RGPs...")
        for s, v in load_rgps(args.rgp_file_mass).items():
            if s not in rgps: rgps[s] = []
            rgps[s] += v
    if args.intersection_mass:
        print("  Loading massiliense defense intersection...")
        for s, v in load_defensefinder(args.intersection_mass).items():
            if s not in df: df[s] = []
            df[s] += v

    # Load optional bolletii data
    if args.rgp_file_boll and args.spots_file_boll:
        print("  Loading bolletii RGPs...")
        for s, v in load_rgps(args.rgp_file_boll).items():
            if s not in rgps: rgps[s] = []
            rgps[s] += v
    if args.intersection_boll:
        print("  Loading bolletii defense intersection...")
        for s, v in load_defensefinder(args.intersection_boll).items():
            if s not in df: df[s] = []
            df[s] += v

    print("\nGrouping patients...")
    # Use union of ALL data sources so strains without JSON tracks are included
    all_known_strains = set(genome_lengths.keys())
    all_known_strains |= set(rgps.keys())
    all_known_strains |= set(df.keys())
    all_known_strains |= set(pl.keys())
    all_known_strains |= set(dp.keys())
    all_known_strains |= set(prophages.keys())
    all_known_strains |= set(genome_stats.keys())
    patient_groups = group_patients(list(all_known_strains))
    multi = sum(1 for v in patient_groups.values() if len(v)>1)
    print(f"  {len(patient_groups)} total patients ({multi} with multiple isolates)")

    trna_data = {}
    if args.trna_proximity:
        print("  Loading abscessus tRNA proximity data...")
        trna_data.update(load_trna_proximity(args.trna_proximity))

    print("  Loading scanner islands and fixed defense data...")
    scanner_islands = load_scanner_islands(args.scanner_islands)
    denovo_defense  = load_denovo_defense(args.denovo_defense)
    island_catalog  = load_island_catalog(getattr(args,"island_catalog",None))
    if args.trna_proximity_mass:
        print("  Loading massiliense tRNA proximity data...")
        trna_data.update(load_trna_proximity(args.trna_proximity_mass))
    if args.trna_proximity_boll:
        print("  Loading bolletii tRNA proximity data...")
        trna_data.update(load_trna_proximity(args.trna_proximity_boll))
    if trna_data:
        print(f"  Total tRNA proximity entries: {len(trna_data)}")

    print("\nBuilding HTML...")
    # Write per-strain JSON files for lazy loading
    data_dir = os.path.join(os.path.dirname(args.outfile), "viewer_data")
    os.makedirs(data_dir, exist_ok=True)

    html, strain_data = build_html(genome_lengths, rgps, df, pl, dp, prophages, genome_stats, patient_groups, trna_data, args.min_evidence, scanner_islands, denovo_defense, island_catalog)

    # Write full strain_data to per-strain JSON files for lazy loading
    import json as _json
    n_written = 0
    print(f"  Writing per-strain data to {data_dir}/...")
    for strain, sdata in strain_data.items():
        strain_file = os.path.join(data_dir, f"{strain}.json")
        with open(strain_file, "w") as f:
            _json.dump(sdata, f)
        n_written += 1
    print(f"  Written {n_written} complete strain JSON files")

    with open(args.outfile, "w") as f:
        f.write(html)
    print(f"\nDone → {args.outfile}  ({len(html)//1024} KB)")


if __name__ == "__main__":
    main()
