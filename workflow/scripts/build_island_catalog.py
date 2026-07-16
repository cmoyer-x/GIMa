"""
build_island_catalog.py

Resolves overlapping genomic islands, assigns stable IDs, and builds
the final catalog with nested island relationships. Output feeds
directly into the patient viewer for nested track rendering.

Overlap resolution rules:
  1. Near-identical (>80% overlap) -> deduplicate, keep higher evidence
  2. Containment (one fully inside other) -> nested parent/child
  3. Partial overlap (<80%, >0%) -> trim at midpoint if <20% of smaller
                                    island, otherwise merge

Output files:
  catalog_groups.tsv   — one row per unique island locus (Mabs_GI_XXX)
  catalog_strains.tsv  — one row per strain-island (GD01_GI_XXX)
  catalog_nesting.tsv  — parent-child relationships

Usage:
    python build_island_catalog.py \\
        --islands  island_predictions/all_islands_combined.tsv \\
        --defense  denovo_defense_intersection.tsv \\
        --amr      amrfinder_results \\
        --gff_dir  all_gffs \\
        --out_dir  island_catalog \\
        --min_ev   3
"""

import os, re, csv, json, argparse
from collections import defaultdict, Counter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--islands",  required=True,
                   help="all_islands_combined.tsv from GIMa")
    p.add_argument("--defense",  default=None,
                   help="denovo_defense_intersection.tsv")
    p.add_argument("--amr_dir",  default=None,
                   help="Directory of AMRFinderPlus per-strain TSVs")
    p.add_argument("--gff_dir",  required=True,
                   help="Directory of Prokka GFF files")
    p.add_argument("--out_dir",  default="island_catalog")
    p.add_argument("--min_ev",   type=int, default=3,
                   help="Minimum evidence lines (default 3)")
    return p.parse_args()


def canonical(s):
    return re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=re.IGNORECASE)


# ── Loading ───────────────────────────────────────────────────────────────────
def load_islands(path, min_ev):
    islands_by_strain = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ev = int(row.get("n_evidence", 0))
            if ev < min_ev: continue
            conf = row.get("confidence","")
            if conf not in ("high","moderate"): continue
            try:
                islands_by_strain[row["strain"]].append({
                    "strain":       row["strain"],
                    "contig":       row["contig"],
                    "start":        int(row["start"]),
                    "end":          int(row["end"]),
                    "length":       int(row["length"]),
                    "rgp_seed":     row.get("rgp_seed",""),
                    "n_evidence":   ev,
                    "evidence":     row.get("evidence",""),
                    "confidence":   conf,
                    "age_estimate": row.get("age_estimate",""),
                    "cai_ratio":    float(row["cai_ratio"]) if row.get("cai_ratio") else 1.0,
                    "trna_flanked": row.get("trna_flanked","No"),
                    "trna_product": row.get("trna_product",""),
                    "has_dr":       row.get("has_dr","No"),
                    "dr_seq":       row.get("dr_seq",""),
                    "mob_types":    row.get("mob_types",""),
                    "isl_gc":       float(row["isl_gc"]) if row.get("isl_gc") else 0,
                    "is_denovo":    not bool(row.get("rgp_seed","")),
                })
            except: continue
    total = sum(len(v) for v in islands_by_strain.values())
    print(f"  Loaded {total:,} islands (ev>={min_ev}) across "
          f"{len(islands_by_strain)} strains")
    return dict(islands_by_strain)


def load_defense_map(path):
    """Returns dict: (strain, approx_coord) -> list of defense subtypes"""
    if not path or not os.path.exists(path): return {}
    dm = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("in_denovo_island") != "Yes": continue
            try:
                strain = row["strain"]
                mid = (int(row["island_start"]) + int(row["island_end"])) // 2
                dm[(strain, mid // 5000)].append(row["subtype"])
            except: continue
    return dm


def load_gff_genes(gff_dir):
    """Returns dict: strain -> list of gene dicts with functional categories."""
    genes_by_strain = {}
    MOBILITY = [r"integrase",r"transposase",r"recombinase",r"relaxase",r"conjugal"]
    TA       = [r"vap[bc]",r"mazef",r"antitoxin",r"toxin.antitoxin",r"higb",r"relbe"]
    METAL    = [r"mercury",r"arsenic",r"copper.resist",r"zinc.resist",r"mer[abc]"]
    EFFLUX   = [r"mmpl",r"mmps",r"efflux",r"multidrug.resist"]
    DEFENSE  = [r"defense",r"anti.phage",r"restriction",r"modification"]
    PHAGE    = [r"capsid",r"tail",r"baseplate",r"terminase",r"phage"]
    REG      = [r"whib",r"sigma.factor",r"transcriptional.regul"]

    def classify(product, gene=""):
        text = f"{product} {gene}".lower()
        for cat, patterns in [
            ("mobility",  MOBILITY), ("ta_system",  TA),
            ("metal",     METAL),    ("efflux",     EFFLUX),
            ("defense",   DEFENSE),  ("phage",      PHAGE),
            ("regulatory",REG),
        ]:
            if any(re.search(p, text, re.I) for p in patterns):
                return cat
        if "hypothetical" in text or "unknown" in text:
            return "hypothetical"
        return "other"

    for fname in os.listdir(gff_dir):
        if not fname.endswith(".gff"): continue
        strain = fname.replace("_prokka.gff","").replace(".gff","")
        genes = []
        try:
            with open(os.path.join(gff_dir, fname)) as f:
                for line in f:
                    if "##FASTA" in line: break
                    if line.startswith("#") or "\tCDS\t" not in line: continue
                    parts = line.strip().split("\t")
                    if len(parts) < 9: continue
                    attrs = {}
                    for a in parts[8].split(";"):
                        if "=" in a:
                            k,v = a.split("=",1)
                            attrs[k.strip()] = v.strip()
                    product = attrs.get("product","").replace("%2C",",")
                    gene    = attrs.get("gene","")
                    genes.append({
                        "contig":    parts[0].replace("gnl|XXX|",""),
                        "start":     int(parts[3])-1,
                        "end":       int(parts[4]),
                        "locus_tag": attrs.get("locus_tag",""),
                        "product":   product,
                        "gene":      gene,
                        "category":  classify(product, gene),
                    })
        except: continue
        genes_by_strain[strain] = genes
    print(f"  Loaded GFF genes for {len(genes_by_strain)} strains")
    return genes_by_strain


# ── Overlap resolution ────────────────────────────────────────────────────────
def resolve_overlaps(islands):
    """
    Resolve overlapping islands within a single strain.
    Returns list of islands with overlap_status and nesting_depth fields.

    Rules:
      1. Near-identical (>80% reciprocal overlap) -> deduplicate
      2. Containment -> nested (parent keeps children list)
      3. Partial overlap (trim or merge)
    """
    if not islands: return []

    # Sort by start, then by length descending (larger islands first)
    islands = sorted(islands, key=lambda x: (x["contig"], x["start"], -x["length"]))

    # Add tracking fields
    for isl in islands:
        isl["overlap_status"] = "unique"
        isl["nesting_depth"]  = 0
        isl["parent_id"]      = None
        isl["children"]       = []
        isl["_keep"]          = True

    n = len(islands)
    for i in range(n):
        if not islands[i]["_keep"]: continue
        for j in range(i+1, n):
            if not islands[j]["_keep"]: continue
            a, b = islands[i], islands[j]

            # Must be on same contig
            if a["contig"] != b["contig"]: continue

            # No overlap
            if a["end"] <= b["start"] or b["end"] <= a["start"]: continue

            # Calculate overlap
            ov_start  = max(a["start"], b["start"])
            ov_end    = min(a["end"],   b["end"])
            ov_len    = ov_end - ov_start
            a_len     = a["end"] - a["start"]
            b_len     = b["end"] - b["start"]
            ov_pct_a  = ov_len / a_len if a_len > 0 else 0
            ov_pct_b  = ov_len / b_len if b_len > 0 else 0

            # Rule 1: Near-identical (>80% of both) -> deduplicate
            if ov_pct_a > 0.80 and ov_pct_b > 0.80:
                # Keep higher evidence, or guided over de novo
                keep_a = (a["n_evidence"] > b["n_evidence"] or
                          (a["n_evidence"] == b["n_evidence"] and
                           bool(a["rgp_seed"]) and not bool(b["rgp_seed"])))
                if keep_a:
                    b["_keep"] = False
                    b["overlap_status"] = "deduplicated"
                else:
                    a["_keep"] = False
                    a["overlap_status"] = "deduplicated"
                continue

            # Rule 2: Containment — b fully inside a
            if ov_pct_b > 0.90 and ov_pct_a < 0.90:
                b["overlap_status"] = "nested_child"
                b["parent_id"]      = id(a)
                a["overlap_status"] = "nested_parent"
                a["children"].append(id(b))
                continue

            # Rule 2b: a fully inside b
            if ov_pct_a > 0.90 and ov_pct_b < 0.90:
                a["overlap_status"] = "nested_child"
                a["parent_id"]      = id(b)
                b["overlap_status"] = "nested_parent"
                b["children"].append(id(a))
                continue

            # Rule 3: Partial overlap
            if ov_pct_a < 0.20 or ov_pct_b < 0.20:
                # Small overlap — trim at midpoint
                midpoint = (ov_start + ov_end) // 2
                if a["start"] < b["start"]:
                    a["end"]    = midpoint
                    b["start"]  = midpoint
                    a["length"] = a["end"] - a["start"]
                    b["length"] = b["end"] - b["start"]
                else:
                    b["end"]    = midpoint
                    a["start"]  = midpoint
                    b["length"] = b["end"] - b["start"]
                    a["length"] = a["end"] - a["start"]
                a["overlap_status"] = "trimmed"
                b["overlap_status"] = "trimmed"
            else:
                # Large partial overlap — merge into the one with more evidence
                if a["n_evidence"] >= b["n_evidence"]:
                    a["start"]  = min(a["start"], b["start"])
                    a["end"]    = max(a["end"],   b["end"])
                    a["length"] = a["end"] - a["start"]
                    a["n_evidence"] = max(a["n_evidence"], b["n_evidence"])
                    a["overlap_status"] = "merged"
                    b["_keep"]  = False
                    b["overlap_status"] = "merged_into"
                else:
                    b["start"]  = min(a["start"], b["start"])
                    b["end"]    = max(a["end"],   b["end"])
                    b["length"] = b["end"] - b["start"]
                    b["n_evidence"] = max(a["n_evidence"], b["n_evidence"])
                    b["overlap_status"] = "merged"
                    a["_keep"]  = False
                    a["overlap_status"] = "merged_into"

    # Assign nesting depths (BFS from top-level parents)
    id_map = {id(isl): isl for isl in islands}

    def assign_depth(isl, depth):
        isl["nesting_depth"] = depth
        for child_id in isl.get("children",[]):
            child = id_map.get(child_id)
            if child: assign_depth(child, depth+1)

    for isl in islands:
        if isl["_keep"] and isl["nesting_depth"] == 0 and not isl["parent_id"]:
            assign_depth(isl, 0)

    # Return only kept islands, sorted
    result = [isl for isl in islands if isl["_keep"]]
    result.sort(key=lambda x: (x["contig"], x["start"]))
    return result


# ── Cargo annotation ──────────────────────────────────────────────────────────
def annotate_cargo(island, genes_by_strain):
    """Count functional gene categories within island boundaries."""
    strain   = island["strain"]
    genes    = genes_by_strain.get(strain) or genes_by_strain.get(canonical(strain)) or []
    contig   = island["contig"]
    isl_start = island["start"]
    isl_end   = island["end"]

    counts = Counter()
    island_genes = []
    for g in genes:
        gc = re.sub(r'_\d+$','',g["contig"])
        ic = re.sub(r'_\d+$','',contig)
        if gc != ic and g["contig"] != contig: continue
        if g["start"] >= isl_start and g["end"] <= isl_end:
            counts[g["category"]] += 1
            island_genes.append(g)

    total = len(island_genes)
    n_hypo = counts["hypothetical"]

    # Determine dominant cargo, prioritizing informative categories.
    # Regulatory genes are near-ubiquitous, so they only win when no
    # informative cargo (defense/phage/mobility/metal/efflux/TA) is present.
    INFORMATIVE = ["defense", "phage", "mobility", "metal", "efflux", "ta_system"]
    informative_counts = {k: counts[k] for k in INFORMATIVE if counts.get(k, 0) > 0}
    if informative_counts:
        dominant = max(informative_counts, key=informative_counts.get)
    elif counts.get("regulatory", 0) > 0:
        dominant = "regulatory"
    elif n_hypo > total * 0.5:
        dominant = "hypothetical"
    else:
        dominant = "unknown"
    if dominant == "ta_system":
        dominant = "ta"

    return {
        "n_genes":        total,
        "n_hypothetical": n_hypo,
        "pct_hypothetical": round(100*n_hypo/total,1) if total > 0 else 0,
        "dominant_cargo": dominant,
        "n_mobility":     counts["mobility"],
        "n_defense":      counts["defense"],
        "n_ta":           counts["ta_system"],
        "n_metal":        counts["metal"],
        "n_efflux":       counts["efflux"],
        "n_phage":        counts["phage"],
        "n_regulatory":   counts["regulatory"],
        "n_other":        counts["other"],
    }


# ── Group assignment ──────────────────────────────────────────────────────────
def assign_groups(all_islands_flat):
    """
    Group islands across strains by approximate locus.
    Islands within 10kb of each other on the same relative chromosomal
    position with the same dominant cargo are the same insertion event.
    Returns dict: island -> group_id
    """
    # Sort by dominant_cargo + approximate coordinate bin
    groups = defaultdict(list)
    for isl in all_islands_flat:
        coord_bin = isl["start"] // 10000
        key = f"{isl['dominant_cargo']}_{coord_bin}"
        groups[key].append(isl)

    # Assign Mabs_GI_XXX IDs ordered by group size (most prevalent first)
    sorted_groups = sorted(groups.items(), key=lambda x: -len(x[1]))
    group_id_map  = {}
    for gi, (key, members) in enumerate(sorted_groups, 1):
        gid = f"Mabs_GI_{gi:03d}"
        for isl in members:
            group_id_map[id(isl)] = gid
    return group_id_map, {f"Mabs_GI_{gi:03d}": members
                          for gi, (_, members) in enumerate(sorted_groups, 1)}


# ── ID assignment ─────────────────────────────────────────────────────────────
def assign_strain_ids(resolved_by_strain):
    """Assign GD01_GI_001 style IDs ordered by chromosomal position."""
    strain_id_map = {}
    for strain, islands in resolved_by_strain.items():
        sorted_isls = sorted(islands, key=lambda x: (x["contig"], x["start"]))
        for gi, isl in enumerate(sorted_isls, 1):
            strain_id_map[id(isl)] = f"{strain}_GI_{gi:03d}"
    return strain_id_map


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("\nLoading data...")
    islands_by_strain = load_islands(args.islands, args.min_ev)
    defense_map       = load_defense_map(args.defense)
    genes_by_strain   = load_gff_genes(args.gff_dir)

    print("\nResolving overlaps per strain...")
    resolved_by_strain = {}
    overlap_stats = Counter()

    for strain, islands in islands_by_strain.items():
        resolved = resolve_overlaps(islands)
        resolved_by_strain[strain] = resolved
        for isl in resolved:
            overlap_stats[isl["overlap_status"]] += 1

    print(f"  Overlap resolution summary:")
    for status, n in overlap_stats.most_common():
        print(f"    {status:20}: {n:,}")

    print("\nAnnotating cargo...")
    for strain, islands in resolved_by_strain.items():
        for isl in islands:
            cargo = annotate_cargo(isl, genes_by_strain)
            isl.update(cargo)

    print("\nAssigning IDs...")
    all_flat = [isl for isls in resolved_by_strain.values() for isl in isls]
    strain_id_map, group_members = assign_groups(all_flat)
    island_id_map = assign_strain_ids(resolved_by_strain)

    for isl in all_flat:
        isl["island_id"] = island_id_map.get(id(isl),"")
        isl["group_id"]  = strain_id_map.get(id(isl),"")

    # ── Write catalog_strains.tsv ─────────────────────────────────────────────
    strain_fields = [
        "island_id","group_id","strain","contig","start","end","length",
        "nesting_depth","overlap_status","parent_id_ref",
        "n_evidence","confidence","evidence","age_estimate","cai_ratio",
        "trna_flanked","trna_product","has_dr","dr_seq","mob_types",
        "isl_gc","rgp_seed","is_denovo",
        "dominant_cargo","n_genes","pct_hypothetical",
        "n_mobility","n_defense","n_ta","n_metal",
        "n_efflux","n_phage","n_regulatory",
    ]

    # Build id -> island_id lookup for parent references
    id_to_strain_id = {id(isl): isl["island_id"] for isl in all_flat}

    strain_path = os.path.join(args.out_dir, "catalog_strains.tsv")
    with open(strain_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=strain_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for strain in sorted(resolved_by_strain.keys()):
            for isl in resolved_by_strain[strain]:
                row = dict(isl)
                row["parent_id_ref"] = id_to_strain_id.get(isl.get("parent_id"),"")
                row["is_denovo"]     = "Yes" if isl.get("is_denovo") else "No"
                writer.writerow(row)
    print(f"  Written: {strain_path}")

    # ── Write catalog_groups.tsv ──────────────────────────────────────────────
    group_fields = [
        "group_id","n_strains","mean_length_kb","dominant_cargo",
        "age_distribution","mean_cai","n_with_defense","n_with_dr",
        "n_with_trna","n_nested_parents","representative_strain",
        "representative_start","representative_end",
    ]
    group_path = os.path.join(args.out_dir, "catalog_groups.tsv")
    with open(group_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=group_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for gid, members in sorted(group_members.items()):
            strains = list(set(m["strain"] for m in members))
            ages    = Counter(m["age_estimate"] for m in members)
            rep     = max(members, key=lambda x: x["n_evidence"])
            writer.writerow({
                "group_id":           gid,
                "n_strains":          len(strains),
                "mean_length_kb":     round(sum(m["length"] for m in members)/
                                            len(members)/1000, 1),
                "dominant_cargo":     rep["dominant_cargo"],
                "age_distribution":   "; ".join(f"{k}:{v}"
                                                for k,v in ages.most_common()),
                "mean_cai":           round(sum(m["cai_ratio"] for m in members)/
                                            len(members), 3),
                "n_with_defense":     sum(1 for m in members if m["n_defense"]>0),
                "n_with_dr":          sum(1 for m in members if m["has_dr"]=="Yes"),
                "n_with_trna":        sum(1 for m in members
                                          if m["trna_flanked"]=="Yes"),
                "n_nested_parents":   sum(1 for m in members
                                          if m["overlap_status"]=="nested_parent"),
                "representative_strain": rep["strain"],
                "representative_start":  rep["start"],
                "representative_end":    rep["end"],
            })
    print(f"  Written: {group_path}")

    # ── Write nesting TSV ─────────────────────────────────────────────────────
    nest_path = os.path.join(args.out_dir, "catalog_nesting.tsv")
    with open(nest_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["parent_id","child_id","strain","parent_start",
                         "parent_end","child_start","child_end",
                         "child_dominant_cargo","child_n_evidence"])
        for isl in all_flat:
            if isl["overlap_status"] == "nested_child" and isl["parent_id"]:
                parent = next((x for x in all_flat if id(x)==isl["parent_id"]),None)
                if parent:
                    writer.writerow([
                        parent["island_id"], isl["island_id"],
                        isl["strain"],
                        parent["start"], parent["end"],
                        isl["start"],   isl["end"],
                        isl["dominant_cargo"], isl["n_evidence"],
                    ])
    print(f"  Written: {nest_path}")

    # ── Write viewer JSON ─────────────────────────────────────────────────────
    # Per-strain JSON for the patient viewer
    viewer_data = {}
    for strain, islands in resolved_by_strain.items():
        viewer_data[strain] = []
        for isl in sorted(islands, key=lambda x: x["start"]):
            viewer_data[strain].append({
                "island_id":     isl["island_id"],
                "group_id":      isl["group_id"],
                "start":         isl["start"],
                "end":           isl["end"],
                "length":        isl["length"],
                "depth":         isl["nesting_depth"],
                "status":        isl["overlap_status"],
                "parent":        id_to_strain_id.get(isl.get("parent_id"),""),
                "n_evidence":    isl["n_evidence"],
                "evidence":      isl["evidence"],
                "confidence":    isl["confidence"],
                "age":           isl["age_estimate"],
                "cai":           isl["cai_ratio"],
                "dominant_cargo":isl["dominant_cargo"],
                "n_genes":       isl["n_genes"],
                "n_defense":     isl["n_defense"],
                "n_mobility":    isl["n_mobility"],
                "n_ta":          isl["n_ta"],
                "trna_flanked":  isl["trna_flanked"],
                "has_dr":        isl["has_dr"],
                "mob_types":     isl["mob_types"],
                "is_denovo":     isl.get("is_denovo", False),
            })

    viewer_json_path = os.path.join(args.out_dir, "island_viewer_data.json")
    with open(viewer_json_path, "w") as f:
        json.dump(viewer_data, f)
    print(f"  Written: {viewer_json_path}")

    # ── Summary ───────────────────────────────────────────────────────────────
    total_kept   = len(all_flat)
    n_top_level  = sum(1 for x in all_flat if x["nesting_depth"]==0)
    n_nested     = sum(1 for x in all_flat if x["nesting_depth"]>0)
    n_parents    = sum(1 for x in all_flat if x["overlap_status"]=="nested_parent")
    cargo_dist   = Counter(x["dominant_cargo"] for x in all_flat)
    age_dist     = Counter(x["age_estimate"]    for x in all_flat)

    print(f"\n{'='*60}")
    print(f"  Total islands in catalog:    {total_kept:,}")
    print(f"  Top-level islands:           {n_top_level:,}")
    print(f"  Nested child islands:        {n_nested:,}")
    print(f"  Islands with nested children:{n_parents:,}")
    print(f"  Unique groups (Mabs_GI_XXX): {len(group_members):,}")
    print(f"  Strains represented:         {len(resolved_by_strain):,}")
    print(f"\nDominant cargo distribution:")
    for cargo, n in cargo_dist.most_common():
        pct = 100*n/total_kept
        print(f"  {cargo:15}: {n:5,} ({pct:.1f}%)")
    print(f"\nAge distribution:")
    for age in ["very_recent","recent","moderate","old"]:
        n = age_dist.get(age,0)
        print(f"  {age:15}: {n:5,} ({100*n/total_kept:.1f}%)" if total_kept else f"  {age:15}: {n:5,}")
    print(f"{'='*60}")
    print(f"\nViewer JSON ready: {viewer_json_path}")
    print(f"Load this in the updated patient viewer with --island_catalog flag")


if __name__ == "__main__":
    main()
