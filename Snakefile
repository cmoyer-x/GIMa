# =============================================================================
# MAbsIslandScanner — Snakemake Workflow
# =============================================================================
#
# Component 2: Genomic island scanning, catalog building, and viewer
#
# Assumes the following are already available (pre-computed):
#   - FASTA files per strain
#   - Prokka GFF files per strain
#   - PPanGGOLiN RGP output (regions_of_genomic_plasticity.tsv)
#   - DefenseFinder output (all_genomes_defense_finder_systems.tsv)
#   - Defense RGP intersection (defense_rgp_intersection.tsv)
#
# Usage:
#   snakemake --configfile config.yaml --cores 8
#
#   Dry run (show what would run):
#   snakemake --configfile config.yaml --cores 8 -n
#
#   With conda environments:
#   snakemake --configfile config.yaml --cores 8 --use-conda
#
#   On a SLURM cluster:
#   snakemake --configfile config.yaml --cluster "sbatch -c {threads} \
#       --mem={resources.mem_mb}M -t {resources.runtime}" \
#       --jobs 50 --use-conda
# =============================================================================

import os
from pathlib import Path


# ── Strain discovery ──────────────────────────────────────────────────────────
# Find all strains that have both FASTA and GFF files

def find_strains(fasta_dir, gff_dir):
    """Discover strains with both FASTA and GFF available."""
    import re

    def canonical(s):
        return re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=re.IGNORECASE)

    # FASTA strains
    fasta_strains = {}
    for f in Path(fasta_dir).glob("*.fasta"):
        strain = f.stem
        fasta_strains[strain] = f
    for f in Path(fasta_dir).glob("*.fa"):
        strain = f.stem
        fasta_strains[strain] = f
    for f in Path(fasta_dir).glob("*.fna"):
        strain = f.stem
        fasta_strains[strain] = f

    # GFF strains (try both with and without _prokka suffix)
    gff_strains = {}
    for f in Path(gff_dir).glob("*_prokka.gff"):
        strain = f.stem.replace("_prokka","")
        gff_strains[strain] = f
    for f in Path(gff_dir).glob("*.gff"):
        if "_prokka" not in f.name:
            strain = f.stem
            gff_strains[strain] = f

    # Intersection: strains with both FASTA and GFF
    # Handle suffix mismatches (GD08_hybrid FASTA + GD08 GFF)
    valid = []
    for fasta_strain, fasta_path in fasta_strains.items():
        # Direct match
        if fasta_strain in gff_strains:
            valid.append(fasta_strain)
            continue
        # Canonical match
        canon = canonical(fasta_strain)
        if canon in gff_strains:
            valid.append(fasta_strain)
            continue

    return sorted(set(valid))


STRAINS = find_strains(
    config["fasta_dir"],
    config["gff_dir"]
)

print(f"Discovered {len(STRAINS)} strains with FASTA + GFF")


# ── Helper: get GFF path for a strain ────────────────────────────────────────
def get_gff(strain):
    import re
    gff_dir = config["gff_dir"]
    canon   = re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', strain, flags=re.IGNORECASE)
    for name in [f"{strain}_prokka.gff", f"{canon}_prokka.gff",
                 f"{strain}.gff",        f"{canon}.gff"]:
        path = os.path.join(gff_dir, name)
        if os.path.exists(path):
            return path
    return os.path.join(gff_dir, f"{canon}_prokka.gff")


def get_fasta(strain):
    fasta_dir = config["fasta_dir"]
    for ext in [".fasta", ".fa", ".fna"]:
        path = os.path.join(fasta_dir, f"{strain}{ext}")
        if os.path.exists(path):
            return path
    return os.path.join(fasta_dir, f"{strain}.fasta")


# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        # Per-strain island predictions
        expand(
            os.path.join(config["out_dir"], "per_strain", "{strain}_islands.tsv"),
            strain=STRAINS
        ),
        # Combined island catalog
        os.path.join(config["out_dir"], "all_islands_combined.tsv"),
        # Fixed defense intersection
        os.path.join(config["out_dir"], "denovo_defense_intersection.tsv"),
        # Final island catalog
        os.path.join(config["out_dir"], "catalog", "catalog_strains.tsv"),
        os.path.join(config["out_dir"], "catalog", "catalog_groups.tsv"),
        os.path.join(config["out_dir"], "catalog", "island_viewer_data.json"),
        # Interactive viewer
        os.path.join(config["out_dir"], "patient_comparison_viewer.html"),


# ── Rule 1: Per-strain island scanning ───────────────────────────────────────
rule scan_genome:
    input:
        fasta = lambda wc: get_fasta(wc.strain),
        gff   = lambda wc: get_gff(wc.strain),
        rgp   = config["rgp_file"],
    output:
        tsv   = os.path.join(config["out_dir"], "per_strain", "{strain}_islands.tsv"),
    params:
        rgp_extend = config.get("rgp_extend", 5000),
        min_ev     = config.get("min_evidence", 2),
        window     = config.get("window_size", 8000),
        step       = config.get("step_size", 2000),
        script     = config["scanner_script"],
    threads: 1
    resources:
        mem_mb   = 2000,   # 2GB per genome — well above typical usage (~500MB)
        runtime  = 20,     # minutes — conservative for 5Mb genome
    conda:
        "envs/mabs_islands.yaml"
    log:
        os.path.join(config["out_dir"], "logs", "scan_{strain}.log"),
    shell:
        """
        python {params.script} \
            --fasta {input.fasta} \
            --gff   {input.gff} \
            --out   {output.tsv} \
            --rgp_guided {input.rgp} \
            --rgp_extend {params.rgp_extend} \
            --min_ev {params.min_ev} \
            --window {params.window} \
            --step   {params.step} \
            > {log} 2>&1
        """


# ── Rule 2: Combine per-strain results ───────────────────────────────────────
rule combine_islands:
    input:
        tsvs = expand(
            os.path.join(config["out_dir"], "per_strain", "{strain}_islands.tsv"),
            strain=STRAINS
        ),
    output:
        combined = os.path.join(config["out_dir"], "all_islands_combined.tsv"),
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 10,
    run:
        import csv, os

        all_rows = []
        header   = None

        for tsv in sorted(input.tsvs):
            if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
                continue
            with open(tsv) as f:
                reader = csv.DictReader(f, delimiter="\t")
                if header is None:
                    header = reader.fieldnames
                for row in reader:
                    all_rows.append(row)

        if not header:
            raise ValueError("No valid island TSV files found")

        with open(output.combined, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_rows)

        conf = {}
        for r in all_rows:
            c = r.get("confidence","")
            conf[c] = conf.get(c,0) + 1

        print(f"Combined {len(all_rows):,} islands from "
              f"{len([t for t in input.tsvs if os.path.getsize(t)>0])} strains")
        for c, n in sorted(conf.items()):
            print(f"  {c}: {n:,}")


# ── Rule 3: Fixed defense intersection ───────────────────────────────────────
rule intersect_defense:
    input:
        islands  = os.path.join(config["out_dir"], "all_islands_combined.tsv"),
        defense  = config["defense_systems"],
        gff_dir  = config["gff_dir"],
    output:
        tsv      = os.path.join(config["out_dir"], "denovo_defense_intersection.tsv"),
        report   = os.path.join(config["out_dir"], "denovo_defense_report.html"),
    params:
        min_conf = config.get("min_confidence", "moderate"),
        script   = config["intersect_script"],
    threads: 1
    resources:
        mem_mb  = 8000,   # GFF index can be large (~3M locus tags)
        runtime = 30,
    conda:
        "envs/mabs_islands.yaml"
    log:
        os.path.join(config["out_dir"], "logs", "intersect_defense.log"),
    shell:
        """
        python {params.script} \
            --islands  {input.islands} \
            --defense  {input.defense} \
            --gff_dir  {input.gff_dir} \
            --outfile  {output.tsv} \
            --report   {output.report} \
            --min_conf {params.min_conf} \
            > {log} 2>&1
        """


# ── Rule 4: Build island catalog ─────────────────────────────────────────────
rule build_catalog:
    input:
        islands = os.path.join(config["out_dir"], "all_islands_combined.tsv"),
        defense = os.path.join(config["out_dir"], "denovo_defense_intersection.tsv"),
        gff_dir = config["gff_dir"],
    output:
        strains = os.path.join(config["out_dir"], "catalog", "catalog_strains.tsv"),
        groups  = os.path.join(config["out_dir"], "catalog", "catalog_groups.tsv"),
        nesting = os.path.join(config["out_dir"], "catalog", "catalog_nesting.tsv"),
        viewer  = os.path.join(config["out_dir"], "catalog", "island_viewer_data.json"),
    params:
        min_ev   = config.get("catalog_min_evidence", 3),
        out_dir  = os.path.join(config["out_dir"], "catalog"),
        script   = config["catalog_script"],
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30,
    conda:
        "envs/mabs_islands.yaml"
    log:
        os.path.join(config["out_dir"], "logs", "build_catalog.log"),
    shell:
        """
        python {params.script} \
            --islands {input.islands} \
            --defense {input.defense} \
            --gff_dir {input.gff_dir} \
            --out_dir {params.out_dir} \
            --min_ev  {params.min_ev} \
            > {log} 2>&1
        """


# ── Rule 5: Build patient viewer ─────────────────────────────────────────────
rule build_viewer:
    input:
        catalog  = os.path.join(config["out_dir"], "catalog", "island_viewer_data.json"),
        islands  = os.path.join(config["out_dir"], "all_islands_combined.tsv"),
        defense  = os.path.join(config["out_dir"], "denovo_defense_intersection.tsv"),
        rgp_file = config["rgp_file"],
        defense_intersection = config["defense_intersection"],
        genome_stats         = config["genome_stats"],
    output:
        html = os.path.join(config["out_dir"], "patient_comparison_viewer.html"),
    params:
        json_dir     = config.get("json_dir", "genome_tracks"),
        padloc       = config.get("padloc", ""),
        dp           = config.get("defense_predictor", ""),
        prophage_dir = config.get("prophage_dir", ""),
        trna         = config.get("trna_proximity", ""),
        trna_mass    = config.get("trna_proximity_mass", ""),
        trna_boll    = config.get("trna_proximity_boll", ""),
        rgp_mass     = config.get("rgp_file_mass", ""),
        spots_mass   = config.get("spots_file_mass", ""),
        int_mass     = config.get("intersection_mass", ""),
        rgp_boll     = config.get("rgp_file_boll", ""),
        spots_boll   = config.get("spots_file_boll", ""),
        int_boll     = config.get("intersection_boll", ""),
        min_ev       = config.get("min_evidence", 1),
        script       = config["viewer_script"],
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 15,
    conda:
        "envs/mabs_islands.yaml"
    log:
        os.path.join(config["out_dir"], "logs", "build_viewer.log"),
    shell:
        """
        python {params.script} \
            --json_dir        {params.json_dir} \
            --intersection    {input.defense_intersection} \
            --rgp_file        {input.rgp_file} \
            --genome_stats    {input.genome_stats} \
            $([ -n "{params.padloc}"       ] && echo "--padloc {params.padloc}") \
            $([ -n "{params.dp}"           ] && echo "--defense_predictor {params.dp}") \
            $([ -n "{params.prophage_dir}" ] && echo "--prophage_dir {params.prophage_dir}") \
            $([ -n "{params.trna}"         ] && echo "--trna_proximity {params.trna}") \
            $([ -n "{params.trna_mass}"    ] && echo "--trna_proximity_mass {params.trna_mass}") \
            $([ -n "{params.trna_boll}"    ] && echo "--trna_proximity_boll {params.trna_boll}") \
            $([ -n "{params.rgp_mass}"     ] && echo "--rgp_file_mass {params.rgp_mass}") \
            $([ -n "{params.spots_mass}"   ] && echo "--spots_file_mass {params.spots_mass}") \
            $([ -n "{params.int_mass}"     ] && echo "--intersection_mass {params.int_mass}") \
            $([ -n "{params.rgp_boll}"     ] && echo "--rgp_file_boll {params.rgp_boll}") \
            $([ -n "{params.spots_boll}"   ] && echo "--spots_file_boll {params.spots_boll}") \
            $([ -n "{params.int_boll}"     ] && echo "--intersection_boll {params.int_boll}") \
            --scanner_islands {input.islands} \
            --denovo_defense  {input.defense} \
            --island_catalog  {input.catalog} \
            --min_evidence    {params.min_ev} \
            --outfile         {output.html} \
            > {log} 2>&1
        """


# ── Rule 6: BLAST validation (optional) ──────────────────────────────────────
rule extract_blast_seqs:
    input:
        defense  = os.path.join(config["out_dir"], "denovo_defense_intersection.tsv"),
        fasta_dir = config["fasta_dir"],
    output:
        fasta    = os.path.join(config["out_dir"], "blast", "unique_high_conf_fixed.fasta"),
        summary  = os.path.join(config["out_dir"], "blast", "blast_summary.tsv"),
    params:
        out_dir  = os.path.join(config["out_dir"], "blast"),
        conf     = config.get("blast_confidence", "high"),
        flank    = config.get("blast_flank", 500),
        script   = config["extract_script"],
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 15,
    conda:
        "envs/mabs_islands.yaml"
    log:
        os.path.join(config["out_dir"], "logs", "extract_blast.log"),
    shell:
        """
        python {params.script} \
            --denovo_defense {input.defense} \
            --fasta_dir      {input.fasta_dir} \
            --out_dir        {params.out_dir} \
            --conf           {params.conf} \
            --flank          {params.flank} \
            > {log} 2>&1
        """


# ── Utility rules ─────────────────────────────────────────────────────────────
rule clean:
    """Remove all outputs — use with care."""
    shell:
        "rm -rf {config[out_dir]}"


rule report:
    """Print summary statistics for completed run."""
    input:
        combined = os.path.join(config["out_dir"], "all_islands_combined.tsv"),
        catalog  = os.path.join(config["out_dir"], "catalog", "catalog_strains.tsv"),
    run:
        import csv
        from collections import Counter

        islands = []
        with open(input.combined) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                islands.append(row)

        conf = Counter(r["confidence"] for r in islands)
        age  = Counter(r["age_estimate"] for r in islands)
        strains = len(set(r["strain"] for r in islands))

        print("\n" + "="*50)
        print("MAbsIslandScanner — Run Report")
        print("="*50)
        print(f"Strains processed:    {strains}")
        print(f"Total islands:        {len(islands):,}")
        print(f"High confidence:      {conf.get('high',0):,}")
        print(f"Moderate confidence:  {conf.get('moderate',0):,}")
        print(f"Low/very low:         {conf.get('low',0)+conf.get('very_low',0):,}")
        print("\nAge distribution:")
        for a in ["very_recent","recent","moderate","old"]:
            print(f"  {a:15}: {age.get(a,0):,}")
        print("="*50)
