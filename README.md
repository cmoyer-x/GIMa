# GIMa

A bioinformatics pipeline for systematic genomic island detection, defense system characterization, and phage infection analysis across *Mycobacterium abscessus* clinical cohorts.

Developed at the University of Pittsburgh (Hatfull Lab) for analysis of *M. abscessus* clinical isolates spanning three subspecies (*abscessus*, *massiliense*, *bolletii*).

---

## Overview

GIMa integrates multiple genomic island prediction and defense system annotation tools into a unified Snakemake pipeline, producing per-patient interactive HTML viewers that link genomic island content to phage infection (EOP) data. It was developed to address questions about how horizontally acquired defense islands shape phage susceptibility in clinical *M. abscessus* infections.

**Key outputs:**
- Catalog of genomic islands across your strains with age estimates, HGT source, and transfer mechanism
- Per-strain defense system calls (DefenseFinder + PADLOC) linked to island coordinates
- Per-strain interactive patient HTML viewers with genome track, defense systems, EOP heatmap, and confidence filtering
- TPP locus integrity assessment for each strain
- TIR domain protein detection via Pfam PF01582 hmmscan

---

## Repository structure

```
GIMa/
├── Snakefile                        # Main pipeline
├── config/
│   └── config.yaml                  # Pipeline configuration
├── workflow/
│   ├── rules/                       # Snakemake rule modules
│   └── scripts/
│       ├── build_island_catalog.py  # Build genomic island catalog across strains
│       ├── run_defensefinder.py     # DefenseFinder wrapper
│       ├── run_padloc.py            # PADLOC wrapper
│       ├── assign_island_age.py     # CAI-based age estimation
│       ├── blast_validate.py        # Cross-strain BLAST validation
│       └── generate_patient_viewers.py  # HTML viewer generation
├── postprocessing/
│   ├── patch_viewers.py             # Collapse fragmentation duplicates in viewers
│   ├── update_viewer_ui.py          # Add confidence toggle, remove legend
│   └── add_eop_to_viewers.py        # Inject EOP heatmap into viewers
├── analysis/
│   ├── tpp_pipeline.sh              # TPP locus BLAST pipeline
│   ├── fetch_tpp_refs.sh            # Fetch TPP reference proteins
│   ├── run_tpp_blast.sh             # BLAST all strains vs TPP locus
│   ├── tpp_merge_eop.py             # Merge TPP status with EOP data
│   └── tir_domain_search.sh         # Pfam PF01582 TIR domain hmmscan
├── results/
│   ├── fixed_defense_catalog_final.tsv      # Raw defense catalog
│   ├── defense_catalog_collapsed.tsv        # Deduplicated catalog
│   └── padloc_all_systems.csv               # PADLOC results all strains
└── README.md
```

---

## Installation

### Dependencies

```bash
# Core pipeline
conda create -n GIMa python=3.10
conda activate GIMa
conda install -c bioconda snakemake prokka hmmer blast
conda install -c conda-forge pandas numpy

# Defense system tools
pip install defensefinder
conda install -c bioconda padloc

# PPanGGOLiN (genomic island detection)
conda install -c bioconda ppanggolin
```

### Clone and configure

```bash
git clone https://github.com/cmoyer-x/GIMa.git
cd GIMa
```

Edit `config/config.yaml` to set:
- `faa_dir`: path to Prokka FAA files for your strains
- `gff_dir`: path to Prokka GFF files
- `genome_dir`: path to assembled FASTA files
- `output_dir`: where results will be written

---

## Running the pipeline

```bash
# Dry run first
snakemake --dry-run --cores 16

# Full run
snakemake --cores 16 --use-conda
```

### Individual steps

**Generate patient HTML viewers:**
```bash
python workflow/scripts/generate_patient_viewers.py \
    --defense_catalog results/fixed_defense_catalog_final.tsv \
    --padloc          results/padloc_all_systems.csv \
    --eop_csv         data/EOP.csv \
    --out_dir         results/patient_viewers
```

---

## Postprocessing viewer scripts

These scripts update the HTML viewers after generation. Run them in order:

### Step 1 — Collapse fragmentation duplicates

Fragmented (draft) assemblies cause single genomic islands to appear as multiple entries with the same catalog group ID. This script collapses them, rescuing age/CAI signals from partial fragments:

```bash
python postprocessing/patch_viewers.py \
    --viewer_dir results/patient_viewers \
    --out_dir    results/patient_viewers_patched
```

**Logic:**
- Groups catalog islands by `group_id` (e.g. `Mabs_GI_001`)
- Prefers `unique` status entries over `trimmed` (contig-edge fragments)
- Takes most recent age estimate from any fragment
- Takes highest confidence from any fragment
- Merges evidence strings
- Adds `has_recent_cai` flag: `True` if any fragment had CAI deviation + recent age
- Adds `n_contig_fragments` and `n_trimmed_fragments` counts

### Step 2 — Update UI (remove legend, add confidence toggle)

```bash
python postprocessing/update_viewer_ui.py \
    --viewer_dir results/patient_viewers_patched \
    --out_dir    results/patient_viewers_ui
```

**Changes:**
- Removes the color legend from the bottom of each viewer
- Adds a **GI confidence** toggle to the controls bar
- Default view shows **high confidence** genomic islands only
- Toggle switches to **high + moderate** for broader investigation

### Step 3 — Add EOP heatmap

```bash
python postprocessing/add_eop_to_viewers.py \
    --viewer_dir results/patient_viewers_ui \
    --eop_csv    data/EOP.csv \
    --out_dir    results/patient_viewers_final
```

**Adds a phage infection panel** below the genome track showing log₁₀ EOP values for each tested phage, grouped by family (Muddy, BPs, ZoeJ, other) and color-coded by infection efficiency.

| Color | EOP range | Interpretation |
|-------|-----------|----------------|
| Bright green | 0 to −1 | Productive infection |
| Light green | −1 to −3 | Low productive |
| Amber | −3 to −5 | Intermediate |
| Red | −5 to −7 | Resistant |
| Dark red | below −7 | Highly resistant |
| Grey | — | Not tested |

---

## TPP locus analysis

The trehalose polyphleate (TPP) biosynthesis locus is the primary receptor for phages Muddy and BPs (Wetzel et al. 2023, *Nature Microbiology*). This pipeline assesses TPP integrity across all sequenced strains.

**Five TPP locus genes (ATCC 19977 reference):**

| Gene | Function |
|------|----------|
| MAB_0939 (Pks) | Polyketide synthase — builds phleic acid chains |
| MAB_0940 (PE) | Transacylase — transfers phleic acids onto DAT |
| MAB_0941 (PapA3) | Acyltransferase — forms diacyltrehalose precursor |
| MAB_0942 (MmpL10) | Transporter — exports DAT across membrane |
| MAB_0943 (FadD23) | Fatty acyl-AMP ligase — activates fatty acid substrate |

```bash
# Fetch reference proteins (uses UniProt REST API)
python3 << 'PYEOF'
import urllib.request, time, re
OUT = "analysis/tpp_reference.faa"
entries = [("B1MB11","MAB_0939_Pks"),("B1MB12","MAB_0940_PE"),
           ("B1MB13","MAB_0941_PapA3"),("B1MB14","MAB_0942_MmpL10"),
           ("B1MB15","MAB_0943_FadD23")]
with open(OUT,'w') as out:
    for acc, label in entries:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        with urllib.request.urlopen(url, timeout=30) as r:
            fasta = r.read().decode()
        lines = fasta.strip().split('\n')
        out.write(f">{label}\n{''.join(lines[1:])}\n")
        time.sleep(0.3)
PYEOF

# Run BLAST pipeline
bash analysis/run_tpp_blast.sh \
    --faa_dir /path/to/prokka_faa_files \
    --out_dir analysis/tpp_results \
    --threads 8

# Merge with EOP data
python analysis/tpp_merge_eop.py \
    --tpp  analysis/tpp_results/tpp_status_per_strain.tsv \
    --eop  data/EOP.csv \
    --out  analysis/tpp_results/merged_tpp_eop.tsv
```

**TPP status classifications:**
- `intact` — all 5 genes present at ≥80% identity and ≥80% query coverage
- `partial` — 3–4 genes intact
- `disrupted` — 1–2 genes found but degraded
- `absent` — no TPP genes detected

---

## TIR domain detection

TIR (Toll/interleukin-1 receptor) domain proteins are standalone anti-phage defense systems not covered by DefenseFinder or PADLOC. This pipeline detects them via Pfam `PF01582`:

```bash
bash analysis/tir_domain_search.sh \
    --faa_dir  /path/to/all_faa_files \
    --out_dir  analysis/tir_results \
    --threads  8
```

Requires `hmmer` (hmmscan). Downloads `PF01582.hmm` from the EBI Pfam API automatically.

**Output:** `tir_results/tir_hits.tsv` — strain, protein ID, e-value, score, domain coordinates.

---

## Biological findings

GIMa was developed and applied to a clinical *M. abscessus* cohort. Detailed biological findings from that work — covering defense island prevalence, subspecies-level phage susceptibility, longitudinal horizontal gene transfer, and anti-restriction protein characterization — will be reported in a forthcoming publication (Moyer et al., in preparation).

---

## Dependencies and versions

| Tool | Version | Purpose |
|------|---------|---------|
| PPanGGOLiN | ≥2.0 | Genomic island detection |
| DefenseFinder | ≥1.3 | Defense system annotation |
| PADLOC | ≥2.0 | Defense system annotation |
| BLAST+ | ≥2.12 | Sequence similarity |
| HMMER | ≥3.3 | TIR/TPP domain search |
| Prokka | ≥1.14 | Genome annotation |
| Python | ≥3.10 | Pipeline scripts |
| pandas | ≥2.0 | Data processing |

---

## Citation

If you use GIMa, please cite:

> Moyer CL, et al. *Genomic island-mediated defense system acquisition shapes phage susceptibility in clinical Mycobacterium abscessus.* (in preparation)

Please also cite the underlying tools: PPanGGOLiN, DefenseFinder, PADLOC, and HMMER.

---

## Contact

Casey Moyer — Research Scientist, University of Pittsburgh  
GitHub: [@cmoyer-x](https://github.com/cmoyer-x)

---

## License

MIT License. See `LICENSE` for details.
