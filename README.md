# MAbsIslandScanner

**Genomic island detection and annotation tool for *Mycobacteroides abscessus* clinical genomics**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Overview

MAbsIslandScanner identifies horizontally acquired genomic islands in *M. abscessus* using five independent evidence signals:

| Signal | Description | Threshold |
|---|---|---|
| GC content deviation | Absolute GC% difference from genome mean | > 1% |
| Codon adaptation index (CAI) | CAI ratio relative to ribosomal protein reference set | ≤ 0.92 |
| Mobility genes | Integrase, transposase, recombinase, XerC, RDF | ≥ 1 gene |
| tRNA proximity | Integration within distance of a tRNA gene | ≤ 15 kb |
| Direct repeats | Flanking direct repeat sequences (insertion scar) | 8–25 bp |

Islands with **≥ 3 evidence lines** are retained at moderate or high confidence.

---

## Key features

- **Cohort-scale batch processing** via Snakemake (validated on 250 clinical genomes)
- **Defense system intersection** with DefenseFinder output
- **CAI-based age estimation** — very_recent / recent / moderate / old
- **Cross-species BLAST validation** of mobility gene sequences
- **Interactive patient longitudinal viewer** (self-contained HTML, no server required)
- **100% sensitivity** on validated defense RGPs across 243 clinical isolates
- **864% more islands** detected vs PPanGGOLiN comparative genomics alone

---

## Dependencies

### Core bioinformatics tools

| Tool | Version tested | Purpose | Install |
|---|---|---|---|
| [Prokka](https://github.com/tseemann/prokka) | 1.14.6 | Genome annotation (GFF + FAA) | `conda install -c bioconda prokka` |
| [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN) | 1.2.74 | Pangenome + RGP detection | `conda install -c bioconda ppanggolin` |
| [DefenseFinder](https://github.com/mdmparis/defense-finder) | 1.2.4 | Defense system annotation | `pip install mdmparis-defense-finder` |
| [HMMER](http://hmmer.org/) | 3.4 | HMM profile search | `conda install -c bioconda hmmer` |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | 2.14.0 | Cross-species validation | `conda install -c bioconda blast` |
| [Snakemake](https://snakemake.readthedocs.io/) | 7.32.4 | Workflow management | `conda install -c bioconda snakemake` |

### Optional tools

| Tool | Version tested | Purpose | Install |
|---|---|---|---|
| [PADLOC](https://github.com/padlocbio/padloc) | 2.0.0 | Additional defense system detection | `conda install -c bioconda padloc` |
| [DefensePredictor](https://github.com/gempasteur/DefensePredictor) | 1.0 | ML-based defense gene prediction | See repo |
| [PhiSpy](https://github.com/linsalrob/PhiSpy) | 4.2.21 | Prophage detection | `pip install phispy` |
| [AntiDefenseFinder](https://github.com/mdmparis/defense-finder) | 1.0 | Anti-defense system detection | Integrated with DefenseFinder |

### Python packages

| Package | Version | Purpose |
|---|---|---|
| `biopython` | ≥ 1.79 | Sequence parsing (FASTA, GFF, GenBank) |
| `numpy` | ≥ 1.21 | Numerical computations (GC, CAI) |
| `pandas` | ≥ 1.3 | TSV/CSV handling |
| `scipy` | ≥ 1.7 | Statistical analysis |
| `matplotlib` | ≥ 3.4 | Visualization (optional) |

---

## Installation

### Recommended: conda environment

```bash
git clone https://github.com/cmoyer-x/MAbsIslandScanner.git
cd MAbsIslandScanner

conda env create -f envs/mabs_islands.yaml
conda activate mabs_islands
pip install -e .
```

### Manual install

```bash
conda create -n mabs_islands python=3.10
conda activate mabs_islands

conda install -c bioconda -c conda-forge \
    prokka=1.14.6 \
    ppanggolin \
    hmmer=3.4 \
    blast=2.14.0 \
    snakemake=7.32.4

pip install mdmparis-defense-finder
pip install biopython numpy pandas scipy
pip install -e .
```

---

## Quick start

### Single genome

```bash
# Step 1: Annotate with Prokka
prokka --outdir prokka_out --prefix GD538 \
       --genus Mycobacteroides --species abscessus \
       GD538.fasta

# Step 2: Run scanner
mabs-scan \
    --fasta GD538.fasta \
    --gff   prokka_out/GD538.gff \
    --out   GD538_islands.tsv \
    --min_ev 3
```

### Full cohort (Snakemake)

```bash
# Edit config.yaml to point to your genome directory
snakemake --configfile config.yaml --cores 8
```

---

## Pipeline overview

```
FASTA files
    ├── Prokka ─────────────────► GFF + FAA
    ├── PPanGGOLiN ─────────────► RGPs (variable islands)
    ├── DefenseFinder ──────────► Defense system coordinates
    └── MAbsIslandScanner
            ├── GC deviation
            ├── CAI (HMMER vs ribosomal HMMs)
            ├── Mobility gene detection
            ├── tRNA proximity
            └── Direct repeat search
                    ▼
            Island TSVs per strain
                    ▼
            build_island_catalog.py ──► Stable IDs (Mabs_GI_001)
                    ▼
            intersect_denovo_defense.py ──► Defense island catalog
                    ▼
            build_patient_viewer.py ──► Interactive HTML viewers
```

---

## Output files

| File | Description |
|---|---|
| `*_islands.tsv` | Per-strain predictions with 5 evidence signals |
| `all_islands_combined.tsv` | Full cohort combined |
| `catalog/catalog_strains.tsv` | Stable strain-level IDs (GDxx_GI_001) |
| `catalog/catalog_groups.tsv` | Cross-strain group IDs (Mabs_GI_001) |
| `denovo_defense_intersection.tsv` | Fixed defense islands |
| `fixed_defense_catalog_final.tsv` | Full catalog with BLAST validation + age estimates |
| `patient_viewers/` | Self-contained HTML viewers per patient |

### Island TSV columns

| Column | Description |
|---|---|
| `island_id` | Stable ID (e.g. GD538_GI_010) |
| `start`, `end`, `length` | Chromosomal coordinates (bp) |
| `n_evidence` | Number of evidence lines (2–5) |
| `confidence` | low / moderate / high |
| `age_estimate` | very_recent / recent / moderate / old |
| `cai_ratio` | CAI relative to host ribosomal proteins |
| `dominant_cargo` | defense / mobility / metal / efflux / regulatory / hypothetical |
| `trna_flanked` | Boolean — integration at tRNA locus |
| `has_dr` | Boolean — direct repeat detected |
| `mob_types` | Mobility gene types detected |
| `blast_validated` | Boolean — cross-species BLAST confirmation |
| `blast_hit_species` | Species with homologous mobility genes |

---

## Performance (250 M. abscessus clinical isolates)

| Metric | Value |
|---|---|
| Strains processed | 312 / 316 across 3 subspecies (98.7%) |
| — M. abscessus | 250 / 253 strains |
| — M. massiliense | 54 / 55 strains |
| — M. bolletii | 8 / 8 strains |
| Total islands detected | 10,314 |
| Mean islands per genome | 34.5 abscessus · 27.8 massiliense · 24.2 bolletii |
| Defense islands (abscessus) | 829 (86 variable + 743 fixed) |
| Increase over PPanGGOLiN alone | 864% |
| BLAST-validated islands | 12 |
| Sensitivity on known defense RGPs | 100% (86/86) |
| Direct repeat detection rate | 75.9% of known defense RGPs |
| tRNA proximity enrichment | 19.8% vs ~2–3% expected |

---

## Tool comparison — GD538

| Feature | IslandViewer 4 | MAbsIslandScanner |
|---|---|---|
| CBASS_II island | 2 fragments (88% + 4 kb orphan) | 1 coherent 65 kb island |
| Hna island | 5 fragments (53% max overlap) | 1 coherent 103 kb island |
| Defense annotation | None | Teal = defense cargo |
| Age estimate | None | very_recent / recent / old |
| Mechanism | None | integrase / transposase / XerC |
| Evidence score | None | 3–5 independent lines |
| Organism calibration | Gram-negative (40–55% GC) | *M. abscessus* (64% GC) |
| Cohort scale | Single genome | 250 genomes · Snakemake |

IslandViewer 4: Bertelli et al. *Nucleic Acids Res* 2017; 45:W30–W35.

---

## Citation

Moyer et al. (in preparation). MAbsIslandScanner: a genomic island detection and annotation tool for *Mycobacteroides abscessus* clinical cohorts. University of Pittsburgh.

Please also cite:
- **PPanGGOLiN**: Gautreau et al. *PLOS Comput Biol* 2020
- **DefenseFinder**: Tesson et al. *Nat Commun* 2022
- **Prokka**: Seemann. *Bioinformatics* 2014; 30(14):2068–2069
- **HMMER**: Eddy. *PLOS Comput Biol* 2011

---

## License

MIT — see [LICENSE](LICENSE)

---

## Contact

Casey Moyer · Research Scientist · University of Pittsburgh  
GitHub: [@cmoyer-x](https://github.com/cmoyer-x)
