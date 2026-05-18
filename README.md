# MAbsIslandScanner

Genomic island detection and annotation tool calibrated for *Mycobacteroides abscessus* clinical genomics.

## Overview

MAbsIslandScanner identifies horizontally acquired genomic islands in *M. abscessus* using five independent evidence signals:

1. **GC content deviation** — regions with GC% differing from genome mean
2. **Codon adaptation index (CAI)** — codon usage foreign to the host
3. **Mobility genes** — integrases, transposases, recombinases, XerC
4. **tRNA proximity** — integration within 15kb of a tRNA gene
5. **Direct repeats** — flanking sequence scars of insertion (8–25bp)

Islands with ≥3 evidence lines are retained. The tool integrates PPanGGOLiN RGP-guided detection for variable islands and de novo sliding-window detection for fixed islands invisible to comparative genomics.

## Key features

- Cohort-scale batch processing via Snakemake (250+ genomes)
- Defense system intersection with DefenseFinder output
- CAI-based age estimation (very_recent / recent / moderate / old)
- Cross-species BLAST validation of mobility genes
- Interactive patient longitudinal viewer (HTML, no server required)
- 100% sensitivity on validated defense RGPs across 243 clinical isolates

## Installation

```bash
conda env create -f envs/mabs_islands.yaml
conda activate mabs_islands
pip install -e .
```

## Quick start

```bash
# Single genome
mabs-scan \
    --fasta SAMPLE.fasta \
    --gff   SAMPLE_prokka.gff \
    --out   SAMPLE_islands.tsv

# Full cohort via Snakemake
snakemake --configfile config.yaml --cores 8
```

## Output

| File | Description |
|---|---|
| `*_islands.tsv` | Per-strain island predictions with evidence lines |
| `all_islands_combined.tsv` | Full cohort combined |
| `catalog_strains.tsv` | Stable island IDs (SAMPLE_GI_001) |
| `catalog_groups.tsv` | Cross-strain group IDs (Mabs_GI_001) |
| `patient_viewers/` | Self-contained HTML viewers per patient |

## Citation

Moyer et al. (in preparation). MAbsIslandScanner: a genomic island detection tool for *Mycobacteroides abscessus* clinical cohorts.

## License

MIT
