# Generating pangenome inputs (PPanGGOLiN RGP) for GIMa

GIMa is **Component 2**: it scans genomic islands, builds the catalog, and
generates viewers. It assumes the per-subspecies PPanGGOLiN outputs already
exist. This document describes **Component 1** — how to generate those inputs
from per-strain Prokka annotations.

GIMa's `Snakefile` header lists the precomputed inputs it requires. The ones
produced here are:

- `<subspecies>_pangenome/rgp_output/regions_of_genomic_plasticity.tsv`
- `<subspecies>_pangenome/rgp_output/spots.tsv`
- `<subspecies>_pangenome/genomes_statistics.tsv`

GIMa runs one pangenome **per subspecies** (`abscessus`, `massiliense`,
`bolletii`) because regions of genomic plasticity (RGPs) are only meaningful
within a set of related genomes; mixing subspecies would inflate the accessory
genome and blur RGP boundaries.

---

## Requirements

- PPanGGOLiN 2.3.0 (in the `GIMa` conda environment)
- Per-strain Prokka output for every genome. GIMa uses the Prokka `.gff`
  (which includes the genome sequence in a `##FASTA` block by default).
- One strain-list file per subspecies (see Step 1). PPanGGOLiN is run
  separately per subspecies, so RGP boundaries are not blurred by mixing
  divergent genomes.

This document assumes you have already run Prokka on each assembly, e.g.:

```bash
prokka --outdir prokka/<strain> --prefix <strain> --cpus 8 <strain>.fasta
```

giving a `prokka/<strain>/<strain>.gff` per genome.

---

## Step 1 — assign each genome to a subspecies

PPanGGOLiN must be run once per subspecies (*abscessus*, *massiliense*,
*bolletii*). Assign each genome by MASH distance to the three subspecies
reference genomes, then write one strain-list file per subspecies
(`abscessus_strains.txt`, `massiliense_strains.txt`, `bolletii_strains.txt`),
one strain name per line.

```bash
mash sketch -o refs \
    abscessus_ref.fasta massiliense_ref.fasta bolletii_ref.fasta
for asm in assemblies/*.fasta; do
    s=$(basename "$asm" .fasta)
    nearest=$(mash dist refs.msh "$asm" \
              | sort -k3,3g | head -1 | cut -f1)
    case "$nearest" in
        *abscessus*)   echo "$s" >> abscessus_strains.txt ;;
        *massiliense*) echo "$s" >> massiliense_strains.txt ;;
        *bolletii*)    echo "$s" >> bolletii_strains.txt ;;
    esac
done
wc -l *_strains.txt
```

Use a correct *bolletii* reference; an incorrect one causes bolletii genomes to
be misassigned to *abscessus* or *massiliense*. If you already have subspecies
assignments from another workflow, just produce the three `*_strains.txt` files
by whatever means and skip this step.

## Step 2 — stage GFFs into a flat per-subspecies directory

Point `PROKKA_DIR` at wherever your per-strain Prokka GFFs live. This example
assumes `prokka/<strain>/<strain>.gff`.

```bash
PROKKA_DIR=prokka
for sub in abscessus massiliense bolletii; do
    mkdir -p pangenome_input/$sub
    while read s; do
        gff="$PROKKA_DIR/$s/$s.gff"
        [ -f "$gff" ] && cp "$gff" "pangenome_input/$sub/$s.gff" \
                      || echo "MISSING: $gff"
    done < ${sub}_strains.txt
    echo "$sub staged: $(ls pangenome_input/$sub/ | wc -l) GFFs"
done
```

Any `MISSING` line indicates a strain whose Prokka GFF was not found at the
expected path — adjust `PROKKA_DIR` or the filename pattern to match your
layout.

## Step 3 — build the PPanGGOLiN annotation list per subspecies

PPanGGOLiN takes a two-column TSV: `<name>\t<absolute path to GFF>`.

```bash
for sub in abscessus massiliense bolletii; do
    for f in pangenome_input/$sub/*.gff; do
        s=$(basename "$f" .gff)
        echo -e "$s\t$(readlink -f "$f")"
    done > ${sub}_gff_list.tsv
    echo "$sub: $(wc -l < ${sub}_gff_list.tsv) genomes"
done
```

## Step 4 — run PPanGGOLiN per subspecies

Four commands per subspecies. `workflow` stops after partitioning; RGP and spot
detection are separate subcommands and must be run explicitly. `rgp` and `spot`
do **not** accept `--cpu`.

```bash
run_pangenome () {
    sub=$1
    ppanggolin workflow --anno ${sub}_gff_list.tsv \
        --output ${sub}_pangenome --cpu 8
    ppanggolin rgp  -p ${sub}_pangenome/pangenome.h5
    ppanggolin spot -p ${sub}_pangenome/pangenome.h5
    ppanggolin write_pangenome -p ${sub}_pangenome/pangenome.h5 \
        --output ${sub}_pangenome/rgp_output --regions --spots
}

run_pangenome bolletii       # ~10 genomes, seconds
run_pangenome massiliense    # ~76 genomes, minutes
```

For the large abscessus set (290 genomes), run in a detachable session so a
dropped connection does not kill it:

```bash
nohup bash -c '
  ppanggolin workflow --anno abscessus_gff_list.tsv \
      --output abscessus_pangenome --cpu 8 && \
  ppanggolin rgp  -p abscessus_pangenome/pangenome.h5 && \
  ppanggolin spot -p abscessus_pangenome/pangenome.h5 && \
  ppanggolin write_pangenome -p abscessus_pangenome/pangenome.h5 \
      --output abscessus_pangenome/rgp_output --regions --spots
' > abscessus_pang.log 2>&1 &

tail -f abscessus_pang.log
```

## Step 5 — verify outputs

```bash
for sub in abscessus massiliense bolletii; do
    echo "=== $sub ==="
    wc -l ${sub}_pangenome/rgp_output/regions_of_genomic_plasticity.tsv
    wc -l ${sub}_pangenome/rgp_output/spots.tsv
done
```

Each `regions_of_genomic_plasticity.tsv` has columns:
`region, genome, contig, genes, first_gene, last_gene, start, stop, length,
coordinates, score, contigBorder, wholeContig`.

These paths match GIMa's `config/config.yaml`
(`rgp_file`, `rgp_file_mass`, `rgp_file_boll`, and the corresponding
`spots_file_*`).

---

## Notes and caveats

- **Small sets partition poorly.** For a subspecies with few genomes (e.g.
  bolletii, n=10), PPanGGOLiN warns that the genome count is too low to robustly
  partition the graph. RGP detection still runs, but RGP calls for very small
  sets are lower confidence than for larger sets (abscessus, massiliense).
- **Assembly fragmentation inflates RGP/island counts.** Fragmented assemblies
  produce contig-border artifacts that can appear as spurious plasticity. The
  RGP output flags `contigBorder`; downstream GIMa postprocessing collapses
  fragmentation duplicates.
- **`workflow` does not run rgp/spot.** This is the most common mistake — the
  `workflow` subcommand ends at partitioning. The `rgp`, `spot`, and
  `write_pangenome --regions --spots` calls are required to produce the files
  GIMa reads.
