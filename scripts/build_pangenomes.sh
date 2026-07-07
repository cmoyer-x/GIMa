#!/usr/bin/env bash
set -euo pipefail

: <<'DOC'
build_pangenomes.sh

Generates the per-subspecies PPanGGOLiN inputs that GIMa (Component 2)
consumes: regions_of_genomic_plasticity.tsv and spots.tsv per subspecies.

Component 1 for GIMa. Groups strains by subspecies from the cohort master,
stages Prokka GFFs, and runs the PPanGGOLiN chain
(workflow -> rgp -> spot -> write_pangenome) once per subspecies.

Prerequisites:
  - PPanGGOLiN 2.3.0 (conda env "GIMa")
  - Per-strain Prokka GFFs at cohort/genomes/<strain>/prokka/<strain>.gff
  - Cohort master CSV with a "subspecies" column

Usage:
  conda activate GIMa
  ./build_pangenomes.sh                         run all three subspecies
  ./build_pangenomes.sh bolletii                run one subspecies
  MASTER=data/mabs_cohort_master.csv ./build_pangenomes.sh

Environment:
  MASTER     path to cohort master CSV   (default: data/mabs_cohort_master.csv)
  GENOMES    per-strain genomes root     (default: cohort/genomes)
  CPU        threads for workflow step   (default: 8)
DOC

MASTER="${MASTER:-data/mabs_cohort_master.csv}"
GENOMES="${GENOMES:-cohort/genomes}"
CPU="${CPU:-8}"

write_strain_lists () {
    python3 - "$MASTER" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], encoding='utf-8-sig')))
def is_mab(r):
    return (r.get('species') or '').strip() in (
        'M. abscessus', 'Mycobacteroides abscessus', 'abscessus')
def sequenced(r):
    return (r.get('sequenced') or '').upper() in ('TRUE', 'YES', '1')
for sub in ('abscessus', 'massiliense', 'bolletii'):
    strains = [r['strain'] for r in rows
               if is_mab(r) and sequenced(r)
               and (r.get('subspecies') or '').strip() == sub]
    open(f'{sub}_strains.txt', 'w').write('\n'.join(strains) + '\n')
    print(f'{sub}: {len(strains)} strains', file=sys.stderr)
PY
}

stage_gffs () {
    local sub=$1
    mkdir -p "pangenome_input/$sub"
    local missing=0
    while read -r s; do
        [ -z "$s" ] && continue
        local gff="$GENOMES/$s/prokka/$s.gff"
        if [ -f "$gff" ]; then
            cp "$gff" "pangenome_input/$sub/$s.gff"
        else
            echo "MISSING: $gff" >&2
            missing=$((missing + 1))
        fi
    done < "${sub}_strains.txt"
    local n
    n=$(ls "pangenome_input/$sub" | wc -l)
    echo "$sub staged: $n GFFs ($missing missing)" >&2
}

build_gff_list () {
    local sub=$1
    : > "${sub}_gff_list.tsv"
    for f in "pangenome_input/$sub"/*.gff; do
        local s
        s=$(basename "$f" .gff)
        printf '%s\t%s\n' "$s" "$(readlink -f "$f")" >> "${sub}_gff_list.tsv"
    done
    echo "$sub: $(wc -l < "${sub}_gff_list.tsv") genomes in list" >&2
}

run_ppanggolin () {
    local sub=$1
    ppanggolin workflow --anno "${sub}_gff_list.tsv" \
        --output "${sub}_pangenome" --cpu "$CPU"
    ppanggolin rgp  -p "${sub}_pangenome/pangenome.h5"
    ppanggolin spot -p "${sub}_pangenome/pangenome.h5"
    ppanggolin write_pangenome -p "${sub}_pangenome/pangenome.h5" \
        --output "${sub}_pangenome/rgp_output" --regions --spots
}

verify () {
    local sub=$1
    local rgp="${sub}_pangenome/rgp_output/regions_of_genomic_plasticity.tsv"
    local spots="${sub}_pangenome/rgp_output/spots.tsv"
    if [ -f "$rgp" ]; then
        echo "$sub RGP:   $(($(wc -l < "$rgp") - 1)) regions" >&2
    else
        echo "$sub RGP:   MISSING" >&2
    fi
    if [ -f "$spots" ]; then
        echo "$sub spots: $(($(wc -l < "$spots") - 1))" >&2
    else
        echo "$sub spots: MISSING" >&2
    fi
}

process_subspecies () {
    local sub=$1
    echo ">>> $sub" >&2
    stage_gffs "$sub"
    build_gff_list "$sub"
    run_ppanggolin "$sub"
    verify "$sub"
}

main () {
    write_strain_lists
    if [ "$#" -ge 1 ]; then
        for sub in "$@"; do
            process_subspecies "$sub"
        done
    else
        process_subspecies bolletii
        process_subspecies massiliense
        process_subspecies abscessus
    fi
    echo "done" >&2
}

main "$@"
