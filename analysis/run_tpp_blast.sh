#!/usr/bin/env bash
# TPP BLAST pipeline — run after fetching reference proteins
# Usage: bash run_tpp_blast.sh --faa_dir <dir> --out_dir <dir> [--threads N]
set -euo pipefail

FAA_DIR="" ; OUT_DIR="" ; THREADS=8
while [[ $# -gt 0 ]]; do
    case "$1" in
        --faa_dir)  FAA_DIR="$2"; shift 2 ;;
        --out_dir)  OUT_DIR="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done
[[ -z "$FAA_DIR" || -z "$OUT_DIR" ]] && { echo "Usage: bash run_tpp_blast.sh --faa_dir <dir> --out_dir <dir>"; exit 1; }

REF="$OUT_DIR/tpp_reference.faa"
[[ ! -f "$REF" ]] && { echo "ERROR: $REF not found. Run fetch_tpp_refs.sh first."; exit 1; }
mkdir -p "$OUT_DIR/blast_results"

echo "[1/3] Building combined FAA..."
COMBINED="$OUT_DIR/all_strains_combined.faa"
> "$COMBINED"
shopt -s nullglob
FAA_FILES=("$FAA_DIR"/*.faa)
[[ ${#FAA_FILES[@]} -eq 0 ]] && { echo "ERROR: No .faa files in $FAA_DIR"; exit 1; }
for faa in "${FAA_FILES[@]}"; do
    strain=$(basename "$faa" .faa)
    awk -v s="$strain" '/^>/{print $0 " STRAIN=" s} !/^>/{print}' "$faa"
done >> "$COMBINED"
echo "    $(grep -c '^>' "$COMBINED") proteins from ${#FAA_FILES[@]} strains"

echo "[2/3] Running BLAST..."
makeblastdb -in "$REF" -dbtype prot -out "$OUT_DIR/tpp_ref_db" -quiet
blastp \
    -query "$COMBINED" -db "$OUT_DIR/tpp_ref_db" \
    -out "$OUT_DIR/blast_results/all_vs_tpp.tsv" \
    -outfmt "6 qseqid sseqid pident length qlen slen qcovs bitscore evalue stitle" \
    -evalue 1e-10 -num_threads "$THREADS" -max_target_seqs 1
echo "    $(wc -l < "$OUT_DIR/blast_results/all_vs_tpp.tsv") hits"

echo "[3/3] Parsing results..."
export OUT_DIR
python3 - << 'PYEOF'
import pandas as pd, numpy as np, os, json

OUT_DIR = os.environ['OUT_DIR']
COLS = ['qseqid','sseqid','pident','length','qlen','slen','qcovs','bitscore','evalue','stitle']
df = pd.read_csv(f"{OUT_DIR}/blast_results/all_vs_tpp.tsv", sep='\t', names=COLS)
df['strain']   = df['qseqid'].str.extract(r'STRAIN=(\S+)')
df['tpp_gene'] = df['sseqid'].str.extract(r'(MAB_\d+c?_\w+)')
df['pident']   = df['pident'].astype(float)
df['qcovs']    = df['qcovs'].astype(float)
df['intact']   = (df['pident'] >= 80) & (df['qcovs'] >= 80)

best = df.sort_values('bitscore', ascending=False).groupby(['strain','tpp_gene'], as_index=False).first()
TPP_GENES = ['MAB_0939_Pks','MAB_0940_PE','MAB_0941_PapA3','MAB_0942_MmpL10','MAB_0943_FadD23']

rows = []
for strain in df['strain'].dropna().unique():
    sub = best[best['strain'] == strain]
    row = {'strain': strain}
    n_found = n_intact = 0
    for gene in TPP_GENES:
        hit = sub[sub['tpp_gene'] == gene]
        if len(hit):
            h = hit.iloc[0]
            row[f'{gene}_pident'] = round(h['pident'], 1)
            row[f'{gene}_qcovs']  = round(h['qcovs'], 1)
            row[f'{gene}_intact'] = bool(h['intact'])
            n_found += 1; n_intact += int(h['intact'])
        else:
            row[f'{gene}_pident'] = None; row[f'{gene}_qcovs'] = None; row[f'{gene}_intact'] = False
    row['genes_found'] = n_found; row['genes_intact'] = n_intact
    row['tpp_status'] = 'intact' if n_intact==5 else 'partial' if n_intact>=3 else 'disrupted' if n_found>=1 else 'absent'
    rows.append(row)

out = pd.DataFrame(rows)
path = f"{OUT_DIR}/tpp_status_per_strain.tsv"
out.to_csv(path, sep='\t', index=False)
print(f"\nTPP status across {len(out)} strains:")
print(out['tpp_status'].value_counts().to_string())
print(f"Saved: {path}")
PYEOF
echo "Done. Output: $OUT_DIR/tpp_status_per_strain.tsv"
