#!/usr/bin/env bash
# =============================================================
# TIR domain search across M. abscessus GD strain FAA files
# Uses Pfam HMM PF01582 (TIR domain)
#
# Usage:
#   bash tir_domain_search.sh \
#       --faa_dir  /path/to/all_faa_files \
#       --out_dir  /path/to/tir_results \
#       --threads  8
# =============================================================
set -euo pipefail

FAA_DIR="" ; OUT_DIR="" ; THREADS=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --faa_dir)  FAA_DIR="$2";  shift 2 ;;
        --out_dir)  OUT_DIR="$2";  shift 2 ;;
        --threads)  THREADS="$2";  shift 2 ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

[[ -z "$FAA_DIR" || -z "$OUT_DIR" ]] && {
    echo "Usage: bash tir_domain_search.sh --faa_dir <dir> --out_dir <dir>"
    exit 1
}

mkdir -p "$OUT_DIR"

# ── Step 1: Download PF01582 HMM from Pfam ───────────────────────────────────
HMM="$OUT_DIR/PF01582_TIR.hmm"
if [[ ! -f "$HMM" ]]; then
    echo "[1/4] Downloading PF01582 (TIR domain) from Pfam..."
    curl -s "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF01582/?annotation=hmm" \
        -o "${HMM}.gz"
    gunzip -f "${HMM}.gz"
    echo "    Downloaded: $HMM"
else
    echo "[1/4] HMM already exists, skipping download."
fi

# Press the HMM database
hmmpress "$HMM" 2>/dev/null || true

# ── Step 2: Build combined FAA with strain tags ───────────────────────────────
echo "[2/4] Building combined FAA..."
COMBINED="$OUT_DIR/all_strains.faa"
> "$COMBINED"
shopt -s nullglob
for faa in "$FAA_DIR"/*.faa; do
    strain=$(basename "$faa" .faa)
    awk -v s="$strain" '/^>/{print $0 " STRAIN=" s} !/^>/{print}' "$faa"
done >> "$COMBINED"
COUNT=$(grep -c "^>" "$COMBINED" || true)
echo "    Combined: $COUNT proteins"

# ── Step 3: Run hmmscan ───────────────────────────────────────────────────────
echo "[3/4] Running hmmscan..."
hmmscan \
    --domtblout "$OUT_DIR/tir_domtbl.txt" \
    --cpu "$THREADS" \
    -E 1e-5 \
    "$HMM" \
    "$COMBINED" \
    > "$OUT_DIR/hmmscan.log" 2>&1
echo "    Done. Results: $OUT_DIR/tir_domtbl.txt"

# ── Step 4: Parse results ─────────────────────────────────────────────────────
echo "[4/4] Parsing results..."
export OUT_DIR
python3 - << 'PYEOF'
import os, re
import pandas as pd

OUT_DIR = os.environ['OUT_DIR']
DOMTBL  = f"{OUT_DIR}/tir_domtbl.txt"

rows = []
with open(DOMTBL) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 23:
            continue
        target   = parts[0]    # protein ID
        query    = parts[3]    # HMM name (PF01582)
        evalue   = float(parts[12])
        score    = float(parts[13])
        dom_start = int(parts[17])
        dom_end   = int(parts[18])

        # Extract strain from STRAIN= tag
        strain_match = re.search(r'STRAIN=(\S+)', target)
        strain = strain_match.group(1) if strain_match else 'unknown'
        protein_id = target.split()[0]

        rows.append({
            'strain':     strain,
            'protein_id': protein_id,
            'evalue':     evalue,
            'score':      score,
            'dom_start':  dom_start,
            'dom_end':    dom_end,
        })

df = pd.DataFrame(rows)
if len(df):
    df = df.sort_values('evalue')
    out_path = f"{OUT_DIR}/tir_hits.tsv"
    df.to_csv(out_path, sep='\t', index=False)
    print(f"\n{'='*50}")
    print(f"TIR domain hits: {len(df)} proteins")
    print(f"Strains with TIR: {df['strain'].nunique()}")
    print(f"\nTop hits:")
    print(df[['strain','protein_id','evalue','score']].head(20).to_string(index=False))
    print(f"\nSaved: {out_path}")
else:
    print("No TIR domain hits found.")
PYEOF

echo ""
echo "Done. Key outputs:"
echo "  $OUT_DIR/tir_hits.tsv      — per-protein TIR hits"
echo "  $OUT_DIR/tir_domtbl.txt    — raw hmmscan domain table"
