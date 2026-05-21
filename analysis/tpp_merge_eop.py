#!/usr/bin/env python3
"""Merge TPP locus status with phage EOP data.

Usage:
    python3 tpp_merge_eop.py \
        --tpp  /path/to/tpp_status_per_strain.tsv \
        --eop  /path/to/EOP.csv \
        --out  /path/to/merged_tpp_eop.tsv
"""
import argparse, pandas as pd, numpy as np

def to_log10_eop(v):
    try:
        f = float(str(v).strip())
        return round(np.log10(f), 1) if f > 0 else None
    except: return None

parser = argparse.ArgumentParser()
parser.add_argument('--tpp', required=True)
parser.add_argument('--eop', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()

tpp = pd.read_csv(args.tpp, sep='\t')
print(f"TPP strains: {len(tpp)}")
print(tpp['tpp_status'].value_counts())

df = pd.read_csv(args.eop, header=None)
strains = list(df.iloc[1].astype(str).values[1:])
phage_rows = {'Muddy': 14, 'Muddy_HRM': 15, 'BPs_HRM10': 8, 'BPs_WT': 12}

eop_data = {}
for phage, ridx in phage_rows.items():
    for s, v in zip(strains, df.iloc[ridx].values[1:]):
        eop = to_log10_eop(v)
        if eop is not None:
            eop_data.setdefault(str(s).strip(), {})[phage] = round(eop, 1)

eop_df = pd.DataFrame.from_dict(eop_data, orient='index').reset_index()
eop_df.columns = ['strain'] + list(phage_rows.keys())

merged = tpp.merge(eop_df, on='strain', how='inner')
merged.to_csv(args.out, sep='\t', index=False)
print(f"\nMerged strains: {len(merged)}")
print("\n--- Muddy EOP by TPP status ---")
for status in ['intact','partial','disrupted','absent']:
    sub = merged[merged['tpp_status']==status]['Muddy'].dropna()
    if len(sub):
        prod = sum(sub >= -3)
        print(f"  {status:12s}: n={len(sub):3d}, mean={sub.mean():.2f}, productive={prod} ({100*prod//len(sub)}%)")
print(f"Saved: {args.out}")
