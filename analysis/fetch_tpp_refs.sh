#!/usr/bin/env bash
# Fetch TPP reference proteins from UniProt (no NCBI/efetch needed)
# Usage: bash fetch_tpp_refs.sh /path/to/tpp_reference.faa
OUT="${1:?Usage: bash fetch_tpp_refs.sh /path/to/tpp_reference.faa}"

python3 << PYEOF
import urllib.request, time
OUT = "$OUT"
entries = [
    ("B1MB11","MAB_0939_Pks"),("B1MB12","MAB_0940_PE"),
    ("B1MB13","MAB_0941_PapA3"),("B1MB14","MAB_0942_MmpL10"),
    ("B1MB15","MAB_0943_FadD23")
]
with open(OUT,'w') as f:
    for acc, label in entries:
        print(f"Fetching {acc} ({label})...", end=' ', flush=True)
        try:
            with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30) as r:
                fasta = r.read().decode()
            lines = fasta.strip().split('\n')
            f.write(f">{label}\n{''.join(lines[1:])}\n")
            print(f"OK ({len(''.join(lines[1:]))} aa)")
            time.sleep(0.3)
        except Exception as e:
            print(f"ERROR: {e}")
print("\nHeaders:")
with open(OUT) as f:
    [print(' ',l.strip()) for l in f if l.startswith('>')]
PYEOF
