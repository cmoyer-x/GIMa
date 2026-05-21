#!/usr/bin/env python3
"""
Patch patient HTML viewers to collapse fragmentation duplicates in catalog_islands,
deduplicate defense system entries (df/padloc), and add has_recent_cai flags.

Handles:
  - Multiple sibling islands with same group_id from contig breaks
  - 'trimmed' status islands (partial fragments at contig edges)
  - Prefers 'unique' status over 'trimmed' as the base entry
  - Rescues 'recent' age signal from any fragment

Usage:
    # Patch ALL viewers (recommended):
    python3 patch_viewers.py \
        --viewer_dir /path/to/patient_viewers \
        --out_dir    /path/to/patient_viewers_patched

    # Patch only defense-catalog-affected strains:
    python3 patch_viewers.py \
        --viewer_dir /path/to/patient_viewers \
        --out_dir    /path/to/patient_viewers_patched \
        --affected_only
"""
import re, json, argparse, shutil, os
from pathlib import Path
from collections import defaultdict, Counter

AFFECTED = {
    'GD172', 'GD177', 'GD18', 'GD205', 'GD206',
    'GD243', 'GD250', 'GD255', 'GD281A', 'GD281B',
    'GD296', 'GD414'
}

AGE_RANK  = {'very_recent': 0, 'recent': 1, 'moderate': 2, 'old': 3, 'nan': 4}
CONF_RANK = {'high': 0, 'moderate': 1, 'low': 2, 'nan': 3}

def best_age(ages):
    return min(ages, key=lambda a: AGE_RANK.get(str(a).lower().strip(), 4))

def best_conf(confs):
    return min(confs, key=lambda c: CONF_RANK.get(str(c).lower().strip(), 3))

def merge_evidence(evs):
    parts = set()
    for e in evs:
        for p in str(e).split(';'):
            p = p.strip()
            if p and p != 'nan':
                parts.add(p)
    return '; '.join(sorted(parts))

def has_recent_cai(entries):
    for e in entries:
        ev  = str(e.get('evidence', '') or e.get('island_evidence', ''))
        age = str(e.get('age', '') or e.get('age_estimate', ''))
        if 'cai' in ev.lower() and age.lower() in ('recent', 'very_recent'):
            return True
    return False

def deduplicate_df(systems):
    """Collapse DefenseFinder systems by type+subtype."""
    groups = defaultdict(list)
    for s in systems:
        key = (s.get('type', ''), s.get('subtype', ''))
        groups[key].append(s)
    out = []
    for key, grp in groups.items():
        if len(grp) == 1:
            out.append(grp[0])
        else:
            # Keep entry furthest from contig edge (start > 0)
            best = max(grp, key=lambda x: x.get('start', 0))
            merged = dict(best)
            merged['n_contig_fragments'] = len(grp)
            out.append(merged)
    return out

def deduplicate_padloc(systems):
    """Collapse PADLOC systems by type+subtype."""
    groups = defaultdict(list)
    for s in systems:
        key = (s.get('type', ''), s.get('subtype', ''))
        groups[key].append(s)
    out = []
    for key, grp in groups.items():
        if len(grp) == 1:
            out.append(grp[0])
        else:
            best = max(grp, key=lambda x: len(str(x.get('proteins', '')).split(',')))
            merged = dict(best)
            merged['n_contig_fragments'] = len(grp)
            out.append(merged)
    return out

def deduplicate_catalog_islands(islands):
    """Collapse catalog_islands by group_id.
    
    - Prefers 'unique' status entries over 'trimmed' (contig-edge fragments)
    - Takes best (most recent) age from any fragment
    - Takes highest confidence from any fragment
    - Merges evidence strings
    - Takes max length, max n_genes, max n_defense, min cai_ratio
    - Adds has_recent_cai and n_contig_fragments fields
    """
    groups = defaultdict(list)
    ungrouped = []
    for isl in islands:
        gid = isl.get('group_id', '').strip()
        if gid:
            groups[gid].append(isl)
        else:
            ungrouped.append(isl)

    out = list(ungrouped)
    for gid, grp in groups.items():
        if len(grp) == 1:
            isl = dict(grp[0])
        else:
            # Prefer unique over trimmed as base
            unique_entries  = [g for g in grp if g.get('status', '') == 'unique']
            trimmed_entries = [g for g in grp if g.get('status', '') == 'trimmed']
            candidates = unique_entries if unique_entries else grp
            base = max(candidates, key=lambda x: x.get('length', 0))
            isl = dict(base)

            # Aggregate across ALL fragments
            isl['age']             = best_age([g.get('age', 'old') for g in grp])
            isl['confidence']      = best_conf([g.get('confidence', 'low') for g in grp])
            isl['evidence']        = merge_evidence([g.get('evidence', '') for g in grp])
            isl['n_evidence']      = len([e for e in isl['evidence'].split(';') if e.strip()])
            isl['length']          = max(g.get('length', 0) for g in grp)
            isl['cai']             = min(float(g.get('cai', 1.0)) for g in grp)
            isl['n_defense']       = max(g.get('n_defense', 0) for g in grp)
            isl['n_mobility']      = max(g.get('n_mobility', 0) for g in grp)
            isl['n_genes']         = max(g.get('n_genes', 0) for g in grp)
            isl['status']          = 'unique' if unique_entries else 'trimmed'
            isl['n_contig_fragments']  = len(grp)
            isl['n_trimmed_fragments'] = len(trimmed_entries)

        isl['has_recent_cai'] = has_recent_cai(grp if len(grp) > 1 else [isl])
        out.append(isl)

    return sorted(out, key=lambda x: x.get('start', 0))

def patch_html(html_path, out_path):
    """Patch a single viewer HTML. Returns (changed, islands_before, islands_after)."""
    with open(html_path) as f:
        content = f.read()

    match = re.search(r'(const STRAIN_DATA = )(\{.*?\})(;\s*\n)', content, re.DOTALL)
    if not match:
        return False, 0, 0

    prefix, json_str, suffix = match.group(1), match.group(2), match.group(3)
    try:
        data = json.loads(json_str)
    except json.JSONDecodeError as e:
        print(f"  ERROR parsing JSON: {e}")
        return False, 0, 0

    changed = False
    islands_before = 0
    islands_after  = 0

    for strain_key in data.keys():
        sd = data[strain_key]

        # Deduplicate df
        if 'df' in sd:
            before = len(sd['df'])
            sd['df'] = deduplicate_df(sd['df'])
            if len(sd['df']) < before:
                changed = True

        # Deduplicate padloc
        if 'padloc' in sd:
            before = len(sd['padloc'])
            sd['padloc'] = deduplicate_padloc(sd['padloc'])
            if len(sd['padloc']) < before:
                changed = True

        # Deduplicate + annotate catalog_islands
        if 'catalog_islands' in sd:
            islands_before += len(sd['catalog_islands'])
            sd['catalog_islands'] = deduplicate_catalog_islands(sd['catalog_islands'])
            islands_after += len(sd['catalog_islands'])
            changed = True  # always true — has_recent_cai added even if no dups

    if not changed:
        return False, 0, 0

    new_json    = json.dumps(data, separators=(',', ':'))
    new_content = content[:match.start()] + prefix + new_json + suffix + content[match.end():]

    with open(out_path, 'w') as f:
        f.write(new_content)

    return True, islands_before, islands_after

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--viewer_dir',    required=True)
    parser.add_argument('--out_dir',       default=None)
    parser.add_argument('--affected_only', action='store_true',
                        help='Only patch the 12 defense-catalog-affected strains')
    args = parser.parse_args()

    viewer_dir = Path(args.viewer_dir)
    out_dir    = Path(args.out_dir) if args.out_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    html_files = sorted(viewer_dir.glob('*.html'))
    print(f"Found {len(html_files)} HTML files")
    print(f"Mode: {'affected strains only' if args.affected_only else 'ALL viewers'}")
    print()

    patched = skipped = total_removed = 0

    for html_path in html_files:
        strain = html_path.stem.replace('_viewer', '').replace('_drilldown', '')

        if args.affected_only and strain not in AFFECTED:
            skipped += 1
            continue

        if not out_dir:
            backup = html_path.with_suffix('.html.bak')
            if not backup.exists():
                shutil.copy2(html_path, backup)

        dest = (out_dir / html_path.name) if out_dir else html_path
        ok, before, after = patch_html(html_path, dest)

        if ok:
            removed = before - after
            total_removed += removed
            if removed > 0:
                print(f"  {strain}: {before} islands -> {after} (-{removed} fragments)")
            patched += 1

    print(f"\n{'='*50}")
    print(f"Viewers patched:          {patched}")
    print(f"Viewers skipped:          {skipped}")
    print(f"Island fragments removed: {total_removed}")
    print(f"has_recent_cai added to all patched viewers")

if __name__ == '__main__':
    main()
