#!/usr/bin/env python3
"""
Inject EOP data into patient HTML viewers and add an EOP heatmap panel.

For each viewer:
  - Adds eop_data to each isolate's STRAIN_DATA entry
  - Adds an EOP heatmap panel below the genome track
  - Shows each phage as a colored bar: green=productive, yellow=intermediate,
    red=resistant, grey=not tested

Usage:
    python3 add_eop_to_viewers.py \
        --viewer_dir /path/to/patient_viewers \
        --eop_csv    /path/to/EOP.csv \
        --out_dir    /path/to/output
"""
import re, json, argparse
import pandas as pd
import numpy as np
from pathlib import Path

# ── Phage groups ──────────────────────────────────────────────────────────────
PHAGE_GROUPS = {
    'Muddy':                'Muddy',
    'Muddy_HRM-GD04':       'Muddy',
    'Muddy_HRM-N0052':      'Muddy',
    'Maco6':                'Muddy',
    'BPsdel33HTH_HRM10':    'BPs',
    'BPsdel33HTH_HRM-GD03': 'BPs',
    'BPsREM1del33':         'BPs',
    'BPsHRM10_pMC09':       'BPs',
    'ZoeJdel43-45':         'ZoeJ',
    'Fionnbharthdel45del47':'ZoeJ',
    'Adephagiadel43del45':  'ZoeJ',
    'CrimD del del':        'other',
    'Itos':                 'other',
    'D29':                  'other',
    'TM4 del rep':          'other',
    'Wildcat':              'other',
    'MissWhite':            'other',
    'Faith1':               'other',
    'Faith1GD69':           'other',
}

# EOP panel HTML — injected once after the canvas
EOP_PANEL_CSS = '''
  /* EOP heatmap panel */
  .eop-panel{margin:8px 0 0 0;background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:10px 14px;font-size:11px}
  .eop-panel-title{font-size:11px;font-weight:600;color:var(--muted);font-family:var(--font-mono);text-transform:uppercase;letter-spacing:.07em;margin-bottom:8px}
  .eop-grid{display:grid;grid-template-columns:90px 1fr;gap:3px 6px;align-items:center}
  .eop-group-label{font-family:var(--font-mono);font-size:10px;color:var(--muted);text-align:right;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .eop-bars{display:flex;flex-wrap:wrap;gap:2px}
  .eop-bar{height:16px;min-width:28px;border-radius:3px;display:flex;align-items:center;justify-content:center;font-size:9px;font-family:var(--font-mono);font-weight:600;color:#fff;cursor:default;position:relative}
  .eop-bar:hover .eop-tooltip{display:block}
  .eop-tooltip{display:none;position:absolute;bottom:calc(100% + 4px);left:50%;transform:translateX(-50%);background:#1a1a2e;color:#fff;font-size:10px;padding:4px 8px;border-radius:4px;white-space:nowrap;z-index:100;pointer-events:none}
  .eop-no-data{font-size:10px;color:var(--muted);font-style:italic}
'''

EOP_PANEL_HTML = '''
<div class="eop-panel" id="eopPanel">
  <div class="eop-panel-title">Phage infection (EOP log₁₀)</div>
  <div id="eopContent"><span class="eop-no-data">Select an isolate to view EOP data</span></div>
</div>
'''

EOP_JS = '''
// ── EOP heatmap ───────────────────────────────────────────────────────────────
const EOP_GROUP_ORDER = ['Muddy', 'BPs', 'ZoeJ', 'other'];
const EOP_GROUP_COLORS = {
  'Muddy': '#185FA5',
  'BPs':   '#534AB7',
  'ZoeJ':  '#0e7a5a',
  'other': '#888780',
};

function eopColor(eop) {
  if (eop === null || eop === undefined) return '#d0d0d0';  // not tested
  if (eop >= -1)  return '#1a9e6e';  // productive (bright green)
  if (eop >= -3)  return '#7ecba1';  // low productive (light green)
  if (eop >= -5)  return '#f59e0b';  // intermediate (amber)
  if (eop >= -7)  return '#ef4444';  // resistant (red)
  return '#991b1b';                  // highly resistant (dark red)
}

function eopLabel(eop) {
  if (eop === null || eop === undefined) return 'NT';
  if (eop === 0)  return '0';
  return eop.toString();
}

function renderEopPanel(isolate) {
  const panel = document.getElementById('eopContent');
  if (!panel) return;

  const data = (STRAIN_DATA[isolate] || {}).eop_data || [];
  if (!data.length) {
    panel.innerHTML = '<span class="eop-no-data">No EOP data for this isolate</span>';
    return;
  }

  // Group by phage family
  const groups = {};
  EOP_GROUP_ORDER.forEach(g => groups[g] = []);
  data.forEach(d => {
    const g = d.group || 'other';
    if (!groups[g]) groups[g] = [];
    groups[g].push(d);
  });

  let html = '<div class="eop-grid">';
  EOP_GROUP_ORDER.forEach(grpName => {
    const entries = groups[grpName] || [];
    if (!entries.length) return;
    const color = EOP_GROUP_COLORS[grpName] || '#888780';
    html += `<div class="eop-group-label" style="color:${color};font-weight:600">${grpName}</div>`;
    html += '<div class="eop-bars">';
    entries.sort((a, b) => (b.eop ?? -9) - (a.eop ?? -9));
    entries.forEach(e => {
      const bg  = eopColor(e.eop);
      const lbl = eopLabel(e.eop);
      html += `<div class="eop-bar" style="background:${bg};min-width:${lbl.length > 2 ? 34 : 28}px">
        ${lbl}
        <div class="eop-tooltip">${e.phage}: ${e.eop !== null ? '10^'+e.eop : 'NT'}</div>
      </div>`;
    });
    html += '</div>';
  });
  html += '</div>';
  panel.innerHTML = html;
}
'''

def to_log10_eop(v):
    try:
        f = float(str(v).strip())
        return round(np.log10(f), 1) if f > 0 else None
    except:
        return None

def build_eop_lookup(eop_csv):
    df = pd.read_csv(eop_csv, header=None)
    strain_row = df.iloc[1].astype(str)
    strains = list(strain_row.values[1:])

    eop_by_isolate = {}
    skip_rows = {6, 7, 18}

    for row_idx in range(8, len(df)):
        phage = str(df.iloc[row_idx, 0]).strip()
        if not phage or phage == 'nan' or row_idx in skip_rows:
            continue
        group = PHAGE_GROUPS.get(phage)
        if group is None:
            continue
        for strain, v in zip(strains, df.iloc[row_idx].values[1:]):
            eop = to_log10_eop(v)
            if eop is not None:
                s = str(strain).strip()
                eop_by_isolate.setdefault(s, []).append({
                    'phage': phage,
                    'group': group,
                    'eop':   float(eop)
                })
    return eop_by_isolate

def patch_viewer(content, eop_by_isolate):
    if 'eopPanel' in content:
        return content, False  # already patched

    # ── 1. Inject EOP CSS into <style> block ─────────────────────────────────
    style_end = content.rfind('</style>')
    if style_end == -1:
        return content, False
    content = content[:style_end] + EOP_PANEL_CSS + '\n</style>' + content[style_end+8:]

    # ── 2. Add EOP data to each isolate in STRAIN_DATA ───────────────────────
    m = re.search(r'(const STRAIN_DATA = )(\{.*?\})(;\s*\n)', content, re.DOTALL)
    if not m:
        return content, False
    data = json.loads(m.group(2))

    injected = 0
    for isolate_key in data:
        eop_entries = eop_by_isolate.get(isolate_key, [])
        data[isolate_key]['eop_data'] = eop_entries
        if eop_entries:
            injected += 1

    new_json = json.dumps(data, separators=(',', ':'))
    content = content[:m.start()] + m.group(1) + new_json + m.group(3) + content[m.end():]

    # ── 3. Inject EOP JS at top-level scope (after PATIENT_GROUPS) ───────────
    pg_match = re.search(r'const PATIENT_GROUPS = \{.*?\};\s*\n', content, re.DOTALL)
    if pg_match:
        pos = pg_match.end()
        content = content[:pos] + EOP_JS + content[pos:]

    # ── 4. Inject EOP panel HTML after canvas / before patient-summary ────────
    anchor = '<div class="patient-summary"'
    if anchor in content:
        content = content.replace(anchor, EOP_PANEL_HTML + '\n' + anchor, 1)

    # ── 5. Hook renderEopPanel into isolate selection ─────────────────────────
    # Find where individual isolate is rendered (when isolate button clicked)
    # The viewer calls renderPatient(pid) which calls renderIsolate(isolate)
    # We need to hook into renderIsolate or wherever single isolate is selected
    old_render_hook = 'function renderPatient(pid) {'
    new_render_hook = '''function renderPatient(pid) {
  // Render EOP for first isolate of this patient
  const firstIsolate = (PATIENT_GROUPS[pid] || [])[0];
  if (firstIsolate) renderEopPanel(firstIsolate);'''

    if old_render_hook in content and 'renderEopPanel' not in content.split('function renderPatient')[0]:
        content = content.replace(old_render_hook, new_render_hook, 1)

    return content, True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--viewer_dir', required=True)
    parser.add_argument('--eop_csv',    required=True)
    parser.add_argument('--out_dir',    required=True)
    args = parser.parse_args()

    print("Building EOP lookup...")
    eop_by_isolate = build_eop_lookup(args.eop_csv)
    print(f"  {len(eop_by_isolate)} isolates with EOP data")

    viewer_dir = Path(args.viewer_dir)
    out_dir    = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    html_files = sorted(viewer_dir.glob('*.html'))
    print(f"Found {len(html_files)} HTML files")

    patched = matched = 0
    for html_path in html_files:
        with open(html_path) as f:
            content = f.read()

        new_content, changed = patch_viewer(content, eop_by_isolate)

        out_path = out_dir / html_path.name
        with open(out_path, 'w') as f:
            f.write(new_content)

        if changed:
            patched += 1

    print(f"\nDone.")
    print(f"  Viewers patched: {patched}/{len(html_files)}")

if __name__ == '__main__':
    main()
