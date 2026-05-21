#!/usr/bin/env python3
"""
Update all patient HTML viewers:
  1. Remove the legend div (depth-aware)
  2. Add GI confidence toggle to controls (default: high only)
  3. Inject GI_CONF_LEVEL at top-level script scope (after PATIENT_GROUPS)
  4. Update the catIslands filter to respect confidence level

Usage:
    python3 update_viewer_ui.py \
        --viewer_dir /path/to/patient_viewers \
        --out_dir    /path/to/patient_viewers_updated
"""
import re, argparse
from pathlib import Path

# ── Confidence toggle HTML ────────────────────────────────────────────────────
CONF_TOGGLE_HTML = '''\n  <div style="display:flex;align-items:center;gap:8px;margin-left:8px">
    <span style="font-size:11px;color:var(--muted);font-family:var(--font-mono);text-transform:uppercase;letter-spacing:.07em;white-space:nowrap">GI confidence:</span>
    <div class="toggle-btn" id="tog-conf-mod"
         onclick="toggleConfidence()"
         style="background:#fff3cc22;border-color:#8a5c00;cursor:pointer">
      <div class="toggle-dot" style="background:#8a5c00"></div>
      <span id="conf-label">High only</span>
    </div>
  </div>'''

# ── Top-level JS — goes right after PATIENT_GROUPS declaration ───────────────
CONF_JS_TOPLEVEL = '''
// ── GI confidence filter (top-level state) ───────────────────────────────────
let GI_CONF_LEVEL = 'high'; // 'high' = high only; 'moderate' = high + moderate

function toggleConfidence() {
  const btn = document.getElementById('tog-conf-mod');
  const lbl = document.getElementById('conf-label');
  if (GI_CONF_LEVEL === 'high') {
    GI_CONF_LEVEL = 'moderate';
    btn.style.background = '#fff3cc99';
    lbl.textContent = 'High + Moderate';
  } else {
    GI_CONF_LEVEL = 'high';
    btn.style.background = '#fff3cc22';
    lbl.textContent = 'High only';
  }
  const pid = document.getElementById('patientSelect').value;
  if (pid) redrawAll(pid);
}
'''

def remove_legend(content):
    """Remove <div class="legend">...</div> using depth-counting."""
    start = content.find('<div class="legend">')
    if start == -1:
        return content, False
    depth = 0
    pos = start
    end = -1
    while pos < len(content):
        if content[pos:pos+4] == '<div':
            depth += 1
        elif content[pos:pos+6] == '</div>':
            depth -= 1
            if depth == 0:
                end = pos + 6
                break
        pos += 1
    if end == -1:
        return content, False
    new_content = content[:start].rstrip('\n') + '\n' + content[end:].lstrip('\n')
    return new_content, True

def add_confidence_toggle(content):
    """Add toggle button to controls bar, before closing </div> of controls."""
    if 'tog-conf-mod' in content:
        return content, False  # already added

    # Insert before <div class="patient-summary"
    anchor = '<div class="patient-summary"'
    if anchor not in content:
        return content, False

    insert_at = content.find(anchor)
    # Walk back to find the </div> that closes the controls div
    close_pos = content.rfind('</div>', 0, insert_at)
    if close_pos == -1:
        return content, False

    new_content = (
        content[:close_pos] +
        CONF_TOGGLE_HTML +
        '\n</div>\n\n' +
        content[insert_at:]
    )
    return new_content, True

def inject_conf_js(content):
    """Inject GI_CONF_LEVEL and toggleConfidence at top-level scope.
    
    Injects after the PATIENT_GROUPS line — guaranteed top-level, before any functions.
    """
    if 'GI_CONF_LEVEL' in content:
        return content, False  # already injected

    # Find end of PATIENT_GROUPS line
    pg_match = re.search(r'const PATIENT_GROUPS = \{.*?\};\s*\n', content, re.DOTALL)
    if not pg_match:
        # Fallback: inject after SPOT_COLORS or DEFENSE_COLORS
        pg_match = re.search(r'const SPOT_COLORS\s*=.*?;\s*\n', content, re.DOTALL)
    if not pg_match:
        return content, False

    inject_pos = pg_match.end()
    new_content = content[:inject_pos] + CONF_JS_TOPLEVEL + content[inject_pos:]
    return new_content, True

def update_catislands_filter(content):
    """Update the catIslands filter to check confidence level."""
    if 'GI_CONF_LEVEL' not in content:
        return content, False  # JS not injected yet

    # Target: the specific catIslands filter line
    old = (
        'const catIslands = (data.catalog_islands || [])\n'
        '        .filter(i => inView(i.start, i.end))\n'
        '        .sort((a,b) => a.start - b.start);'
    )
    new = (
        'const catIslands = (data.catalog_islands || [])\n'
        '        .filter(i => {\n'
        '          if (!inView(i.start, i.end)) return false;\n'
        '          const conf = (i.confidence || \'\').toLowerCase();\n'
        '          if (GI_CONF_LEVEL === \'high\') return conf === \'high\';\n'
        '          return conf === \'high\' || conf === \'moderate\';\n'
        '        })\n'
        '        .sort((a,b) => a.start - b.start);'
    )
    if old not in content:
        return content, False

    return content.replace(old, new, 1), True

def patch_viewer(content):
    changed = False
    content, c = remove_legend(content);         changed = changed or c
    content, c = add_confidence_toggle(content); changed = changed or c
    content, c = inject_conf_js(content);        changed = changed or c
    content, c = update_catislands_filter(content); changed = changed or c
    return content, changed

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--viewer_dir', required=True)
    parser.add_argument('--out_dir',    required=True)
    args = parser.parse_args()

    viewer_dir = Path(args.viewer_dir)
    out_dir    = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    html_files = sorted(viewer_dir.glob('*.html'))
    print(f"Found {len(html_files)} HTML files")

    counts = dict(patched=0, legend=0, toggle=0, js=0, filter=0)

    for html_path in html_files:
        with open(html_path) as f:
            content = f.read()

        had = {
            'legend': '<div class="legend">' in content,
            'toggle': 'tog-conf-mod' not in content,
            'js':     'GI_CONF_LEVEL' not in content,
            'filter': 'filter(i => inView(i.start, i.end))' in content,
        }

        new_content, changed = patch_viewer(content)

        out_path = out_dir / html_path.name
        with open(out_path, 'w') as f:
            f.write(new_content)

        if changed: counts['patched'] += 1
        if had['legend'] and '<div class="legend">' not in new_content: counts['legend'] += 1
        if had['toggle'] and 'tog-conf-mod' in new_content:             counts['toggle'] += 1
        if had['js']     and 'GI_CONF_LEVEL' in new_content:            counts['js']     += 1
        if had['filter'] and 'filter(i => inView' not in new_content:   counts['filter'] += 1

    print(f"\nDone.")
    print(f"  Files processed:      {len(html_files)}")
    print(f"  Files changed:        {counts['patched']}")
    print(f"  Legend removed:       {counts['legend']}")
    print(f"  Toggle added:         {counts['toggle']}")
    print(f"  Top-level JS added:   {counts['js']}")
    print(f"  catIslands updated:   {counts['filter']}")

if __name__ == '__main__':
    main()
