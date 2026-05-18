#!/usr/bin/env python3
"""
generate_patient_html.py

Generates a self-contained single-patient HTML viewer with all data
embedded inline. No server, no viewer_data/ folder, no fetch calls.
Opens in any browser from any location.

Usage:
    python generate_patient_html.py \
        --patient   GD233 \
        --out       GD233_viewer.html \
        --data_dir  mabs_island_results/viewer_data \
        --html_dir  mabs_island_results

    # Generate all patients at once
    python generate_patient_html.py \
        --all \
        --out_dir   patient_viewers \
        --data_dir  mabs_island_results/viewer_data \
        --html_dir  mabs_island_results
"""

import os, re, json, argparse, csv


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--patient",  default=None,
                   help="Patient ID to generate (e.g. GD233)")
    p.add_argument("--all",      action="store_true",
                   help="Generate HTML for every patient")
    p.add_argument("--out",      default=None,
                   help="Output HTML file (single patient mode)")
    p.add_argument("--out_dir",  default="patient_viewers",
                   help="Output directory (--all mode, default: patient_viewers)")
    p.add_argument("--data_dir", default="mabs_island_results/viewer_data",
                   help="Directory of per-strain JSON files")
    p.add_argument("--html_dir", default="mabs_island_results",
                   help="Directory containing patient_comparison_viewer.html")
    p.add_argument("--template", default=None,
                   help="Path to patient_comparison_viewer.html (overrides --html_dir)")
    return p.parse_args()


def canonical(s):
    return re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', s, flags=re.IGNORECASE)


def load_patient_groups(html_content):
    """Extract PATIENT_GROUPS from the viewer HTML."""
    m = re.search(r'const PATIENT_GROUPS\s*=\s*(\{.*?\});', html_content,
                  re.DOTALL)
    if not m:
        raise ValueError("Could not find PATIENT_GROUPS in HTML template")
    return json.loads(m.group(1))


def load_strain_data(strain, data_dir):
    """Load per-strain JSON. Returns {} if not found."""
    # Try exact name first
    path = os.path.join(data_dir, f"{strain}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    # Try canonical
    canon = canonical(strain)
    path2 = os.path.join(data_dir, f"{canon}.json")
    if os.path.exists(path2):
        with open(path2) as f:
            return json.load(f)
    # Try scanning for any file that canonicalises to same name
    if os.path.exists(data_dir):
        for fname in os.listdir(data_dir):
            if canonical(fname.replace(".json","")) == canon:
                with open(os.path.join(data_dir, fname)) as f:
                    return json.load(f)
    return {
        "length": 5000000,
        "rgps": [], "scanner_islands": [], "fixed_defense": [],
        "catalog_islands": [], "df": [], "padloc": [], "dp": [],
        "prophage": [], "asm": []
    }


def load_fixed_blast_annotations(data_dir):
    """Load optional island/source BLAST annotations from the results folder."""
    base_dir = os.path.dirname(os.path.abspath(data_dir))
    path = os.path.join(base_dir, "fixed_defense_catalog_final.tsv")
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def merge_unique(values):
    out = []
    seen = set()
    for value in values:
        value = (value or "").strip()
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def add_blast_annotations(strain, data, blast_rows):
    """Attach matched BLAST/source calls to catalog islands by coordinate overlap."""
    islands = data.get("catalog_islands") or []
    if not islands or not blast_rows:
        return data

    strain_keys = {strain, canonical(strain)}
    rows = [r for r in blast_rows if r.get("strain") in strain_keys or canonical(r.get("strain", "")) in strain_keys]
    if not rows:
        return data

    for isl in islands:
        try:
            istart, iend = int(isl.get("start", 0)), int(isl.get("end", 0))
        except (TypeError, ValueError):
            continue
        matches = []
        for row in rows:
            try:
                rstart, rend = int(float(row.get("island_start") or 0)), int(float(row.get("island_end") or 0))
            except (TypeError, ValueError):
                continue
            overlap = max(0, min(iend, rend) - max(istart, rstart))
            if overlap <= 0:
                continue
            denom = max(1, min(iend - istart, rend - rstart))
            if overlap / denom >= 0.45:
                matches.append(row)
        if not matches:
            continue

        isl["blast_systems"] = "; ".join(merge_unique(m.get("subtype") for m in matches))
        isl["blast_hit_species"] = "; ".join(merge_unique(m.get("blast_hit_species") for m in matches))
        isl["blast_validated"] = "Yes" if any(m.get("blast_validated") == "Yes" for m in matches) else "No"
        isl["transfer_mechanism"] = "; ".join(merge_unique(m.get("transfer_mechanism") for m in matches))
        isl["blast_accession"] = "; ".join(merge_unique(m.get("blast_accession") for m in matches))
        isl["blast_description"] = "; ".join(merge_unique(m.get("blast_description") for m in matches))
        isl["blast_evalue"] = "; ".join(merge_unique(m.get("blast_evalue") for m in matches))
        isl["source_subspecies"] = "; ".join(merge_unique(m.get("subspecies") for m in matches))
    return data


SAFARI_CANVAS_POLYFILL = r"""
// Safari compatibility: older Safari versions do not implement
// CanvasRenderingContext2D.roundRect(), which otherwise stops drawing.
if (!CanvasRenderingContext2D.prototype.roundRect) {
  CanvasRenderingContext2D.prototype.roundRect = function(x, y, w, h, r) {
    r = Math.min(r || 0, Math.abs(w) / 2, Math.abs(h) / 2);
    this.moveTo(x + r, y);
    this.lineTo(x + w - r, y);
    this.quadraticCurveTo(x + w, y, x + w, y + r);
    this.lineTo(x + w, y + h - r);
    this.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
    this.lineTo(x + r, y + h);
    this.quadraticCurveTo(x, y + h, x, y + h - r);
    this.lineTo(x, y + r);
    this.quadraticCurveTo(x, y, x + r, y);
    this.closePath();
  };
}
"""


GI_TABLE_CSS = r"""

  .gi-table-section{border-top:1px solid var(--border);padding:.75rem 1.25rem 1rem;background:#fff}
  .gi-table-head{display:flex;align-items:center;justify-content:space-between;gap:1rem;margin-bottom:.5rem}
  .gi-table-title{font-family:var(--font-mono);font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}
  .gi-table-count{font-family:var(--font-mono);font-size:11px;color:var(--muted)}
  .gi-table-controls{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-bottom:.55rem}
  .gi-table-controls select,.gi-table-controls input{font-family:var(--font-mono);font-size:11px;border:1px solid var(--border);border-radius:4px;background:#f8fafd;color:#283149;padding:5px 8px}
  .gi-table-controls select{min-width:180px}
  .gi-table-controls input{min-width:220px;max-width:340px}
  .gi-filter-chip{font-family:var(--font-mono);font-size:10px;color:#5a6380;background:#eaecf0;border-radius:3px;padding:3px 6px}
  .gi-table-scroll{max-height:360px;overflow:auto;border:1px solid var(--border);border-radius:6px;background:#fff}
  table.gi-table{width:100%;border-collapse:collapse;font-family:var(--font-mono);font-size:11px;min-width:1180px}
  .gi-table th{position:sticky;top:0;background:#f8fafd;color:#4d5874;text-align:left;font-weight:600;border-bottom:1px solid var(--border);padding:7px 8px;z-index:1}
  .gi-table td{border-bottom:1px solid #edf0f5;padding:6px 8px;vertical-align:top;color:#283149}
  .gi-table tr:last-child td{border-bottom:0}
  .gi-table .num{text-align:right;white-space:nowrap}
  .gi-table .coords{white-space:nowrap}
  .gi-table .muted{color:#7a8297}
  .gi-pill{display:inline-block;border-radius:3px;padding:1px 5px;font-weight:600;white-space:nowrap}
  .gi-pill.high{background:#d4f2e8;color:#0e7a5a}
  .gi-pill.moderate{background:#fff3cc;color:#8a5c00}
  .gi-pill.low,.gi-pill.very_low{background:#eaecf0;color:#5a6380}
  .gi-pill.hgt{background:#fde8e4;color:#a82c1a}
  .gi-pill.within{background:#eaecf0;color:#5a6380}
  .gi-table .wide{min-width:180px}
"""


GI_TABLE_JS = r"""

function htmlEscape(value) {
  return String(value == null || value === '' ? '—' : value)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')
    .replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}
function fmtBp(value) {
  const n = Number(value || 0);
  return n ? n.toLocaleString() : '0';
}
function fmtKb(value) {
  const n = Number(value || 0);
  return n ? (n/1000).toFixed(n >= 100000 ? 0 : 1) : '0';
}
function compactText(value, maxLen) {
  const s = String(value || '').trim();
  if (!s) return '—';
  return s.length > maxLen ? s.slice(0, maxLen - 1) + '…' : s;
}
function inferSource(isl) {
  return isl.blast_hit_species || isl.source_species || isl.source_subspecies ||
         isl.blast_taxon || isl.top_blast_hit || isl.best_hit_species || '';
}
function inferTransfer(isl) {
  const value = isl.transfer_mechanism || isl.source_class || isl.origin_class || '';
  if (value) return value;
  const source = inferSource(isl).toLowerCase();
  const cargo = String(isl.dominant_cargo || '').toLowerCase();
  const evidence = String(isl.evidence || '').toLowerCase();
  if (source.indexOf('phage') >= 0 || cargo.indexOf('phage') >= 0) return 'phage/prophage-like';
  if (source && source.indexOf('abscessus only') < 0) return 'possible HGT';
  if (evidence.indexOf('mobility_gene') >= 0) return 'mobile element';
  return '';
}
function islandSearchText(isl, source, transfer) {
  return [
    isl.island_id, isl.group_id, isl.start, isl.end, isl.length, isl.confidence,
    isl.age || isl.age_estimate, isl.dominant_cargo, isl.mob_types, isl.evidence,
    isl.trna_flanked, isl.has_dr, isl.blast_validated, source, transfer,
    isl.blast_systems, isl.blast_description, isl.blast_accession
  ].join(' ').toLowerCase();
}
function applyIslandTableFilter(sectionId) {
  const section = document.getElementById(sectionId);
  if (!section) return;
  const mode = section.querySelector('.gi-filter-mode').value;
  const query = section.querySelector('.gi-filter-query').value.trim().toLowerCase();
  const rows = Array.from(section.querySelectorAll('tbody tr'));
  let visibleCount = 0;
  rows.forEach(row => {
    let show = true;
    if (mode === 'high') show = row.dataset.confidence === 'high';
    else if (mode === 'hgt') show = row.dataset.transfer.indexOf('hgt') >= 0;
    else if (mode === 'external_hgt') show = row.dataset.transfer.indexOf('hgt') >= 0 && row.dataset.source.indexOf('abscessus only') < 0;
    else if (mode === 'blast') show = row.dataset.blast === 'yes';
    else if (mode === 'phage') show = row.dataset.source.indexOf('phage') >= 0 || row.dataset.cargo.indexOf('phage') >= 0 || row.dataset.transfer.indexOf('phage') >= 0;
    else if (mode === 'defense') show = row.dataset.defense === 'true';
    else if (mode === 'mobility') show = row.dataset.mobility === 'true';
    else if (mode === 'recent') show = row.dataset.age.indexOf('recent') >= 0;
    if (show && query) show = row.dataset.search.indexOf(query) >= 0;
    row.style.display = show ? '' : 'none';
    if (show) visibleCount++;
  });
  const count = section.querySelector('.gi-table-count');
  if (count) count.textContent = visibleCount + ' of ' + rows.length + ' islands';
}
function renderIslandTable(strain) {
  const data = STRAIN_DATA[strain] || {};
  const islands = (data.catalog_islands || []).slice().sort((a,b)=>(a.start||0)-(b.start||0));
  if (!islands.length) {
    return '<div class="gi-table-section"><div class="gi-table-head"><div class="gi-table-title">Genomic island summary</div></div><div class="no-data">No catalog genomic islands for this isolate</div></div>';
  }
  const sectionId = 'gi_table_' + strain.replace(/[^a-zA-Z0-9]/g,'_');
  const rows = islands.map(isl => {
    const confidence = htmlEscape(isl.confidence || '');
    const transfer = inferTransfer(isl);
    const transferClass = /hgt/i.test(transfer) ? 'hgt' : 'within';
    const source = inferSource(isl);
    const search = htmlEscape(islandSearchText(isl, source, transfer));
    const defense = Number(isl.n_defense || 0) > 0 || isl.has_defense === true || isl.has_defense === 'true';
    const mobility = Number(isl.n_mobility || 0) > 0 || String(isl.mob_types || '').trim() !== '';
    return `<tr
      data-confidence="${htmlEscape(String(isl.confidence || '').toLowerCase())}"
      data-transfer="${htmlEscape(String(transfer || '').toLowerCase())}"
      data-source="${htmlEscape(String(source || '').toLowerCase())}"
      data-cargo="${htmlEscape(String(isl.dominant_cargo || '').toLowerCase())}"
      data-age="${htmlEscape(String(isl.age || isl.age_estimate || '').toLowerCase())}"
      data-blast="${htmlEscape(String(isl.blast_validated || '').toLowerCase())}"
      data-defense="${defense ? 'true' : 'false'}"
      data-mobility="${mobility ? 'true' : 'false'}"
      data-search="${search}">
      <td>${htmlEscape(isl.island_id || '')}</td>
      <td>${htmlEscape(isl.group_id || '')}</td>
      <td class="coords">${fmtBp(isl.start)}–${fmtBp(isl.end)}</td>
      <td class="num">${fmtKb(isl.length || ((isl.end||0)-(isl.start||0)))} kb</td>
      <td><span class="gi-pill ${confidence}">${confidence}</span></td>
      <td>${htmlEscape(isl.age || isl.age_estimate || '')}</td>
      <td>${htmlEscape(isl.dominant_cargo || '')}</td>
      <td class="num">${htmlEscape(isl.n_genes || '')}</td>
      <td class="num">${htmlEscape(isl.n_defense || 0)}</td>
      <td class="num">${htmlEscape(isl.n_mobility || 0)}</td>
      <td>${htmlEscape(isl.mob_types || '')}</td>
      <td>${htmlEscape(isl.n_evidence || 0)} <span class="muted">${htmlEscape(compactText(isl.evidence, 42))}</span></td>
      <td>${htmlEscape(isl.trna_flanked || '')}</td>
      <td>${htmlEscape(isl.has_dr || '')}</td>
      <td>${htmlEscape(isl.blast_validated || '')}</td>
      <td class="wide">${htmlEscape(compactText(source, 70))}</td>
      <td>${transfer ? `<span class="gi-pill ${transferClass}">${htmlEscape(compactText(transfer, 36))}</span>` : '—'}</td>
      <td class="wide">${htmlEscape(compactText(isl.blast_systems || isl.blast_description || '', 80))}</td>
    </tr>`;
  }).join('');
  return `<div class="gi-table-section" id="${sectionId}">
    <div class="gi-table-head">
      <div class="gi-table-title">Genomic island summary</div>
      <div class="gi-table-count">${islands.length} islands</div>
    </div>
    <div class="gi-table-controls">
      <select class="gi-filter-mode" onchange="applyIslandTableFilter('${sectionId}')">
        <option value="all">All islands</option>
        <option value="high">High confidence</option>
        <option value="hgt">Any HGT call</option>
        <option value="external_hgt">External/non-abscessus HGT</option>
        <option value="blast">BLAST validated</option>
        <option value="phage">Phage/prophage-like</option>
        <option value="defense">Defense cargo</option>
        <option value="mobility">Mobility cargo</option>
        <option value="recent">Recent/very recent</option>
      </select>
      <input class="gi-filter-query" type="search" placeholder="Search islands, species, cargo..." oninput="applyIslandTableFilter('${sectionId}')">
      <span class="gi-filter-chip">table only</span>
    </div>
    <div class="gi-table-scroll">
      <table class="gi-table">
        <thead><tr>
          <th>Island</th><th>Group</th><th>Coords</th><th>Length</th><th>Conf</th><th>Age</th>
          <th>Cargo</th><th>Genes</th><th>Defense</th><th>Mobility</th><th>Mob types</th>
          <th>Evidence</th><th>tRNA</th><th>DR</th><th>BLAST</th><th>Likely source</th><th>Transfer</th><th>Source note</th>
        </tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
  </div>`;
}
"""


def make_generated_html_safari_safe(html):
    """Patch template-derived viewer JS for Safari and empty-prophage patients."""
    if ".gi-table-section" not in html:
        html = html.replace("</style>", GI_TABLE_CSS + "\n</style>", 1)
    if "CanvasRenderingContext2D.prototype.roundRect" not in html:
        html = html.replace("<script>\n", "<script>\n" + SAFARI_CANVAS_POLYFILL + "\n", 1)
    if "function renderIslandTable(strain)" not in html:
        html = html.replace("// ── Tooltip", GI_TABLE_JS + "\n// ── Tooltip", 1)

    # In the template, fixed-defense and genomic-island drawing accidentally sit
    # inside the "has prophage rows" block. Patients with no prophage then lose
    # GI/fixed-defense tracks entirely. Keep the prophage drawing conditional,
    # but always compute the downstream row positions.
    html = html.replace(
        "  // Prophage track\n"
        "  if (visible.pp && pp.length>0) {\n"
        "    const PP_Y = GENOME_Y+GENOME_H+(visible.rgp?RGP_H+8:4);",
        "  // Prophage track\n"
        "  const PP_Y = GENOME_Y+GENOME_H+(visible.rgp?RGP_H+8:4);\n"
        "  if (visible.pp && pp.length>0) {"
    )
    html = html.replace(
        "    });\n\n"
        "    // ── Fixed defense islands track",
        "    });\n"
        "  }\n\n"
        "  // ── Fixed defense islands track"
    )
    html = html.replace(
        "    ctx.fillStyle='#7c2d92';ctx.font='9px IBM Plex Mono,monospace';\n"
        "    ctx.textAlign='right';\n"
        "    const PP_Y2 = GENOME_Y+GENOME_H+(visible.rgp?RGP_H+8:4);\n"
        "    ctx.fillText('Phage',PAD-4,PP_Y2+PP_H/2+3);\n"
        "  }",
        "  if (visible.pp) {\n"
        "    ctx.fillStyle='#7c2d92';ctx.font='9px IBM Plex Mono,monospace';\n"
        "    ctx.textAlign='right';\n"
        "    ctx.fillText('Phage',PAD-4,PP_Y+PP_H/2+3);\n"
        "  }"
    )
    html = html.replace(
        "    // ── Nested genomic island catalog tracks ─────────────────────────────────\n"
        "    // Each nesting depth gets its own row. Top-level (depth 0) are full height.\n"
        "    // Nested children are shorter and offset downward, visually inside the parent.\n"
        "    const ISL_ROW_H  = 12;\n"
        "    const ISL_GAP    = 3;\n"
        "    const ISL_BASE_Y = FIX_Y + (visible.fix ? FIX_H + 6 : 0);",
        "  // ── Nested genomic island catalog tracks ─────────────────────────────────\n"
        "  // Each nesting depth gets its own row. Top-level (depth 0) are full height.\n"
        "  // Nested children are shorter and offset downward, visually inside the parent.\n"
        "  const ISL_BASE_Y = FIX_Y + (visible.fix ? FIX_H + 6 : 0);"
    )
    html = html.replace(
        """  // Estimate height
  const dfR  = visible.df     ? Math.min(df.length,5)  : 0;
  const plR  = visible.padloc ? Math.min(pl.length,5)  : 0;
  const dpR  = visible.dp     ? Math.min(dp.length,5)  : 0;
  const ppR  = visible.pp     ? 1 : 0;
  const aboveH = (dfR+plR+dpR)*( ARROW_H+3) + (dfR>0?TRACK_GAP:0)+(plR>0?TRACK_GAP:0)+(dpR>0?TRACK_GAP:0) + 24;
  // Dynamic height: estimate lanes needed from island count and genome length
  const catIslands  = data.catalog_islands || [];
  const ISL_ROW_H   = 12;
  const ISL_GAP     = 3;
  // Estimate max lanes: islands / (genome_length / avg_island_width_px)
  // Conservative: assume up to 4 lanes for dense genomes, 1 for sparse
  const islDensity  = catIslands.length / Math.max(1, (data.length||5000000) / 50000);
  const estLanes    = visible.isl && catIslands.length > 0
    ? Math.min(5, Math.max(1, Math.ceil(islDensity)))
    : 0;
  // Reserve minimum height for island track even when lane count is unknown
  const islTrackH   = visible.isl
    ? Math.max(4 * (ISL_ROW_H + ISL_GAP) + 8, estLanes * (ISL_ROW_H + ISL_GAP) + (estLanes > 0 ? 8 : 0))
    : 0;
  const belowH = (visible.rgp?RGP_H+6:0) + (visible.pp?PP_H+6:0) +
                 (visible.fix?8+4:0) + islTrackH + 24;
  const H = Math.max(120, aboveH + GENOME_H + belowH);""",
        """  // Viewport-dependent layout sizing. Tracks are stacked by the features
  // visible in the current zoom window, so dense GI regions stay inside frame.
  const trackW0 = W - PAD - 8;
  const visibleBp0 = maxLen / vs.zoom;
  const panClamped0 = Math.max(0, Math.min(vs.pan, maxLen - visibleBp0));
  function preInView(start, stop) { return stop >= panClamped0 && start <= panClamped0 + visibleBp0; }
  function preBpToX(bp) { return PAD + (bp - panClamped0) / (visibleBp0 / trackW0); }
  function preBpWidth(bp) { return bp / (visibleBp0 / trackW0); }
  function countPackedLanes(items, startKey, stopKey, minWidth) {
    const laneEnds = [];
    const visibleItems = (items || [])
      .filter(d => preInView(d[startKey], d[stopKey]))
      .sort((a,b) => a[startKey] - b[startKey]);
    visibleItems.forEach(d => {
      const x = preBpToX(d[startKey]);
      const w = Math.max(minWidth, preBpWidth(d[stopKey] - d[startKey]));
      let lane = 0;
      while (laneEnds[lane] !== undefined && laneEnds[lane] > x - 2) lane++;
      laneEnds[lane] = x + w + 2;
    });
    return laneEnds.length;
  }

  const dfR  = visible.df     ? Math.min(countPackedLanes(df, 'start', 'stop', 4), 5) : 0;
  const plR  = visible.padloc ? Math.min(countPackedLanes(pl, 'start', 'stop', 4), 5) : 0;
  const dpR  = visible.dp     ? Math.min(countPackedLanes(dp, 'start', 'stop', 4), 5) : 0;
  const aboveH = (dfR+plR+dpR)*(ARROW_H+3) + (dfR>0?TRACK_GAP:0)+(plR>0?TRACK_GAP:0)+(dpR>0?TRACK_GAP:0) + 24;

  const catIslands  = data.catalog_islands || [];
  const ISL_ROW_H   = 12;
  const ISL_GAP     = 3;
  const GI_LANE_LIMIT = vs.zoom < 1.5 ? 6 : vs.zoom < 3 ? 10 : vs.zoom < 8 ? 16 : 28;
  const rawIslandLanes = visible.isl ? countPackedLanes(catIslands, 'start', 'end', 4) : 0;
  const visibleIslandLanes = Math.min(rawIslandLanes, GI_LANE_LIMIT);
  const islTrackH = visible.isl && visibleIslandLanes > 0
    ? visibleIslandLanes * (ISL_ROW_H + ISL_GAP) + 12
    : 0;
  const belowH = (visible.rgp?RGP_H+6:0) + (visible.pp?PP_H+6:0) +
                 (visible.fix?8+4:0) + islTrackH + 24;
  const H = Math.max(120, aboveH + GENOME_H + belowH);"""
    )
    html = html.replace(
        """        const laneEnds = [];   // laneEnds[lane] = pixel x where last island ends
        const LANE_PAD = 2;    // minimum pixel gap between islands in same lane
        const islLanes = catIslands.map(isl => {
          const x = bpToX(isl.start);
          const w = Math.max(4, bpWidth(isl.end - isl.start));
          // Find lowest lane with room
          let lane = 0;
          while (laneEnds[lane] !== undefined && laneEnds[lane] > x - LANE_PAD) {
            lane++;
          }
          laneEnds[lane] = x + w;
          return lane;
        });

        const maxLane = Math.max(...islLanes);""",
        """        const laneEnds = [];   // laneEnds[lane] = pixel x where last island ends
        const LANE_PAD = 2;    // minimum pixel gap between islands in same lane
        const overflowLane = Math.max(0, GI_LANE_LIMIT - 1);
        let overflowCount = 0;
        const islLanes = catIslands.map(isl => {
          const x = bpToX(isl.start);
          const w = Math.max(4, bpWidth(isl.end - isl.start));
          // Find lowest lane with room, but fold excess density into an overflow lane.
          let lane = 0;
          while (lane < overflowLane && laneEnds[lane] !== undefined && laneEnds[lane] > x - LANE_PAD) {
            lane++;
          }
          if (lane >= GI_LANE_LIMIT || (lane === overflowLane && laneEnds[lane] !== undefined && laneEnds[lane] > x - LANE_PAD)) {
            lane = overflowLane;
            overflowCount++;
          }
          laneEnds[lane] = Math.max(laneEnds[lane] || 0, x + w + LANE_PAD);
          return lane;
        });

        const maxLane = Math.min(overflowLane, Math.max(...islLanes));"""
    )
    html = html.replace(
        "          ctx.fillText(lane===0 ? 'GI' : 'GI+'+lane, PAD-2, labelY);",
        "          ctx.fillText(lane===0 ? 'GI' : (overflowCount && lane===maxLane ? 'GI+'+lane+'+' : 'GI+'+lane), PAD-2, labelY);"
    )
    html = html.replace(
        '      <div class="minimap-wrap"><canvas id="${mid}" class="minimap" style="height:28px"></canvas></div>\n'
        '      <div class="track-wrap"><canvas id="${cid}"></canvas></div>',
        '      <div class="minimap-wrap"><canvas id="${mid}" class="minimap" style="height:28px"></canvas></div>\n'
        '      <div class="track-wrap"><canvas id="${cid}"></canvas></div>\n'
        '      ${renderIslandTable(strain)}'
    )

    # Avoid optional chaining and Array.flatMap in generated standalone files so
    # the viewer tolerates older Safari builds a little better.
    html = html.replace(
        "Math.max(...isolates.map(s => STRAIN_DATA[s]?.length || 5e6))",
        "Math.max(...isolates.map(s => (STRAIN_DATA[s] && STRAIN_DATA[s].length) || 5e6))"
    )
    html = html.replace(
        "Math.max(...isolates.map(s=>STRAIN_DATA[s]?.length||5e6))",
        "Math.max(...isolates.map(s=>(STRAIN_DATA[s]&&STRAIN_DATA[s].length)||5e6))"
    )
    html = html.replace(
        "new Set(isolates.flatMap(s=>(STRAIN_DATA[s]?.df||[]).map(d=>d.subtype))).size",
        "new Set([].concat(...isolates.map(s=>((STRAIN_DATA[s]&&STRAIN_DATA[s].df)||[]).map(d=>d.subtype)))).size"
    )
    html = html.replace(
        "(STRAIN_DATA[s]?.prophage||[]).length",
        "((STRAIN_DATA[s]&&STRAIN_DATA[s].prophage)||[]).length"
    )
    html = html.replace(
        "[...(STRAIN_DATA[d.strain]?.df||[]),...(STRAIN_DATA[d.strain]?.prophage||[])]",
        "[...((STRAIN_DATA[d.strain]&&STRAIN_DATA[d.strain].df)||[]),...((STRAIN_DATA[d.strain]&&STRAIN_DATA[d.strain].prophage)||[])]"
    )

    return html


def generate_patient_html(patient_id, strains, template_html, data_dir):
    """
    Generate a self-contained HTML for one patient.
    Embeds all strain data inline — no external fetch needed.
    """
    # Load data for all isolates of this patient
    blast_rows = load_fixed_blast_annotations(data_dir)
    strain_data = {}
    for strain in strains:
        strain_data[strain] = add_blast_annotations(
            strain, load_strain_data(strain, data_dir), blast_rows)

    # Build minimal patient groups (just this patient)
    patient_groups = {patient_id: strains}

    # Serialize
    strain_data_json  = json.dumps(strain_data)
    patient_json      = json.dumps(patient_groups)

    # Extract JS/CSS from template (everything between <script> and </script>
    # and <style> and </style> tags for the main viewer logic)
    # Replace the lazy-loading STRAIN_DATA with embedded data
    html = template_html

    # 1. Replace empty STRAIN_DATA with embedded data
    sd_replacement = f'const STRAIN_DATA = {strain_data_json};'
    html = re.sub(
        r'const STRAIN_DATA\s*=\s*\{\};.*?// populated on demand via fetch',
        lambda m: sd_replacement,
        html, flags=re.DOTALL
    )

    # 2. Replace PATIENT_GROUPS with single-patient groups
    pg_replacement = f'const PATIENT_GROUPS = {patient_json};'
    html = re.sub(
        r'const PATIENT_GROUPS\s*=\s*\{.*?\};',
        lambda m: pg_replacement,
        html, flags=re.DOTALL
    )

    # 3. Remove DATA_DIR and _fetchCache lines (no longer needed)
    html = re.sub(r"^\s*const DATA_DIR\s*=.*?;\s*(?://.*)?\n", "", html,
                  flags=re.MULTILINE)
    html = re.sub(r"^\s*let\s+_fetchCache\s*=.*?;\s*(?://.*)?\n", "", html,
                  flags=re.MULTILINE)

    # 4. Replace async fetchStrainData with sync version that just returns
    #    from the already-embedded STRAIN_DATA
    old_fetch = re.search(
        r'// Fetch strain data from per-strain JSON.*?^async function redrawAll',
        html, re.DOTALL | re.MULTILINE
    )
    if old_fetch:
        html = html[:old_fetch.start()] + \
               """// Data is embedded inline — no fetch needed
function fetchStrainData(strain) {
  return Promise.resolve(STRAIN_DATA[strain] || {});
}

async function redrawAll""" + \
               html[old_fetch.end() - len("async function redrawAll")
                    + len("async function redrawAll"):]

    # 5. Update the title to show patient ID
    html = re.sub(
        r'<title>.*?</title>',
        f'<title>Patient {patient_id} — M. abscessus Genomic Islands</title>',
        html
    )

    # 6. Add a note in the header
    html = re.sub(
        r'(PATIENT LONGITUDINAL VIEWER.*?serial isolate defense &amp; prophage comparison)',
        f'\\1 — Patient {patient_id}',
        html
    )

    # 7. Pre-select this patient in the dropdown
    html = re.sub(
        r"document\.addEventListener\('DOMContentLoaded'.*?}\);",
        f"""document.addEventListener('DOMContentLoaded', function() {{
  const sel = document.getElementById('patientSelect');
  if (sel) {{
    sel.value = '{patient_id}';
    sel.dispatchEvent(new Event('change'));
  }}
}});""",
        html, flags=re.DOTALL
    )

    html = make_generated_html_safari_safe(html)

    return html


def main():
    args = parse_args()

    # Find template
    template_path = args.template or os.path.join(
        args.html_dir, "patient_comparison_viewer.html")
    if not os.path.exists(template_path):
        print(f"ERROR: Template not found at {template_path}")
        print("Use --template to specify the path explicitly")
        return

    print(f"Loading template: {template_path}")
    with open(template_path) as f:
        template_html = f.read()

    # Load patient groups
    patient_groups = load_patient_groups(template_html)
    print(f"Found {len(patient_groups)} patients in template")

    if args.all:
        # Generate all patients
        os.makedirs(args.out_dir, exist_ok=True)
        n_generated = 0
        for pid, strains in sorted(patient_groups.items()):
            out_path = os.path.join(args.out_dir, f"{pid}_viewer.html")
            print(f"  Generating {pid} ({len(strains)} isolates)...", end=" ")
            try:
                html = generate_patient_html(
                    pid, strains, template_html, args.data_dir)
                with open(out_path, "w") as f:
                    f.write(html)
                size_kb = os.path.getsize(out_path) / 1024
                print(f"{size_kb:.0f}KB")
                n_generated += 1
            except Exception as e:
                print(f"ERROR: {e}")
        print(f"\nGenerated {n_generated} patient HTML files in {args.out_dir}/")

    elif args.patient:
        # Single patient
        pid = args.patient
        if pid not in patient_groups:
            # Try prefix match
            matches = [k for k in patient_groups if k.startswith(pid)]
            if len(matches) == 1:
                pid = matches[0]
                print(f"Matched patient: {pid}")
            elif len(matches) > 1:
                print(f"Multiple matches for '{args.patient}': {matches}")
                print("Please be more specific")
                return
            else:
                print(f"Patient '{args.patient}' not found")
                print(f"Available: {sorted(patient_groups.keys())[:10]}...")
                return

        strains   = patient_groups[pid]
        out_path  = args.out or f"{pid}_viewer.html"

        print(f"Generating {pid} ({len(strains)} isolates)...")
        html = generate_patient_html(pid, strains, template_html, args.data_dir)
        with open(out_path, "w") as f:
            f.write(html)

        size_kb = os.path.getsize(out_path) / 1024
        print(f"Written: {out_path} ({size_kb:.0f}KB)")
        print(f"Open with: open {out_path}")

    else:
        print("Specify --patient GD233 or --all")
        print(f"Available patients: {sorted(patient_groups.keys())[:10]}...")


if __name__ == "__main__":
    main()
