import os
import csv
import re
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from collections import Counter, defaultdict
from scipy.stats import kruskal


BG    = '#fdf4ff'
DARK  = '#4a044e'
MID   = '#9d4edd'
GRID  = '#e9d5ff'
SPINE = '#c4b5fd'

DCC_COLORS = {
    'DCC1':'#6b21a8','DCC2':'#7c3aed','DCC3':'#9333ea','DCC4':'#a855f7',
    'DCC5':'#c084fc','DCC6':'#d946ef','DCC7':'#ec4899','Non-DCC':'#9ca3af',
}
CARGO_COLORS = {
    'regulatory':'#6b21a8','efflux':'#ec4899','metal':'#f59e0b',
    'mobility':  '#0ea5e9','defense':'#10b981','hypothetical':'#94a3b8',
    'unknown':   '#e2e8f0','ta_system':'#fde68a',
}
DCC_LIST   = ['DCC1','DCC2','DCC3','DCC4','DCC5','DCC7','Non-DCC']
CARGO_LIST = ['regulatory','efflux','metal','mobility','defense','hypothetical']


def b(v): return str(v).strip().upper() == 'TRUE'
def i(v):
    try: return int(v)
    except: return 0

def base_strain(s):
    return re.sub(r'[_\-](hybrid|WGS|wgs|short).*', '', s,
                  flags=re.IGNORECASE).strip()

def style(ax):
    ax.set_facecolor(BG)
    ax.yaxis.grid(True, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    ax.spines['left'].set_color(SPINE)
    ax.spines['bottom'].set_color(SPINE)
    ax.tick_params(colors=MID, labelsize=9)


def load_data(islands_tsv, master_csv):
    with open(islands_tsv, newline='', encoding='utf-8-sig') as f:
        islands = list(csv.DictReader(f, delimiter='\t'))
    with open(master_csv, newline='', encoding='utf-8-sig') as f:
        master = {r['strain']: r for r in csv.DictReader(f)}

    for isl in islands:
        strain = base_strain(isl['strain'])
        m = master.get(strain, {})
        isl['_strain'] = strain
        isl['_dcc']    = m.get('dcc', '')
        isl['_subsp']  = m.get('subspecies', '')
        isl['_seq']    = m.get('sequenced', '')

    all_seq = [r for r in islands if b(r['_seq'])]
    # KEY FILTER: only islands containing at least one defense gene
    def_isl = [r for r in all_seq if i(r['n_defense']) > 0]
    sequenced_strains = {s: r for s, r in master.items()
                         if b(r.get('sequenced', ''))}

    return all_seq, def_isl, sequenced_strains


def generate_island_landscape(islands_tsv, master_csv, output_path, dpi=200):
    print('Loading data...')
    all_seq, seq_isl, sequenced_strains = load_data(islands_tsv, master_csv)

    n_total          = len(all_seq)
    n_shown         = len(seq_isl)
    n_strains        = len(set(r['_strain'] for r in seq_isl))
    total_def_genes  = sum(i(r['n_defense']) for r in seq_isl)

    print(f'  All islands (sequenced):        {n_total}')
    print(f'  Defense-containing islands:     {n_def} '
          f'({round(n_def/n_total*100,1)}%)')
    print(f'  Strains with >=1 defense island:{n_strains}')
    print(f'  Total defense gene instances:   {total_def_genes}')

    # Per-strain island counts by DCC
    strain_islands    = defaultdict(list)
    for isl in seq_isl:
        strain_islands[isl['_strain']].append(isl)

    dcc_island_counts = defaultdict(list)
    for strain, row in sequenced_strains.items():
        dcc = row.get('dcc', '')
        if dcc not in DCC_LIST: continue
        dcc_island_counts[dcc].append(len(strain_islands.get(strain, [])))

    dcc_plot  = [d for d in DCC_LIST if dcc_island_counts.get(d)]
    vals_list = [dcc_island_counts[d] for d in dcc_plot]

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 13), facecolor=BG)
    fig.patch.set_facecolor(BG)
    gs  = GridSpec(3, 3, figure=fig,
                   left=0.07, right=0.97, top=0.92, bottom=0.06,
                   hspace=0.45, wspace=0.35,
                   height_ratios=[1.3, 1.2, 1.2])

    fig.text(0.5, 0.975,
             'Defense-Carrying Genomic Island Landscape — '
             'M. abscessus Clinical Cohort',
             ha='center', va='top', fontsize=14,
             fontweight='bold', color=DARK)
    fig.text(0.5, 0.955,
             (f'{n_shown} defense-containing islands '
              f'({round(n_shown/n_total*100,1)}% of all {n_total} islands)  ·  '
              f'{n_strains} strains  ·  {total_def_genes} defense gene instances'
              if defense_only else
              f'{n_shown} total genomic islands  ·  '
              f'{n_strains} strains  ·  {total_def_genes} defense gene instances'),
             ha='center', va='top', fontsize=9, color=MID)

    # ── Panel A: violin by DCC ────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.set_facecolor(BG)
    bp = ax1.violinplot(vals_list, positions=range(len(dcc_plot)),
                         showmedians=True, showextrema=False)
    for pc, dcc in zip(bp['bodies'], dcc_plot):
        pc.set_facecolor(DCC_COLORS[dcc]); pc.set_alpha(0.7)
    bp['cmedians'].set_color(DARK); bp['cmedians'].set_linewidth(2)

    np.random.seed(42)
    for xi, (dcc, vals) in enumerate(zip(dcc_plot, vals_list)):
        jitter = np.random.uniform(-0.12, 0.12, len(vals))
        ax1.scatter(np.full(len(vals), xi)+jitter, vals,
                    color=DCC_COLORS[dcc], alpha=0.3, s=14, zorder=3)
        ax1.text(xi, max(vals)+0.2,
                 f'μ={round(np.mean(vals),1)}\nn={len(vals)}',
                 ha='center', va='bottom', fontsize=8.5,
                 color=DARK, fontweight='bold', linespacing=1.3)

    stat, p = kruskal(*vals_list)
    p_label = 'p < 0.0001' if p < 0.0001 else f'p = {round(p,4)}'
    ax1.set_xticks(range(len(dcc_plot)))
    ax1.set_xticklabels(dcc_plot, fontsize=11, color=DARK, fontweight='bold')
    ax1.set_ylabel('Defense-Carrying Islands per Strain', fontsize=10,
                   color=DARK, fontweight='bold')
    ax1.set_title(f'Defense Island Burden by DCC Lineage  '
                  f'(Kruskal-Wallis {p_label})',
                  fontsize=11, color=DARK, fontweight='bold', pad=8)
    style(ax1)
    ax1.tick_params(colors=MID, labelsize=10)

    # ── Panel B: defense genes per island ─────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.set_facecolor(BG)
    def_per_isl = Counter(i(r['n_defense']) for r in seq_isl)
    max_def     = max(def_per_isl.keys())
    counts      = [def_per_isl.get(k, 0) for k in range(1, max_def+1)]
    bars2 = ax2.bar(range(1, max_def+1), counts,
                    color=MID, linewidth=0, zorder=3, width=0.75, alpha=0.85)
    for bar, val in zip(bars2, counts):
        if val > 0:
            ax2.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
                     str(val), ha='center', va='bottom', fontsize=8,
                     color=DARK, fontweight='bold')
    ax2.set_xlabel('Defense Genes per Island', fontsize=9.5,
                   color=DARK, fontweight='bold')
    ax2.set_ylabel('Number of Islands', fontsize=9.5,
                   color=DARK, fontweight='bold')
    ax2.set_title('Defense Genes per\nDefense Island',
                  fontsize=11, color=DARK, fontweight='bold', pad=8)
    ax2.set_xticks(range(1, max_def+1))
    ax2.set_ylim(0, max(counts)*1.2 if not defense_only else 700)
    style(ax2)
    mean_def = total_def_genes / n_shown if n_shown else 0
    ax2.text(0.97, 0.97,
             f'Mean: {round(mean_def,2)}\nMax: {max_def}\nTotal: {total_def_genes}',
             transform=ax2.transAxes, ha='right', va='top', fontsize=9,
             color=DARK,
             bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                       edgecolor=SPINE, alpha=0.9))

    # ── Panel C: cargo by DCC ─────────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, :2])
    ax3.set_facecolor(BG)
    x       = np.arange(len(dcc_plot))
    bottoms = np.zeros(len(dcc_plot))
    for cargo in CARGO_LIST:
        vals = []
        for dcc in dcc_plot:
            dcc_isl = [r for r in seq_isl if r['_dcc']==dcc]
            cargo_n = sum(1 for r in dcc_isl if r['dominant_cargo']==cargo)
            vals.append(cargo_n/len(dcc_isl)*100 if dcc_isl else 0)
        ax3.bar(x, vals, 0.65, bottom=bottoms,
                color=CARGO_COLORS.get(cargo, '#e2e8f0'),
                label=cargo, linewidth=0, zorder=3, alpha=0.9)
        bottoms += np.array(vals)
    ax3.set_xticks(x)
    ax3.set_xticklabels(dcc_plot, fontsize=11, color=DARK, fontweight='bold')
    ax3.set_ylabel('% of Defense Islands', fontsize=10,
                   color=DARK, fontweight='bold')
    ax3.set_title('Defense Island Cargo Composition by DCC\n'
                  '(dominant cargo — most defense islands are within '
                  'regulatory islands)',
                  fontsize=11, color=DARK, fontweight='bold', pad=8)
    ax3.set_ylim(0, 108)
    style(ax3)
    ax3.tick_params(colors=MID, labelsize=10)
    leg3 = ax3.legend(fontsize=9, frameon=True, edgecolor=GRID,
                      facecolor='white', loc='upper right', ncol=3)
    for t in leg3.get_texts(): t.set_color(DARK)

    # ── Panel D: size distribution ────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.set_facecolor(BG)
    lengths = np.array([i(r['length'])/1000 for r in seq_isl
                        if i(r['length']) > 0])
    ax4.hist(lengths, bins=30, color=MID, alpha=0.85, linewidth=0, zorder=3)
    ax4.axvline(x=np.mean(lengths), color=DARK, linewidth=1.8,
                linestyle='--',
                label=f'Mean: {round(np.mean(lengths),1)} kb')
    ax4.axvline(x=np.median(lengths), color='#ec4899', linewidth=1.8,
                linestyle=':',
                label=f'Median: {round(np.median(lengths),1)} kb')
    ax4.set_xlabel('Island Length (kb)', fontsize=10,
                   color=DARK, fontweight='bold')
    ax4.set_ylabel('Number of Islands', fontsize=10,
                   color=DARK, fontweight='bold')
    ax4.set_title('Defense Island Size\nDistribution',
                  fontsize=11, color=DARK, fontweight='bold', pad=8)
    style(ax4)
    leg4 = ax4.legend(fontsize=8.5, frameon=True, edgecolor=GRID,
                      facecolor='white')
    for t in leg4.get_texts(): t.set_color(DARK)

    # ── Panel E: top catalog groups ───────────────────────────────────────────
    ax5 = fig.add_subplot(gs[2, :2])
    ax5.set_facecolor(BG)
    group_counts = Counter(r['group_id'] for r in seq_isl)
    top_groups   = [g for g,_ in group_counts.most_common(15)]
    top_counts   = [group_counts[g] for g in top_groups]
    top_cargo, top_avg_def = [], []
    for g in top_groups:
        g_isl = [r for r in seq_isl if r['group_id']==g]
        c = Counter(r['dominant_cargo'] for r in g_isl)
        top_cargo.append(c.most_common(1)[0][0] if c else 'unknown')
        top_avg_def.append(
            round(np.mean([i(r['n_defense']) for r in g_isl]), 1))

    bar_colors = [CARGO_COLORS.get(c, '#94a3b8') for c in top_cargo]
    bars5 = ax5.barh(range(len(top_groups)), top_counts,
                      color=bar_colors, linewidth=0, height=0.65,
                      zorder=3, alpha=0.9)
    for bar, val, avg_d in zip(bars5, top_counts, top_avg_def):
        pct = round(val/n_strains*100, 0)
        ax5.text(bar.get_width()+1,
                 bar.get_y()+bar.get_height()/2,
                 f'{val} ({pct:.0f}% of strains)  avg {avg_d} def genes',
                 va='center', fontsize=8.5, color=DARK)
    ax5.set_yticks(range(len(top_groups)))
    ax5.set_yticklabels(top_groups, fontsize=9, color=DARK,
                         fontweight='bold', family='monospace')
    ax5.set_xlabel('Number of Occurrences Across Cohort', fontsize=10,
                   color=DARK, fontweight='bold')
    ax5.set_title('Top 15 Most Widespread Defense-Carrying Catalog Islands',
                  fontsize=11, color=DARK, fontweight='bold', pad=8)
    ax5.set_xlim(0, max(top_counts)*1.45)
    ax5.xaxis.grid(True, color=GRID, linewidth=0.7, zorder=0)
    ax5.set_axisbelow(True)
    for sp in ['top','right']: ax5.spines[sp].set_visible(False)
    ax5.spines['left'].set_color(SPINE)
    ax5.spines['bottom'].set_color(SPINE)
    ax5.tick_params(colors=MID, labelsize=9)
    leg_patches = [mpatches.Patch(color=CARGO_COLORS[c], label=c)
                   for c in CARGO_LIST]
    leg5 = ax5.legend(handles=leg_patches, fontsize=8.5, frameon=True,
                      edgecolor=GRID, facecolor='white', loc='lower right',
                      title='Dominant cargo', title_fontsize=9)
    leg5.get_title().set_color(DARK)

    # ── Panel F: confidence + tRNA ────────────────────────────────────────────
    ax6 = fig.add_subplot(gs[2, 2])
    ax6.set_facecolor(BG)
    dcc_high_pct, dcc_trna_pct = [], []
    for dcc in dcc_plot:
        dcc_isl = [r for r in seq_isl if r['_dcc']==dcc]
        if not dcc_isl: continue
        dcc_high_pct.append(
            sum(1 for r in dcc_isl if r['confidence']=='high')
            / len(dcc_isl) * 100)
        dcc_trna_pct.append(
            sum(1 for r in dcc_isl if r['trna_flanked']=='Yes')
            / len(dcc_isl) * 100)
    x6, w6 = np.arange(len(dcc_plot)), 0.35
    ax6.bar(x6-w6/2, dcc_high_pct, w6, label='High confidence',
            color='#6b21a8', linewidth=0, zorder=3, alpha=0.9)
    ax6.bar(x6+w6/2, dcc_trna_pct, w6, label='tRNA-flanked',
            color='#0ea5e9', linewidth=0, zorder=3, alpha=0.9)
    for xi, (h, t) in enumerate(zip(dcc_high_pct, dcc_trna_pct)):
        ax6.text(xi-w6/2, h+1, f'{h:.0f}%', ha='center', va='bottom',
                 fontsize=8, color='#6b21a8', fontweight='bold')
        ax6.text(xi+w6/2, t+1, f'{t:.0f}%', ha='center', va='bottom',
                 fontsize=8, color='#0ea5e9', fontweight='bold')
    ax6.set_xticks(x6)
    ax6.set_xticklabels(dcc_plot, fontsize=8.5, color=DARK,
                         fontweight='bold', rotation=30, ha='right')
    ax6.set_ylabel('% of Defense Islands', fontsize=9.5,
                   color=DARK, fontweight='bold')
    ax6.set_title('Island Quality by DCC\n(high confidence + tRNA-flanking)',
                  fontsize=11, color=DARK, fontweight='bold', pad=8)
    ax6.set_ylim(0, 115)
    style(ax6)
    leg6 = ax6.legend(fontsize=9, frameon=True, edgecolor=GRID,
                      facecolor='white')
    for t in leg6.get_texts(): t.set_color(DARK)

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor=BG, edgecolor='none')
    plt.close()
    print(f'Saved: {output_path}')


def main():
    parser = argparse.ArgumentParser(
        description='Generate defense-carrying genomic island landscape figure.'
    )
    parser.add_argument('--islands_tsv', default='catalog_strains.tsv',
                        help='MAbsIslandScanner catalog TSV')
    parser.add_argument('--master_csv',  default='mabs_cohort_master.csv')
    parser.add_argument('--output',      default='genomic_island_landscape.png')
    parser.add_argument('--dpi',         type=int, default=200)
    parser.add_argument('--defense_only', action='store_true', default=True,
                        help='Only show islands with n_defense >= 1 (default: True)')
    parser.add_argument('--all_islands',  action='store_true',
                        help='Show all genomic islands regardless of defense content')
    args = parser.parse_args()
    if args.all_islands:
        args.defense_only = False

    generate_island_landscape(args.islands_tsv, args.master_csv,
                               args.output,
                               defense_only=args.defense_only,
                               dpi=args.dpi)


if __name__ == '__main__':
    main()
