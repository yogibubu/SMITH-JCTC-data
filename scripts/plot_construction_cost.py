#!/usr/bin/env python3
"""Plot the measurements reported in SI Table S16."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures' / 'construction_cost.pdf'
# System, family, atoms, candidates, build time (s), one B evaluation (s),
# peak Python allocation (MiB). Values are the SI Table S16 medians.
ROWS = [
    ('Benzene', 'PAH', 12, 30, 0.131, 0.0054, 0.23),
    ('Naphthalene', 'PAH', 18, 48, 0.139, 0.0093, 0.28),
    ('Anthracene', 'PAH', 24, 66, 0.243, 0.0139, 0.51),
    ('Phenanthrene', 'PAH', 24, 66, 0.211, 0.0131, 0.52),
    ('Pyrene', 'PAH', 26, 78, 0.256, 0.0162, 0.59),
    ('Coronene', 'PAH', 36, 120, 0.722, 0.0343, 1.08),
    ('Cubane', 'Ring/cage', 16, 54, 0.400, 0.0139, 0.45),
    ('Norbornane', 'Ring/cage', 19, 52, 0.099, 0.0120, 0.29),
    ('Cyclooctane', 'Ring/cage', 24, 66, 0.142, 0.0201, 0.53),
    ('Camphor', 'Ring/cage', 27, 76, 0.145, 0.0151, 0.44),
    ('Spiro', 'Ring/cage', 31, 87, 0.253, 0.0225, 0.75),
    ('Ferrocene D5h', 'Special', 21, 651, 3.324, 0.0642, 2.83),
    ('Ferrocene D5d', 'Special', 21, 656, 3.182, 0.0653, 2.59),
]
COLORS = {'PAH': '#315f88', 'Ring/cage': '#668d65', 'Special': '#ad6b35'}
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8,
    'axes.linewidth': 0.75, 'pdf.fonttype': 42,
})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.75), constrained_layout=True)
for family in COLORS:
    rows = [r for r in ROWS if r[1] == family]
    xs = [r[2] for r in rows]
    ax1.scatter(xs, [r[4] for r in rows], color=COLORS[family], marker='o', s=31, zorder=3)
    ax1.scatter(xs, [r[5] for r in rows], facecolors='none', edgecolors=COLORS[family], marker='s', s=30, linewidths=1.1, zorder=3)
    ax2.scatter(xs, [r[6] for r in rows], color=COLORS[family], marker='o', s=31, zorder=3)
for ax in (ax1, ax2):
    ax.set_xlim(9, 39)
    ax.set_xticks([12, 18, 24, 30, 36])
    ax.set_xlabel('Number of atoms')
    ax.grid(axis='y', color='#e4e7e8', linewidth=0.6)
    ax.spines[['top', 'right']].set_visible(False)
ax1.set_yscale('log')
ax1.set_ylim(0.003, 6)
ax1.set_ylabel('Elapsed time (s)')
ax2.set_ylim(0, 3.25)
ax2.set_ylabel('Peak Python allocation (MiB)')
ax1.text(0.02, 0.98, '(a)', transform=ax1.transAxes, va='top', fontweight='bold')
ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes, va='top', fontweight='bold')
ax1.annotate('ferrocene\n651–656 candidates', (21, 3.2), xytext=(23, 1.6), fontsize=7,
             arrowprops={'arrowstyle': '-', 'lw': 0.65, 'color': '#777777'})
legend = [
    Line2D([0], [0], marker='o', color='none', markerfacecolor='#555555', markeredgecolor='#555555', markersize=5, label='Construction'),
    Line2D([0], [0], marker='s', color='none', markerfacecolor='none', markeredgecolor='#555555', markersize=5, label='One B evaluation'),
]
ax1.legend(handles=legend, loc='lower right', fontsize=7, frameon=False)
legend2 = [Line2D([0], [0], marker='o', color='none', markerfacecolor=color,
                  markeredgecolor=color, markersize=5, label=family)
           for family, color in COLORS.items()]
ax2.legend(handles=legend2, loc='upper right', fontsize=7, frameon=False)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, bbox_inches='tight')
fig.savefig(OUT.with_suffix('.png'), dpi=240, bbox_inches='tight')
plt.close(fig)
print(OUT)
