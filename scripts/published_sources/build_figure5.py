"""Figure 5 (UK Biobank) redrawn in the house style from the archived SVG sources.

Panels B-D are replotted from the exact values stored in the SVGs; panel A is
the archived network drawing with its node colours mapped to the palette.
"""
import sys, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2] / 'src'))
from correlated_traits.style import BLUE, ORANGE, GREEN, DIVERGING, panel_letter, plt, save
from matplotlib.colors import Normalize
from matplotlib.image import imread
S = sys.argv[1]
bars = json.load(open(S + '/panelB.json'))
M = np.load(S + '/panelC.npy')
D = json.load(open(S + '/panelD.json'))
traits = ['num moderate workout', 'num vigorous workout', 'trans fat', 'total sugar', 'cholesterol',
          'smoke', 'waist circum', 'triglycerides', 'hdl', 'ldl', 'HbA1c', 'diastolic blood pressure',
          'systolic blood pressure', 'hypertension', 'cardiomyopathy', 'rheumatic heart disease',
          'pulmonary', 'cerebrovascular', 'arteries', 'ischemic']
assert M.shape == (20, 20)
print('check hdl-triglycerides', round(M[8, 7], 3), ' HbA1c-waist', round(M[10, 6], 3))

fig = plt.figure(figsize=(13, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], hspace=0.28, wspace=0.42,
                      left=0.05, right=0.97, top=0.96, bottom=0.07)

# A: network (recoloured archive drawing)
ax = fig.add_subplot(gs[0, 0])
net = imread(S + '/fig5/network_house.png')
ink = np.where((net[..., :3] < 0.98).any(2))
ax.imshow(net[ink[0].min():ink[0].max() + 1, ink[1].min():ink[1].max() + 1])
ax.set_axis_off()
cell = ax.get_position()
fig.text(cell.x0 - 0.03, cell.y1 + 0.005, 'A', fontsize=20, fontweight='bold', va='top', ha='left')

# B: pairwise-correlation counts per NMR trait
ax = fig.add_subplot(gs[0, 1])
ax.bar(np.arange(1, len(bars) + 1), bars, width=1.0, color=BLUE, linewidth=0)
ax.set_xlim(0, len(bars) + 1); ax.set_ylim(0, 62)
ax.set_xlabel('NMR trait index'); ax.set_ylabel(r'Number of pairwise correlations $> 0.8$')
panel_letter(ax, 'B')

# C: correlation heatmap of the 19 clinical and lifestyle traits plus HbA1c
ax = fig.add_subplot(gs[1, 0])
im = ax.imshow(M, cmap=DIVERGING, norm=Normalize(-1, 1), aspect='equal')
ax.set_xticks(range(20)); ax.set_xticklabels(traits, rotation=65, ha='right', fontsize=9)
ax.set_yticks(range(20)); ax.set_yticklabels(traits, fontsize=9)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.ax.set_ylim(M.min(), 1); cb.set_label('Pearson correlation')
panel_letter(ax, 'C', x=-0.42)

# D: AUC-ROC as tag helper traits are added
ax = fig.add_subplot(gs[1, 1])
n, auc = np.array(D['line2d_28']['x']).round().astype(int), np.array(D['line2d_28']['y'])
base, hba1c = D['line2d_29']['y'][0], D['hba1c']
ax.axhline(base, color=ORANGE, ls='--', lw=1.8, label='Baseline: T2D ~ Covariates + 37 SNPs', zorder=2)
ax.axhline(hba1c, color=GREEN, ls='--', lw=1.8, label='T2D ~ Covariates + HbA1c', zorder=2)
ax.plot(n, auc, '-o', color=BLUE, lw=1.8, ms=5, mew=0, label='T2D ~ Baseline + Tag Helper Traits', zorder=3)
for k in (1, 2, 3, 5, 10, 20, 30):
    ax.annotate('%.3f' % auc[k], (k, auc[k]), xytext=(6, -12) if k < 5 else (0, 8),
                textcoords='offset points', ha='left' if k < 5 else 'center', fontsize=10, color=BLUE)
ax.text(40, base - 0.012, '%.3f' % base, ha='right', va='top', fontsize=11, color=ORANGE, fontweight='bold')
ax.text(40, hba1c - 0.012, '%.3f' % hba1c, ha='right', va='top', fontsize=11, color=GREEN, fontweight='bold')
ax.set_xlim(0, 40); ax.set_ylim(0.6, 1.0)
ax.set_xlabel('Number of added tag helper traits'); ax.set_ylabel('AUC-ROC')
ax.grid(axis='y', color='0.88', ls='--', lw=0.8); ax.set_axisbelow(True)
ax.legend(loc='lower left', bbox_to_anchor=(0.04, 0.02))
panel_letter(ax, 'D')
print('annotations', {k: round(float(auc[k]), 3) for k in (1, 2, 3, 5, 10, 20, 30)}, 'base', round(base, 3), 'hba1c', round(hba1c, 3))
print(save(fig, S + '/fig5/Figure5_house'))
