"""Supplementary Fig. 5 (tag helper traits vs HbA1c) redrawn in the house style.

Values are the dot positions recovered from the published figure (pixel
precision ~0.0003 in correlation units); labels are transcribed from it.
"""
import sys, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2] / 'src'))
from correlated_traits.style import BLUE, ORANGE, plt, save
S = sys.argv[1]
d = json.load(open(S + '/extract.json'))
grid = sorted(x for x in d['gridlines'] if 1400 < x < 3200)
assert len(grid) == 8
a, b = np.polyfit(grid, np.arange(-0.2, 0.51, 0.1), 1)
values = np.array([a * x + b for _, x in d['centroids']])      # top to bottom
labels = """helper_total_sugar
helper_cardiomyopathy
nmr_Concentration of Small HDL Particles
nmr_Apolipoprotein B to Apolipoprotein A1 ratio
nmr_Acetate
nmr_Acetone
helper_num_moderate_workout
helper_trans_fat
nmr_Triglycerides to Total Lipids in Large VLDL percentage
nmr_Histidine
helper_cholesterol
helper_smoke
nmr_Cholesteryl Esters to Total Lipids in Medium LDL percentage
nmr_Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio
nmr_Phospholipids to Total Lipids in Large LDL percentage
helper_pulmonary
helper_num_vigorous_workout
nmr_Lactate
nmr_Citrate
nmr_Phospholipids to Total Lipids in IDL percentage
nmr_Phenylalanine
nmr_Concentration of Medium HDL Particles
nmr_Phospholipids to Total Lipids in Chylomicrons and Extremely Large VLDL percentage
helper_chronic_rheumatic_heart_disease
nmr_Cholesterol to Total Lipids in Chylomicrons and Extremely Large VLDL percentage
nmr_Creatinine
nmr_Glutamine
helper_cerebrovascular
helper_arteries
nmr_Glycine
nmr_Cholesteryl Esters to Total Lipids in Very Large HDL percentage
helper_systolic_blood_pressure
nmr_Albumin
nmr_Cholesteryl Esters in IDL
nmr_Degree of Unsaturation
nmr_Average Diameter for LDL Particles
nmr_Free Cholesterol to Total Lipids in Very Large VLDL percentage
nmr_Triglycerides in Very Small VLDL
nmr_Valine
nmr_Cholesteryl Esters in Very Large HDL
nmr_Cholesterol to Total Lipids in Small LDL percentage
helper_ischemic
nmr_Spectrometer-corrected alanine
nmr_Linoleic Acid to Total Fatty Acids percentage
nmr_Glycoprotein Acetyls
helper_hypertension
helper_waist_circum
nmr_Cholesteryl Esters to Total Lipids in Very Small VLDL percentage
nmr_Glucose-lactate""".split('\n')
assert len(labels) == len(values) == 49
order = np.argsort(np.abs(values))
print('|r| ascending top->bottom?', bool((order == np.arange(49)).all()), 'first', values[:3].round(3), 'last', values[-3:].round(3))
y = np.arange(49)[::-1]
fig, ax = plt.subplots(figsize=(10, 11))
ax.hlines(y, 0, values, color=BLUE, lw=2, zorder=2)
ax.plot(values, y, 'o', color=ORANGE, ms=7, mew=0, zorder=3)
ax.axvline(0, color='0.3', lw=0.8, zorder=1)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=10)
ax.set_xticks(np.arange(-0.2, 0.51, 0.1)); ax.set_xlim(-0.26, 0.52); ax.set_ylim(-1, 49)
ax.grid(axis='x', color='0.85', ls='--', lw=0.8); ax.set_axisbelow(True)
ax.set_xlabel('Correlation with HbA1c'); ax.set_ylabel('Helper features')
fig.tight_layout()
print(save(fig, S + '/FigureS5_tag_helpers'))
json.dump(dict(zip(labels, values.round(4).tolist())), open(S + '/values.json', 'w'), indent=1)
