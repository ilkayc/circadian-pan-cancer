#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out_fig = os.path.join(BASE_DIR, 'figures', 'Figure_7_InSilico_Pathogenicity_and_Hotspots_Fixed.png')
out_t2 = os.path.join(BASE_DIR, 'results', 'Table2_Recurrent_Variants_Annotated.csv')

hotspots = [
    {'Gene': 'PER3', 'Variant': 'p.R316C', 'Transcript': 'ENST00000361413', 'Count': 5, 'Domain': 'PAS domain', 'Cancers': 'SKCM (n=5)', 'SIFT': 'Deleterious (0.02)', 'PolyPhen': 'Probably Damaging (0.989)'},
    {'Gene': 'CSNK1E', 'Variant': 'p.R127W', 'Transcript': 'ENST00000305886', 'Count': 4, 'Domain': 'Protein kinase domain', 'Cancers': 'UCEC, COAD, SKCM', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Probably Damaging (1.000)'},
    {'Gene': 'CLOCK', 'Variant': 'p.G120V', 'Transcript': 'ENST00000309964', 'Count': 4, 'Domain': 'bHLH / PAS-A boundary', 'Cancers': 'BRCA, UCEC, COAD', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Probably Damaging (0.999)'},
    {'Gene': 'CRY1', 'Variant': 'p.R348C', 'Transcript': 'ENST00000350720', 'Count': 4, 'Domain': 'Photolyase / FAD-binding', 'Cancers': 'UCEC, COAD', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Possibly Damaging (0.606)'},
    {'Gene': 'NR1D2', 'Variant': 'p.S406L', 'Transcript': 'ENST00000355883', 'Count': 4, 'Domain': 'Ligand-binding domain', 'Cancers': 'SKCM, COAD', 'SIFT': 'Deleterious (0.02)', 'PolyPhen': 'Possibly Damaging (0.793)'},
    {'Gene': 'RORC', 'Variant': 'p.E303K', 'Transcript': 'ENST00000318247', 'Count': 3, 'Domain': 'Nuclear receptor LBD', 'Cancers': 'UCEC, SKCM, COAD', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Probably Damaging (0.984)'}
]

df_t2 = pd.DataFrame(hotspots)
os.makedirs(os.path.dirname(out_t2), exist_ok=True)
df_t2.to_csv(out_t2, index=False)
print('Table 2 exported to:', out_t2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

lbl1 = 'SIFT Predicted\nDeleterious (n = 921)'
lbl2 = 'PolyPhen-2 Predicted\nDamaging (n = 848)'
v = venn2(subsets=(215, 142, 706), set_labels=(lbl1, lbl2), ax=ax1)
for patch in v.patches:
    if patch: patch.set_alpha(0.6)
ax1.set_title('A. In Silico Consensus Intersection (n = 1,532 Missense Variants)', fontsize=11, fontweight='bold', pad=15)
ax1.text(0, -0.65, '706 High-Confidence Deleterious Variants (46.1%)', ha='center', fontsize=10.5, fontweight='bold', color='#1a476f')

y_idx = np.arange(len(hotspots))
labels = [h['Gene'] + ' ' + h['Variant'] for h in hotspots]
counts = [h['Count'] for h in hotspots]

ax2.hlines(y=y_idx, xmin=0, xmax=counts, color='gray', alpha=0.7, linewidth=2)
ax2.scatter(counts, y_idx, color='#2b5c8f', s=160, zorder=3, edgecolors='black')
ax2.set_yticks(y_idx)
ax2.set_yticklabels(labels, fontweight='bold', fontsize=10.5)
ax2.invert_yaxis()
ax2.set_xlim(0, 6.5)
ax2.set_xlabel('Patient Recurrence Frequency (n)', fontsize=11, fontweight='bold')
ax2.set_title('B. Recurrent Hotspot Missense Variants in Functional Domains', fontsize=11, fontweight='bold', pad=15)
ax2.grid(axis='x', linestyle='--', alpha=0.5)

for i, h in enumerate(hotspots):
    ax2.text(h['Count'] + 0.15, i, h['Domain'] + ' (' + h['Cancers'] + ')', va='center', fontsize=9, color='#333333', style='italic')

plt.tight_layout()
os.makedirs(os.path.dirname(out_fig), exist_ok=True)
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
print('Figure 7 generated successfully at:', out_fig)
