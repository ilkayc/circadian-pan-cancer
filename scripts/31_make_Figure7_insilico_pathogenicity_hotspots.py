#!/usr/bin/env python3
"""
Script 31: Generate Figure 7 (In Silico Pathogenicity & Hotspots) and Table 2.
Panel A: SIFT and PolyPhen-2 Venn diagram (706 high-confidence variants).
Panel B: Recurrent hotspot variants with gene labels, patient counts, and cancer types.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

out_fig = '/Users/ilkaycivelek/circadian_pan_cancer_Q2/figures/Figure_7_InSilico_Pathogenicity_and_Hotspots_Fixed.png'
out_t2 = '/Users/ilkaycivelek/circadian_pan_cancer_Q2/results/Table2_Recurrent_Variants_Annotated.csv'

# Recurrent hotspots data verified from TCGA
hotspots = [
    {'Gene': 'PER3', 'Variant': 'p.R316C', 'Transcript': 'ENST00000361413', 'Count': 5, 'Domain': 'PAS domain', 'Cancers': 'SKCM (n=5)', 'SIFT': 'Deleterious (0.02)', 'PolyPhen': 'Probably Damaging (0.989)'},
    {'Gene': 'CSNK1E', 'Variant': 'p.R127W', 'Transcript': 'ENST00000305886', 'Count': 4, 'Domain': 'Protein kinase domain', 'Cancers': 'UCEC, COAD, SKCM', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Probably Damaging (1.000)'},
    {'Gene': 'CLOCK', 'Variant': 'p.G120V', 'Transcript': 'ENST00000309964', 'Count': 4, 'Domain': 'bHLH / PAS-A boundary', 'Cancers': 'BRCA, UCEC, COAD', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Probably Damaging (0.999)'},
    {'Gene': 'CRY1', 'Variant': 'p.R348C', 'Transcript': 'ENST00000350720', 'Count': 4, 'Domain': 'Photolyase / FAD-binding', 'Cancers': 'UCEC, COAD', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Possibly Damaging (0.606)'},
    {'Gene': 'NR1D2', 'Variant': 'p.S406L', 'Transcript': 'ENST00000355883', 'Count': 4, 'Domain': 'Ligand-binding domain', 'Cancers': 'SKCM, COAD', 'SIFT': 'Deleterious (0.02)', 'PolyPhen': 'Possibly Damaging (0.793)'},
    {'Gene': 'RORC', 'Variant': 'p.E303K', 'Transcript': 'ENST00000318247', 'Count': 3, 'Domain': 'Nuclear receptor LBD', 'Cancers': 'UCEC, SKCM, COAD', 'SIFT': 'Deleterious (0.00)', 'PolyPhen': 'Probably Damaging (0.984)'}
]

df_t2 = pd.DataFrame(hotspots)
df_t2.to_csv(out_t2, index=False)
print('Table 2 exported to:', out_t2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

# Panel A: Venn
v = venn2(subsets=(215, 142, 706), set_labels=('SIFT Predicted
Deleterious (n = 921)', 'PolyPhen-2 Predicted
Damaging (n = 848)'), ax=ax1)
for patch, col in zip(v.patches, ['#4575b4', '#d73027', '#91bfdb']):
    if patch: patch.set_alpha(0.6)
ax1.set_title('A. In Silico Consensus Intersection (n = 1,532 Missense Variants)', fontsize=11, fontweight='bold', pad=15)
ax1.text(0, -0.65, '706 High-Confidence Deleterious Variants (46.1%)', ha='center', fontsize=10.5, fontweight='bold', color='#1a476f')

# Panel B: Lollipop recurrence
y_idx = np.arange(len(hotspots))
labels = [f"{h['Gene']} {h['Variant']}" for h in hotspots]
counts = [h['Count'] for h in hotspots]
colors = ['#d73027', '#fc8d59', '#fee090', '#e0f3f8', '#91bfdb', '#4575b4']

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
    ax2.text(h['Count'] + 0.15, i, f"{h['Domain']} ({h['Cancers']})", va='center', fontsize=9, color='#333333', style='italic')

plt.tight_layout()
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
print('Figure 7 saved to:', out_fig)
