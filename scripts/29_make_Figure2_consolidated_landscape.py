#!/usr/bin/env python3
"""
Script 29: Generate Consolidated Figure 2 (Panels A & B)
Replaces overlapping Figures 2-4 per Reviewer 1 request.
Panel A: 2,503 total somatic mutations and cohort frequencies across 16 circadian genes.
Panel B: Functional consequence spectrum of 1,744 non-silent mutations.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

mut_file = '/Users/ilkaycivelek/results_circadian/circadian_all_mutations_extracted.csv'
out_fig = '/Users/ilkaycivelek/circadian_pan_cancer_Q2/figures/Figure_2_Consolidated_Mutational_Landscape.png'

df = pd.read_csv(mut_file, low_memory=False)
N_TOTAL_SAMPLES = 10239

gene_order = df['Hugo_Symbol'].value_counts().index.tolist()
gene_counts = df['Hugo_Symbol'].value_counts()[gene_order]
gene_pcts = (gene_counts / N_TOTAL_SAMPLES) * 100

# Non-silent categories
non_silent_df = df[df['Variant_Classification'] != 'Silent'].copy()
cats = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'In_Frame_Del', 'Other']
def map_cat(c):
    if c in cats: return c
    if 'Splice' in c: return 'Splice_Site'
    if 'Nonstop' in c or 'Start' in c: return 'Other'
    return 'Other'
non_silent_df['Class_Group'] = non_silent_df['Variant_Classification'].apply(map_cat)

ct = pd.crosstab(non_silent_df['Hugo_Symbol'], non_silent_df['Class_Group']).reindex(gene_order).fillna(0)
for c in cats:
    if c not in ct.columns: ct[c] = 0
ct = ct[cats]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5), dpi=300)

# Panel A
y_pos = np.arange(len(gene_order))
bars = ax1.barh(y_pos, gene_counts, color='#2b5c8f', edgecolor='none', height=0.7)
ax1.set_yticks(y_pos)
ax1.set_yticklabels(gene_order, fontweight='bold', fontsize=11)
ax1.invert_yaxis()
ax1.set_xlabel('Total Somatic Mutation Count (n)', fontsize=11, fontweight='bold')
ax1.set_title('A. Overall Somatic Mutation Burden Across 16 Clock Genes (n = 2,503)', fontsize=12, fontweight='bold', pad=12)
ax1.grid(axis='x', linestyle='--', alpha=0.5)

for bar, cnt, pct in zip(bars, gene_counts, gene_pcts):
    ax1.text(bar.get_width() + 8, bar.get_y() + bar.get_height()/2, f'{cnt} ({pct:.1f}%)', va='center', fontsize=9.5, color='#1a252f')
ax1.set_xlim(0, max(gene_counts) * 1.25)

# Panel B
palette = {'Missense_Mutation': '#2ca02c', 'Nonsense_Mutation': '#d62728', 'Frame_Shift_Del': '#1f77b4',
           'Frame_Shift_Ins': '#ff7f0e', 'Splice_Site': '#9467bd', 'In_Frame_Del': '#8c564b', 'Other': '#7f7f7f'}
labels = {'Missense_Mutation': 'Missense', 'Nonsense_Mutation': 'Nonsense', 'Frame_Shift_Del': 'Frame Shift Del',
          'Frame_Shift_Ins': 'Frame Shift Ins', 'Splice_Site': 'Splice Site', 'In_Frame_Del': 'In Frame Del', 'Other': 'Other Non-silent'}

bottom = np.zeros(len(gene_order))
for cat in cats:
    values = ct[cat].values
    ax2.barh(y_pos, values, left=bottom, label=labels[cat], color=palette[cat], height=0.7)
    bottom += values

ax2.set_yticks(y_pos)
ax2.set_yticklabels([])
ax2.invert_yaxis()
ax2.set_xlabel('Non-Silent Mutation Count (n = 1,744)', fontsize=11, fontweight='bold')
ax2.set_title('B. Functional Consequence Spectrum of Non-Silent Mutations', fontsize=12, fontweight='bold', pad=12)
ax2.legend(loc='lower right', frameon=True, fontsize=9.5)
ax2.grid(axis='x', linestyle='--', alpha=0.5)

plt.tight_layout()
os.makedirs(os.path.dirname(out_fig), exist_ok=True)
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
print('Figure 2 saved to:', out_fig)
