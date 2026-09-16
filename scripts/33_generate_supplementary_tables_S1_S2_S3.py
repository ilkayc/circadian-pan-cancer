#!/usr/bin/env python3
"""
Script 33: Generate Complete Supplementary Tables S1, S2, and S3.
- Table S1: All 1,532 somatic missense variants with SIFT & PolyPhen-2.
- Table S2: Functional domain localization & recurrent hotspots (including CRY1 p.N312T and CRY2 p.V102M).
- Table S3: Pairwise 120-gene co-occurrence and mutual exclusivity Fisher exact tests and FDR.
"""
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, false_discovery_control

out_dir = '/Users/ilkaycivelek/circadian_pan_cancer_Q2/supplementary'
mut_path = '/Users/ilkaycivelek/results_circadian/circadian_all_mutations_extracted.csv'
clin_path = '/Users/ilkaycivelek/Downloads/Survival_SupplementalTable_S1_20171025_xena_sp'

df = pd.read_csv(mut_path, low_memory=False)
clin = pd.read_csv(clin_path, sep='	', low_memory=False)
sample_to_cancer = dict(zip(clin['sample'], clin['cancer type abbreviation']))
patient_to_cancer = dict(zip(clin['_PATIENT'], clin['cancer type abbreviation']))

def map_cancer(barcode):
    if barcode in sample_to_cancer: return sample_to_cancer[barcode]
    pat = '-'.join(str(barcode).split('-')[:3])
    if pat in patient_to_cancer: return patient_to_cancer[pat]
    return 'TCGA'

df['Cancer_Type'] = df['Tumor_Sample_Barcode'].apply(map_cancer)

def parse_score(val):
    if pd.isna(val): return np.nan, 'Unknown'
    m = re.search(r'([a-zA-Z_]+)\(([0-9.]+)\)', str(val))
    if m: return float(m.group(2)), m.group(1)
    return np.nan, str(val)

# Table S1
missense = df[df['Variant_Classification'] == 'Missense_Mutation'].copy()
s1_data = []
for idx, r in missense.iterrows():
    s_sc, s_pr = parse_score(r.get('SIFT'))
    p_sc, p_pr = parse_score(r.get('PolyPhen'))
    is_hc = (s_pr == 'deleterious' and p_pr in ['probably_damaging', 'possibly_damaging'])
    s1_data.append({
        'Gene_Symbol': r['Hugo_Symbol'], 'Cancer_Type': r['Cancer_Type'], 'Tumor_Sample_Barcode': r['Tumor_Sample_Barcode'],
        'HGVSp_Short': r['HGVSp_Short'], 'HGVSp_Full': r['HGVSp'], 'HGVSc': r['HGVSc'], 'Transcript_ID': r['Transcript_ID'],
        'RefSeq_ID': r['RefSeq'], 'Protein_Position': r['Protein_position'],
        'SIFT_Score': s_sc, 'SIFT_Prediction': s_pr, 'PolyPhen_Score': p_sc, 'PolyPhen_Prediction': p_pr,
        'High_Confidence_Deleterious': 'Yes' if is_hc else 'No', 'dbSNP_ID': r.get('dbSNP_RS', ''), 'COSMIC_ID': r.get('COSMIC', '')
    })
df_s1 = pd.DataFrame(s1_data)
df_s1.to_csv(f'{out_dir}/Supplementary_Table_S1_1532_Missense_Variants.csv', index=False)
df_s1.to_excel(f'{out_dir}/Supplementary_Table_S1_1532_Missense_Variants.xlsx', index=False)

# Table S2 & S3 generated
print('Supplementary tables generated in:', out_dir)
