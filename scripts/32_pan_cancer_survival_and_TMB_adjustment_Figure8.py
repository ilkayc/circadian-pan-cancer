#!/usr/bin/env python3
"""
Script 32: Pan-Cancer Clinical Survival and Multivariable TMB-Adjusted Cox Modeling (Table 3 & Figure 8)
Evaluates TCGA-CDR overall survival (OS), univariable Cox proportional hazards HRs,
and multivariable Cox models adjusted for log1p-transformed tumor mutational burden [log(TMB + 1)].
"""
import os
import gzip
import collections
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(BASE_DIR, 'data')
res_dir = os.path.join(BASE_DIR, 'results')
os.makedirs(res_dir, exist_ok=True)

print('Calculating per-sample total mutational burden (TMB)...')
tmb_dict = collections.Counter()
with gzip.open(os.path.join(data_dir, 'mc3.v0.2.8.PUBLIC.xena.gz'), 'rt') as f:
    header = f.readline().strip().split('	')
    s_idx = header.index('sample')
    for line in f:
        s = line.split('	')[s_idx][:12]
        tmb_dict[s] += 1

# Clock mutant patients
df_clock = pd.read_csv(os.path.join(data_dir, 'circadian_all_mutations_extracted.csv'))
clock_patients = set(df_clock['Tumor_Sample_Barcode'].str[:12].unique())

# Clinical survival data
df_surv = pd.read_csv(os.path.join(data_dir, 'TCGA_PanCancer_CDR_Survival.tsv'), sep='	')
df_surv['patient'] = df_surv['_PATIENT'].str[:12]
df_surv['clock_mut'] = df_surv['patient'].apply(lambda x: 1 if x in clock_patients else 0)
df_surv['TMB'] = df_surv['patient'].map(tmb_dict).fillna(0)
df_surv['log1p_TMB'] = np.log1p(df_surv['TMB'])

cohort_names = {
    'UCEC': 'Uterine Corpus Endometrial Carcinoma',
    'BLCA': 'Bladder Urothelial Carcinoma',
    'STAD': 'Stomach Adenocarcinoma',
    'LUAD': 'Lung Adenocarcinoma',
    'HNSC': 'Head and Neck Squamous Cell Carcinoma',
    'BRCA': 'Breast Invasive Carcinoma',
    'LUSC': 'Lung Squamous Cell Carcinoma',
    'COAD': 'Colon Adenocarcinoma',
    'SKCM': 'Skin Cutaneous Melanoma',
    'KIRC': 'Kidney Renal Clear Cell Carcinoma',
    'PAN-CANCER': 'Combined TCGA Pan-Cancer Cohort'
}

table3_rows = []

for c, full_name in cohort_names.items():
    if c == 'PAN-CANCER':
        sub = df_surv[(df_surv['OS.time'] > 0) & (df_surv['OS'].notna()) & (df_surv['TMB'] > 0)].copy()
    else:
        sub = df_surv[(df_surv['cancer type abbreviation'] == c) & (df_surv['OS.time'] > 0) & (df_surv['OS'].notna()) & (df_surv['TMB'] > 0)].copy()
    
    if len(sub) == 0: continue
    
    # Univariable Cox
    cph_uni = CoxPHFitter()
    cph_uni.fit(sub[['OS.time', 'OS', 'clock_mut']], duration_col='OS.time', event_col='OS')
    u_hr = np.exp(cph_uni.params_['clock_mut'])
    u_low = np.exp(cph_uni.confidence_intervals_.loc['clock_mut'].iloc[0])
    u_high = np.exp(cph_uni.confidence_intervals_.loc['clock_mut'].iloc[1])
    u_p = cph_uni.summary.loc['clock_mut', 'p']
    
    # Multivariable Cox with log1p(TMB)
    cph_multi = CoxPHFitter()
    cph_multi.fit(sub[['OS.time', 'OS', 'clock_mut', 'log1p_TMB']], duration_col='OS.time', event_col='OS')
    m_hr = np.exp(cph_multi.params_['clock_mut'])
    m_low = np.exp(cph_multi.confidence_intervals_.loc['clock_mut'].iloc[0])
    m_high = np.exp(cph_multi.confidence_intervals_.loc['clock_mut'].iloc[1])
    m_p = cph_multi.summary.loc['clock_mut', 'p']
    
    mut_n = int(sub['clock_mut'].sum())
    wt_n = len(sub) - mut_n
    
    table3_rows.append({
        'Cohort': c,
        'Primary Cancer Type': full_name,
        'Total (N)': len(sub),
        'Clock-Mutant (n)': mut_n,
        'Wild-Type (n)': wt_n,
        'Unadjusted OS HR (95% CI)': f'{u_hr:.2f} ({u_low:.2f} - {u_high:.2f})',
        'Unadjusted OS P': f'{u_p:.4f}' if u_p >= 0.0001 else '< 0.0001',
        'TMB-Adjusted OS HR (95% CI)': f'{m_hr:.2f} ({m_low:.2f} - {m_high:.2f})',
        'TMB-Adjusted OS P': f'{m_p:.4f}' if m_p >= 0.0001 else '< 0.0001'
    })

df_t3 = pd.DataFrame(table3_rows)
out_csv = os.path.join(res_dir, 'Table 3.csv')
out_xlsx = os.path.join(res_dir, 'Table 3.xlsx')
df_t3.to_csv(out_csv, index=False)
df_t3.to_excel(out_xlsx, index=False)
print('Script 32 executed successfully! Table 3 saved to:', out_csv)
