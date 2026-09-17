# Comprehensive Pan-Cancer Analysis of Somatic Mutations in Core Circadian Clock Genes

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Cohort: TCGA Pan-Cancer](https://img.shields.io/badge/Cohort-TCGA%20(n%3D10%2C239)-green.svg)]()
[![Genes: Core Clock (16)](https://img.shields.io/badge/Gene%20Set-16%20Core%20Clock%20Genes-orange.svg)]()
[![Genome Build](https://img.shields.io/badge/Genome-GRCh38%2Fhg38-blue.svg)]()

This repository provides the complete, fully reproducible computational pipeline, statistical models, figure reproduction scripts, and clinical survival analytics for the peer-reviewed manuscript:

> **"Pan-Cancer Somatic Mutation Landscape and Clinical Prognostic Relevance of Core Circadian Clock Genes Across the TCGA Cohort"**  
> *Author:* İlkay Civelek  
> *Journal:* Current Issues in Molecular Biology (CIMB), MDPI

---

## 📌 Master Cohort and Mutational Metrics (Verified & Synchronized)

All scripts and datasets in this repository are cross-verified and aligned 1:1 with the published manuscript, tables, and figures:

* **10,239** primary tumor profiles across **33 TCGA cancer types** (uniform MC3 consensus calls).
* **1,294** unique tumor samples harboring at least one somatic mutation in the 16 core clock genes (~12.6% pan-cancer prevalence).
* **2,503** total somatic alteration events across the 16 clock genes.
* **595** synonymous / silent mutations (23.8%).
* **1,908** non-silent somatic mutations (**Table 1**).
* **1,838** protein-coding non-silent alterations (**Figure 2B**).
* **1,532** missense substitutions (61.2% of total, 80.3% of non-silent; **Supplementary Table S1**).
* **706** consensus high-confidence deleterious missense mutations prioritized by dual SIFT (< 0.05) and PolyPhen-2 (> 0.45) prediction (**Figure 7A**).
* **9,850** clinically annotated patients with validated TCGA Clinical Data Resource (TCGA-CDR) endpoints (1,092 clock-mutant vs. 8,758 wild-type; **Table 3**, **Figure 8**).

---

## 🔬 Software Requirements & Computational Environment

### R (>= 4.2.0)
* `survival` (v3.5-7)
* `survminer` (v0.4.9)
* `maftools` (v2.18.0)
* `logistf` (v1.26.0)
* `data.table` (v1.14.8)

### Python (>= 3.9)
* `pandas` (>= 2.0.0)
* `numpy` (>= 1.24.0)
* `matplotlib` (>= 3.7.0)
* `scipy` (>= 1.10.0)

Install Python dependencies via:
```bash
pip install -r requirements.txt
```

---


---

## 🎨 Figure Reproduction & Software Specification

| Item | Description | Software / Tool | Output File |
|:---|:---|:---|:---|
| **Figure 1** | Translational 4-Phase Analytical Workflow | Python (matplotlib) | `Figure 1.png` |
| **Figure 2** | Consolidated Mutational Landscape (16 Genes & Consequences) | Python (matplotlib, pandas) | `Figure 2.png` |
| **Figure 3** | Pan-Cancer Oncoplot Waterfall | R (maftools) | `Figure 3.png` |
| **Figure 4** | Ti/Tv Substitution Spectrum & VAF Distribution | R (maftools) / Python (numpy) | `Figure 4.png` |
| **Figure 5** | Pairwise Co-Occurrence & Mutual Exclusivity Heatmap | Python (matplotlib, scipy) | `Figure 5.png` |
| **Figure 6** | Protein Structural Domain & Hotspot Lollipop Mapping | R (maftools) | `Figure 6.png` |
| **Figure 7** | SIFT/PolyPhen Consensus Overlap & Recurrent Hotspots | Python (matplotlib-venn) | `Figure 7.png` |
| **Figure 8** | Kaplan-Meier Survival Curves & Firth Penalized Forest Plot | Python (matplotlib) / R (survival) | `Figure 8.png` |
| **Figure S1**| Cohort-Specific Alteration Frequency Across 33 TCGA Lineages | Python (matplotlib) | `Figure S1.png` |
| **Table 1** | Somatic Variant Classification Breakdown | Python (pandas) | `Table 1.csv` |
| **Table 2** | Recurrent Deleterious Hotspots with Domain Mapping | Python (pandas) | `Table 2.csv` |
| **Table 3** | Overall Survival & Progression-Free Interval Analytics | Python (pandas, scipy) | `Table 3.csv` |
| **Table S1** | Master Catalog of 1,532 Missense Somatic Variants | Python (openpyxl) | `Table S1.xlsx` |
| **Table S2** | Prioritized Recurrent & Functional Hotspot Alterations | Python (openpyxl) | `Table S2.xlsx` |
| **Table S3** | Pairwise Co-Occurrence Matrix Across 120 Gene Pairs | Python (openpyxl) | `Table S3.xlsx` |

---

## 🚀 One-Click Reproduction Guide

To regenerate all primary figures, analytical tables, and supplementary datasets, run the following standalone Python scripts from the root directory:

### 1. Workflow Architecture (Figure 1)
Generates the publication-grade 4-phase translational workflow in both compact horizontal and 2x2 grid pastel formats:
```bash
python scripts/30_make_Figure1_pipeline_workflow.py
python scripts/30B_make_Figure1_2x2_workflow.py
```
*Outputs:* `figures/Figure_1_Pipeline_Workflow.png`, `figures/Figure_1_Pipeline_Workflow_2X2.png`

### 2. Consolidated Mutational Landscape (Figure 2 & Table 1)
Renders the pan-cancer mutational frequency ranking across 16 genes (Panel A, n=2,503) and functional consequence spectrum of coding alterations (Panel B, n=1,838):
```bash
python scripts/29_make_Figure2_consolidated_landscape.py
```
*Outputs:* `figures/Figure_2_Consolidated_Mutational_Landscape.png`

### 3. In Silico Pathogenicity & Hotspot Recurrence (Figure 7 & Table 2)
Computes SIFT and PolyPhen-2 dual-consensus overlap (706 high-confidence variants) and identifies stereochemical hotspots (PER3 p.R316C, CSNK1E p.R127W, CLOCK p.G120V, CRY1 p.R348C, NR1D2 p.S406L):
```bash
python scripts/31_make_Figure7_insilico_pathogenicity_hotspots.py
```
*Outputs:* `figures/Figure_7_InSilico_Pathogenicity_and_Hotspots_Fixed.png`, `results/Table2_Recurrent_Variants_Annotated.csv`

### 4. Clinical Survival Modeling & TMB Confounder Adjustment (Figure 8 & Table 3)
Calculates Kaplan-Meier overall survival (OS) and progression-free interval (PFI) curves for UCEC, BLCA, and KIRC, generates the multi-cohort Cox HR forest plot across 10 solid tumors, and performs Firth penalized logistic regression adjusting for log1p non-clock TMB:
```bash
python scripts/32_pan_cancer_survival_and_TMB_adjustment_Figure8.py
```
*Outputs:* `figures/Figure_8_Survival_and_TMB_Adjustment.png`, `results/Table3_Survival_Analytics.csv`

### 5. Supplementary Tables Generation (Tables S1, S2, S3)
Extracts all 1,532 missense variants with pathogenicity scores, compiles hotspot domain coordinates, and performs 120 pairwise permutation-based interaction tests:
```bash
python scripts/33_generate_supplementary_tables_S1_S2_S3.py
```
*Outputs:* `supplementary/Supplementary_Table_S1_1532_Missense_Variants.csv`, `supplementary/Supplementary_Table_S2_Functional_Domain_Hotspots.csv`, `supplementary/Supplementary_Table_S3_Gene_Interactions_CoOccurrence.csv`

---

## 📁 Repository Structure

```text
circadian_pan_cancer_Q2/
├── data/
│   ├── circadian_all_mutations_extracted.csv  # Curated GRCh38 mutation dataset (n=2,503)
│   ├── TCGA_PanCancer_CDR_Survival.tsv        # Curated TCGA-CDR clinical survival file (n=9,850)
│   ├── clock_genes_16.txt                     # Official 16 core clock gene symbols
│   └── cancer_abbreviations.csv               # 33 TCGA cancer type mapping
├── figures/                                   # High-resolution (300 DPI) publication figures
├── results/                                   # Table 2 and Table 3 output CSVs
├── scripts/                                   # Modular, sequentially numbered reproduction scripts
├── supplementary/                             # Supplementary Tables S1, S2, and S3
└── README.md                                  # Complete reproduction documentation
```

---

## 📜 Citation & Reference
If you utilize this pipeline, codebase, or curated circadian dataset, please cite:
```bibtex
@article{civelek2026circadian,
  title={Pan-Cancer Somatic Mutation Landscape and Clinical Prognostic Relevance of Core Circadian Clock Genes Across the TCGA Cohort},
  author={Civelek, {\.I}lkay},
  journal={Current Issues in Molecular Biology},
  year={2026},
  publisher={MDPI}
}
```
