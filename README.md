# Comprehensive Pan-Cancer Analysis of Somatic Mutations in Core Circadian Clock Genes

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Cohort: TCGA Pan-Cancer](https://img.shields.io/badge/Cohort-TCGA%20(n%3D10%2C239)-green.svg)]()
[![Genes: Core Clock (16)](https://img.shields.io/badge/Gene%20Set-16%20Core%20Clock%20Genes-orange.svg)]()

This repository contains the end-to-end computational pipeline, statistical models, figure reproduction scripts, and clinical survival analytics for the manuscript:

> **"Pan-Cancer Somatic Mutation Landscape and Clinical Prognostic Relevance of Core Circadian Clock Genes Across the TCGA Cohort"**  
> *Author:* İlkay Civelek

---

## 🔬 Pipeline Overview & Execution Workflow

All scripts are located in the  directory and are numbered sequentially according to analytical phase:

1. **Quality Control & Mutation Calling:**
   -  - Quality control of MC3 somatic mutation calls.
   -  - Extraction of 16 core clock gene coordinates.
   - – - Patient-level QC, denominator validation, and mutational prevalence estimation.

2. **Landscape & Oncoplot Analytics:**
   - – - Gene-level cancer enrichment, consequence distribution, and burden adjustments.
   - – - Pairwise co-occurrence / mutual exclusivity analyses.
   -  - Consolidated Pan-Cancer mutational landscape (Panels A & B, n=2,503 variants).

3. **Structural Mapping & In Silico Pathogenicity:**
   -  - Domain-level mapping for TIMELESS, PER1, CLOCK, and ARNTL.
   -  - SIFT/PolyPhen-2 consensus intersection (706 high-confidence variants) and hotspot recurrence (Table 2).

4. **Clinical Survival & TMB Confounder Modeling:**
   -  - TCGA-CDR overall survival (OS), progression-free interval (PFI) Kaplan-Meier modeling, univariable Cox proportional hazards across 10 major solid tumors, and Firth bias-reduced penalized logistic regression adjusting for non-clock somatic TMB (Figure 8, Table 3).

5. **Supplementary Data Generation:**
   -  - Generates Supplementary Tables S1 (1,532 missense variants), S2 (recurrent & domain hotspots), and S3 (120 pairwise gene interactions).

---

## 📊 Repository Structure



## 🛠️ Requirements & Environment

- **R (>= 4.2.0):** , , , , 
- **Python (>= 3.9):** , , , , , , 
