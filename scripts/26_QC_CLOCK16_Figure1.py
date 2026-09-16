import os
import pandas as pd
import numpy as np

PROJECT = os.path.expanduser("~/circadian_pan_cancer_Q2")
QC_DIR = os.path.expanduser("~/Downloads/CLOCK16_FINAL/QC")
os.makedirs(QC_DIR, exist_ok=True)

REPORT = os.path.join(QC_DIR, "Figure1_QC_report.txt")

lines = []

def log(x=""):
    print(x)
    lines.append(str(x))

def save_qc(df, name):
    path = os.path.join(QC_DIR, name)
    df.to_csv(path, sep="\t", index=False)
    log(f"Saved: {path}")

log("=" * 90)
log("CLOCK16 FIGURE 1 — FINAL SOURCE QC")
log("=" * 90)

# ============================================================
# A — PREVALENCE
# ============================================================

A_path = os.path.join(
    PROJECT,
    "results/clock16/CLOCK16_final_cancer_prevalence.tsv"
)

A = pd.read_csv(A_path, sep="\t")

A_qc = A.copy()

log("\n[FIGURE 1A — PREVALENCE]")
log(f"Rows: {len(A)}")
log(f"Total patients: {A['N'].sum()}")
log(f"Cancer types: {A['Cancer_Type'].nunique()}")
log(f"Total CLOCK16-mutated patients: {A['CLOCK16_Mutated'].sum()}")

if "Prevalence_Percent" in A.columns:
    log(
        f"Prevalence range: "
        f"{A['Prevalence_Percent'].min():.6f} - "
        f"{A['Prevalence_Percent'].max():.6f}%"
    )

save_qc(A_qc, "Figure1A_QC.tsv")

# ============================================================
# B — CANCER ASSOCIATION
# ============================================================

B_path = os.path.join(
    PROJECT,
    "results/clock16/REAL_CLOCK16_cancer_association_ALL_v2.tsv"
)

B = pd.read_csv(B_path, sep="\t")

required_B = [
    "phenotype",
    "cancer_type",
    "total_primary_patients",
    "mutated_patients",
    "nonmutated_patients",
    "prevalence_percent",
    "odds_ratio",
    "lower_95CI",
    "upper_95CI",
    "fisher_p",
    "FDR"
]

missing_B = [x for x in required_B if x not in B.columns]

if missing_B:
    raise ValueError(f"Figure 1B missing columns: {missing_B}")

B = B[B["phenotype"].astype(str) == "Any_CLOCK16"].copy()

B["Direction"] = np.where(
    B["odds_ratio"] > 1,
    "Enriched",
    "Depleted"
)

log("\n[FIGURE 1B — CANCER ASSOCIATION]")
log(f"Rows: {len(B)}")
log(f"Cancer types: {B['cancer_type'].nunique()}")
log(f"FDR < 0.05: {(B['FDR'] < 0.05).sum()}")
log(f"FDR >= 0.05: {(B['FDR'] >= 0.05).sum()}")
log(f"Enriched OR > 1: {(B['odds_ratio'] > 1).sum()}")
log(f"Depleted OR < 1: {(B['odds_ratio'] < 1).sum()}")
log(
    f"OR range: "
    f"{B['odds_ratio'].min():.6f} - "
    f"{B['odds_ratio'].max():.6f}"
)

save_qc(B, "Figure1B_QC.tsv")

# ============================================================
# C — GENE LEVEL FIRTH
# ============================================================

C_path = os.path.join(
    PROJECT,
    "results/clock16/CLOCK16_gene_level_pooled_FIRTH_LR.tsv"
)

C = pd.read_csv(C_path, sep="\t")

required_C = [
    "Gene",
    "N",
    "Mutated_Patients",
    "Wild_Type_Patients",
    "OR_per_log1p_non_CLOCK16_mutation",
    "CI95_Low",
    "CI95_High",
    "Firth_LR_P",
    "FDR"
]

missing_C = [x for x in required_C if x not in C.columns]

if missing_C:
    raise ValueError(f"Figure 1C missing columns: {missing_C}")

log("\n[FIGURE 1C — GENE LEVEL FIRTH]")
log(f"Genes: {len(C)}")
log(f"Unique genes: {C['Gene'].nunique()}")
log(f"FDR < 0.05: {(C['FDR'] < 0.05).sum()}")
log(f"FDR >= 0.05: {(C['FDR'] >= 0.05).sum()}")
log(
    f"OR range: "
    f"{C['OR_per_log1p_non_CLOCK16_mutation'].min():.6f} - "
    f"{C['OR_per_log1p_non_CLOCK16_mutation'].max():.6f}"
)

save_qc(C, "Figure1C_QC.tsv")

# ============================================================
# D — FUNCTIONAL FIRTH
# ============================================================

D_path = os.path.join(
    PROJECT,
    "results/clock16/CLOCK16_functional_CORRECTED_pooled_FIRTH_LR.tsv"
)

D = pd.read_csv(D_path, sep="\t")

required_D = [
    "Endpoint",
    "N",
    "Cases",
    "Controls",
    "OR_per_log1p_non_CLOCK16_mutation",
    "CI95_Low",
    "CI95_High",
    "Firth_LR_P",
    "FDR"
]

missing_D = [x for x in required_D if x not in D.columns]

if missing_D:
    raise ValueError(f"Figure 1D missing columns: {missing_D}")

log("\n[FIGURE 1D — FUNCTIONAL]")
log(f"Endpoints: {len(D)}")
log(f"FDR < 0.05: {(D['FDR'] < 0.05).sum()}")
log(f"FDR >= 0.05: {(D['FDR'] >= 0.05).sum()}")

save_qc(D, "Figure1D_QC.tsv")

# ============================================================
# CROSS-PANEL CONSISTENCY
# ============================================================

log("\n" + "=" * 90)
log("CROSS-PANEL CONSISTENCY")
log("=" * 90)

Ns = {
    "A_total": int(A["N"].sum()),
    "A_types": int(len(A)),
    "C": int(C["N"].iloc[0]),
    "D": int(D["N"].iloc[0])
}

log(f"Panel A total N: {Ns['A_total']}")
log(f"Panel A cancer types: {Ns['A_types']}")
log(f"Panel C N: {Ns['C']}")
log(f"Panel D N: {Ns['D']}")

if Ns["A_total"] == Ns["C"] == Ns["D"]:
    log("PASS: A/C/D total cohort N values are consistent.")
else:
    log("WARNING: A/C/D total cohort N values are NOT identical.")

# Expected final values based on current locked results

checks = {
    "A cancer types == 32": len(A) == 32,
    "A CLOCK16 mutated == 1165": int(A["CLOCK16_Mutated"].sum()) == 1165,
    "B cancer types == 32": len(B) == 32,
    "B FDR significant == 23": int((B["FDR"] < 0.05).sum()) == 23,
    "B enriched == 12": int((B["odds_ratio"] > 1).sum()) == 12,
    "B depleted == 20": int((B["odds_ratio"] < 1).sum()) == 20,
    "C genes == 16": len(C) == 16,
    "C FDR significant == 16": int((C["FDR"] < 0.05).sum()) == 16,
    "D endpoints == 3": len(D) == 3,
    "D FDR significant == 3": int((D["FDR"] < 0.05).sum()) == 3,
}

log("\nFINAL EXPECTED-VALUE CHECKS")

all_pass = True

for name, result in checks.items():
    status = "PASS" if result else "FAIL"
    log(f"{status}: {name}")
    if not result:
        all_pass = False

log("\n" + "=" * 90)

if all_pass:
    log("FINAL STATUS: PASS")
    log("Figure 1 source tables are internally consistent with the locked results.")
else:
    log("FINAL STATUS: FAIL")
    log("One or more expected-value checks failed. DO NOT lock Figure 1 yet.")

log("=" * 90)

with open(REPORT, "w") as f:
    f.write("\n".join(lines))

print("\nQC REPORT:")
print(REPORT)
print("\nQC COMPLETE.")
