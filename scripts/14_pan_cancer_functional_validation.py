#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path("results/clock16")
DEN = Path("results/denominator")

PRIMARY = DEN / "primary_samples_with_patient_id.tsv"
MC3_PRIMARY = DEN / "primary_tumors_in_MC3.tsv"
MUTATIONS = BASE / "REAL_primary_CLOCK16_mutations.tsv.gz"
PATIENT_EFFECT = BASE / "REAL_primary_CLOCK16_patient_effect_class.tsv"
PATIENT_BURDEN = BASE / "REAL_CLOCK16_patient_functional_burden.tsv"
GENE_CANCER = BASE / "REAL_CLOCK16_gene_cancer_enrichment.tsv"


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, sep="\t")


print("=" * 72)
print("CLOCK16 FUNCTIONAL VALIDATION + DATA INTEGRITY")
print("=" * 72)

# ------------------------------------------------------------
# 1. PRIMARY SAMPLE DENOMINATOR
# ------------------------------------------------------------

primary = read_tsv(PRIMARY)

print("\n[1] PRIMARY SAMPLE DENOMINATOR")
print("File:", PRIMARY)
print("Rows:", len(primary))
print("Columns:", list(primary.columns))

required_primary = {"sample", "patient_id"}
missing = required_primary - set(primary.columns)

if missing:
    raise ValueError(
        f"Primary denominator missing required columns: {sorted(missing)}"
    )

print("Unique samples:", primary["sample"].nunique())
print("Unique patients:", primary["patient_id"].nunique())

if "sample_type" in primary.columns:
    print("Sample types:")
    print(primary["sample_type"].value_counts(dropna=False))

if "cancer_type" in primary.columns:
    print("Cancer types:", primary["cancer_type"].nunique())


# ------------------------------------------------------------
# 2. MC3 PRIMARY SAMPLES
# ------------------------------------------------------------

mc3 = read_tsv(MC3_PRIMARY)

print("\n[2] MC3 PRIMARY SAMPLE SET")
print("File:", MC3_PRIMARY)
print("Rows:", len(mc3))
print("Columns:", list(mc3.columns))

print("Unique values by column:")
for c in mc3.columns:
    print(f"  {c}: {mc3[c].nunique(dropna=False)}")


# ------------------------------------------------------------
# 3. CLOCK16 MUTATIONS
# ------------------------------------------------------------

mut = read_tsv(MUTATIONS)

print("\n[3] REAL CLOCK16 MUTATIONS")
print("Rows:", len(mut))
print("Columns:", list(mut.columns))

if "patient_id" not in mut.columns:
    raise ValueError("CLOCK16 mutation file lacks patient_id")

if "gene" not in mut.columns and "Gene" not in mut.columns:
    raise ValueError("CLOCK16 mutation file lacks gene column")

gene_col = "gene" if "gene" in mut.columns else "Gene"

print("Unique mutation patients:", mut["patient_id"].nunique())
print("Unique CLOCK genes:", mut[gene_col].nunique())

print("\nMutation counts by gene:")
print(mut[gene_col].value_counts().to_string())


# ------------------------------------------------------------
# 4. PATIENT EFFECT CLASS
# ------------------------------------------------------------

effect = read_tsv(PATIENT_EFFECT)

print("\n[4] PATIENT EFFECT CLASS")
print("Rows:", len(effect))
print("Columns:", list(effect.columns))

required_effect = {
    "patient_id",
    "Has_Nonsynonymous",
    "Has_Silent",
    "Has_Other",
    "Patient_Class",
}

missing = required_effect - set(effect.columns)

if missing:
    raise ValueError(
        f"Effect-class file missing columns: {sorted(missing)}"
    )

print("Unique patients:", effect["patient_id"].nunique())

print("\nPatient classes:")
print(effect["Patient_Class"].value_counts(dropna=False).to_string())


# ------------------------------------------------------------
# 5. FUNCTIONAL BURDEN
# ------------------------------------------------------------

burden = read_tsv(PATIENT_BURDEN)

print("\n[5] FUNCTIONAL BURDEN")
print("Rows:", len(burden))
print("Columns:", list(burden.columns))

required_burden = {
    "patient_id",
    "Nonsynonymous_Events",
    "Nonsynonymous_Genes",
    "Strong_Damaging_Events",
    "Damaging_Events",
    "Strong_Damaging_Genes",
    "Damaging_Genes",
}

missing = required_burden - set(burden.columns)

if missing:
    raise ValueError(
        f"Functional burden file missing columns: {sorted(missing)}"
    )

print("Unique patients:", burden["patient_id"].nunique())

print("\nBurden summary:")
print(
    burden[
        [
            "Nonsynonymous_Events",
            "Nonsynonymous_Genes",
            "Strong_Damaging_Events",
            "Damaging_Events",
            "Strong_Damaging_Genes",
            "Damaging_Genes",
        ]
    ].describe().round(3).to_string()
)


# ------------------------------------------------------------
# 6. PATIENT-LEVEL INTEGRITY
# ------------------------------------------------------------

print("\n[6] PATIENT-LEVEL INTEGRITY")

primary_patients = set(primary["patient_id"].dropna())
mutation_patients = set(mut["patient_id"].dropna())
effect_patients = set(effect["patient_id"].dropna())
burden_patients = set(burden["patient_id"].dropna())

print("Primary patients:", len(primary_patients))
print("Mutation patients:", len(mutation_patients))
print("Effect-class patients:", len(effect_patients))
print("Burden patients:", len(burden_patients))

print(
    "Mutation patients contained in primary denominator:",
    len(mutation_patients & primary_patients),
    "/",
    len(mutation_patients),
)

print(
    "Effect patients contained in primary denominator:",
    len(effect_patients & primary_patients),
    "/",
    len(effect_patients),
)

print(
    "Burden patients contained in primary denominator:",
    len(burden_patients & primary_patients),
    "/",
    len(burden_patients),
)


# ------------------------------------------------------------
# 7. MERGE PATIENT-LEVEL FUNCTIONAL DATA
# ------------------------------------------------------------

merged = (
    effect.merge(
        burden,
        on="patient_id",
        how="outer",
        indicator="_effect_burden_merge",
    )
)

print("\n[7] EFFECT + BURDEN MERGE")
print("Rows:", len(merged))

print(
    merged["_effect_burden_merge"]
    .value_counts(dropna=False)
    .to_string()
)

merged.drop(columns=["_effect_burden_merge"], inplace=True)


# ------------------------------------------------------------
# 8. FUNCTIONAL FRACTIONS
# ------------------------------------------------------------

n_patients = len(merged)

if n_patients > 0:

    summary = {
        "Patients_in_effect_or_burden": n_patients,
        "Nonsynonymous_patients": int(
            (merged["Has_Nonsynonymous"] == True).sum()
        ),
        "Silent_only_patients": int(
            (merged["Patient_Class"] == "silent_only").sum()
        ),
        "Other_only_patients": int(
            (merged["Patient_Class"] == "other_only").sum()
        ),
        "Strong_damaging_patients": int(
            (merged["Strong_Damaging_Genes"].fillna(0) > 0).sum()
        ),
        "Damaging_patients": int(
            (merged["Damaging_Genes"].fillna(0) > 0).sum()
        ),
    }

    print("\nFunctional patient counts:")
    for k, v in summary.items():
        print(f"{k}: {v}")


# ------------------------------------------------------------
# 9. GENE × CANCER ENRICHMENT VALIDATION
# ------------------------------------------------------------

gc = read_tsv(GENE_CANCER)

print("\n[8] GENE × CANCER ENRICHMENT")
print("Rows:", len(gc))
print("Columns:", list(gc.columns))

required_gc = {
    "Cancer_Type",
    "Gene",
    "Cancer_Patients",
    "Gene_Mutated_Cancer_Patients",
    "Odds_Ratio",
    "Fisher_P",
    "Direction",
    "FDR",
}

missing = required_gc - set(gc.columns)

if missing:
    raise ValueError(
        f"Gene-cancer enrichment missing columns: {sorted(missing)}"
    )

print("Genes:", gc["Gene"].nunique())
print("Cancer types:", gc["Cancer_Type"].nunique())

sig = gc[gc["FDR"] < 0.05].copy()

print("Significant associations:", len(sig))

print("\nDirection:")
print(sig["Direction"].value_counts().to_string())


# ------------------------------------------------------------
# 10. SAVE VALIDATION SUMMARY
# ------------------------------------------------------------

out = BASE / "PAN_CANCER_FUNCTIONAL_VALIDATION.tsv"

validation = pd.DataFrame(
    [
        ["Primary_samples", len(primary)],
        ["Primary_patients", primary["patient_id"].nunique()],
        ["CLOCK16_mutation_patients", mut["patient_id"].nunique()],
        ["Effect_class_patients", effect["patient_id"].nunique()],
        ["Functional_burden_patients", burden["patient_id"].nunique()],
        ["CLOCK16_genes", mut[gene_col].nunique()],
        ["Gene_cancer_rows", len(gc)],
        ["Significant_gene_cancer_associations", len(sig)],
        [
            "Significant_enriched_associations",
            int((sig["Direction"] == "Enriched").sum()),
        ],
        [
            "Significant_depleted_associations",
            int((sig["Direction"] == "Depleted").sum()),
        ],
    ],
    columns=["Metric", "Value"],
)

validation.to_csv(out, sep="\t", index=False)

print("\nSaved:")
print(out)

print("\n" + "=" * 72)
print("VALIDATION COMPLETED")
print("=" * 72)
