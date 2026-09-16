import pandas as pd
import numpy as np
from pathlib import Path

PROJECT = Path.home() / "circadian_pan_cancer_Q2"

INFILE = (
    PROJECT
    / "results/clock16/REAL_CLOCK16_gene_cancer_enrichment.tsv"
)

OUTFILE = (
    PROJECT
    / "results/clock16/UCEC_independent_recurrent_gene_associations.tsv"
)

# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INFILE, sep="\t")

print("\n========================================")
print("UCEC-INDEPENDENT RECURRENCE ANALYSIS")
print("========================================")

print("Input rows:", len(df))
print("Genes:", df["Gene"].nunique())
print("Cancer types:", df["Cancer_Type"].nunique())

# ============================================================
# REMOVE UCEC
# ============================================================

ucec_name = "uterine corpus endometrioid carcinoma"

df_no_ucec = df[
    df["Cancer_Type"].str.lower() != ucec_name.lower()
].copy()

print("\nUCEC rows removed:", len(df) - len(df_no_ucec))
print("Remaining rows:", len(df_no_ucec))
print("Remaining cancer types:", df_no_ucec["Cancer_Type"].nunique())

assert len(df_no_ucec) == 31 * 16
assert df_no_ucec["Cancer_Type"].nunique() == 31
assert df_no_ucec["Gene"].nunique() == 16

# ============================================================
# SIGNIFICANT ASSOCIATIONS
# ============================================================

sig = df_no_ucec[
    (df_no_ucec["FDR"] < 0.05)
].copy()

sig["Direction"] = sig["Direction"].astype(str)

# ============================================================
# RECURRENT PER-GENE ASSOCIATIONS
# ============================================================

summary = (
    sig.groupby("Gene")
    .agg(
        Significant_Cancer_Count=("Cancer_Type", "nunique"),
        Enriched_Cancer_Count=(
            "Direction",
            lambda x: (x == "Enriched").sum()
        ),
        Depleted_Cancer_Count=(
            "Direction",
            lambda x: (x == "Depleted").sum()
        ),
        Best_FDR=("FDR", "min"),
        Max_OR=("Odds_Ratio", "max"),
    )
    .reset_index()
)

# Add all 16 genes
all_genes = pd.DataFrame({
    "Gene": sorted(df_no_ucec["Gene"].unique())
})

summary = all_genes.merge(
    summary,
    on="Gene",
    how="left"
)

for col in [
    "Significant_Cancer_Count",
    "Enriched_Cancer_Count",
    "Depleted_Cancer_Count",
]:
    summary[col] = summary[col].fillna(0).astype(int)

summary["Best_FDR"] = summary["Best_FDR"].fillna(1.0)
summary["Max_OR"] = summary["Max_OR"].fillna(1.0)

summary = summary.sort_values(
    [
        "Significant_Cancer_Count",
        "Enriched_Cancer_Count",
        "Best_FDR"
    ],
    ascending=[False, False, True]
)

# ============================================================
# RECURRENT ASSOCIATIONS
# ============================================================

# At least 2 independent cancer types
recurrent = summary[
    summary["Significant_Cancer_Count"] >= 2
].copy()

# ============================================================
# SAVE
# ============================================================

summary.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

# ============================================================
# PRINT
# ============================================================

print("\n========================================")
print("UCEC-INDEPENDENT SIGNIFICANT RESULTS")
print("========================================")

print("Total significant gene × cancer associations:",
      len(sig))

print("\nSignificant associations by direction:")
print(sig["Direction"].value_counts())

print("\n========================================")
print("GENE-LEVEL RECURRENCE")
print("========================================")

print(summary.to_string(index=False))

print("\n========================================")
print("RECURRENT GENES")
print("(>=2 independent cancer types)")
print("========================================")

if len(recurrent) == 0:
    print("No gene has significant associations in >=2 cancer types.")
else:
    print(recurrent.to_string(index=False))

# ============================================================
# DETAILED RECURRING ASSOCIATIONS
# ============================================================

print("\n========================================")
print("DETAILS OF RECURRENT GENES")
print("========================================")

for gene in recurrent["Gene"]:

    sub = sig[
        sig["Gene"] == gene
    ].sort_values("FDR")

    print(f"\n--- {gene} ---")

    print(
        sub[
            [
                "Cancer_Type",
                "Cancer_Patients",
                "Gene_Mutated_Cancer_Patients",
                "Cancer_Prevalence_Percent",
                "Other_Prevalence_Percent",
                "Odds_Ratio",
                "Fisher_P",
                "FDR",
                "Direction",
            ]
        ].to_string(index=False)
    )

print("\n========================================")
print("SAVED")
print("========================================")

print(OUTFILE)
