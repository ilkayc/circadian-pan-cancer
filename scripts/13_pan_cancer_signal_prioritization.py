import pandas as pd
import numpy as np

INFILE = "results/clock16/REAL_CLOCK16_gene_cancer_enrichment.tsv"

OUT_GENE = "results/clock16/PAN_CANCER_GENE_PRIORITIZATION.tsv"
OUT_CANCER = "results/clock16/PAN_CANCER_CANCER_PRIORITIZATION.tsv"
OUT_ASSOC = "results/clock16/PAN_CANCER_STRONGEST_ASSOCIATIONS.tsv"
OUT_MATRIX = "results/clock16/PAN_CANCER_GENE_CANCER_OR_MATRIX.tsv"

print("\n========================================")
print("PAN-CANCER SIGNAL PRIORITIZATION")
print("========================================")

df = pd.read_csv(INFILE, sep="\t")

print("Input rows:", len(df))
print("Columns:", list(df.columns))

# ------------------------------------------------------------
# CHECK
# ------------------------------------------------------------

required = [
    "Cancer_Type",
    "Gene",
    "Cancer_Patients",
    "Gene_Mutated_Cancer_Patients",
    "Cancer_Prevalence_Percent",
    "Other_Prevalence_Percent",
    "Odds_Ratio",
    "Fisher_P",
    "Direction",
    "FDR"
]

missing = [x for x in required if x not in df.columns]

if missing:
    raise ValueError(
        f"Missing required columns: {missing}"
    )

print("Genes:", df["Gene"].nunique())
print("Cancer types:", df["Cancer_Type"].nunique())

# ------------------------------------------------------------
# REMOVE UCEC FOR INDEPENDENT ANALYSIS
# ------------------------------------------------------------

df = df[
    df["Cancer_Type"] !=
    "uterine corpus endometrioid carcinoma"
].copy()

print("\nUCEC removed")
print("Remaining rows:", len(df))
print("Remaining cancer types:", df["Cancer_Type"].nunique())

# ------------------------------------------------------------
# SIGNIFICANT ASSOCIATIONS
# ------------------------------------------------------------

sig = df[df["FDR"] < 0.05].copy()

print("\n========================================")
print("SIGNIFICANT ASSOCIATIONS")
print("========================================")

print("Total significant:", len(sig))

print("\nDirection:")
print(sig["Direction"].value_counts())

# ------------------------------------------------------------
# SIGNAL SCORE
# ------------------------------------------------------------

sig["minus_log10_FDR"] = (
    -np.log10(sig["FDR"].clip(lower=1e-300))
)

sig["log2_OR"] = np.log2(
    sig["Odds_Ratio"].clip(lower=1e-6)
)

sig["Signal_Score"] = (
    sig["minus_log10_FDR"] *
    sig["log2_OR"].abs()
)

# ------------------------------------------------------------
# GENE LEVEL
# ------------------------------------------------------------

gene_summary = (
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

        Mean_Enriched_OR=(
            "Odds_Ratio",
            lambda x:
            x[
                sig.loc[x.index, "Direction"] ==
                "Enriched"
            ].mean()
            if (
                sig.loc[x.index, "Direction"] ==
                "Enriched"
            ).any()
            else np.nan
        ),

        Strongest_Signal=("Signal_Score", "max")
    )
    .reset_index()
)

gene_summary["Enrichment_Consistency"] = (
    gene_summary["Enriched_Cancer_Count"] -
    gene_summary["Depleted_Cancer_Count"]
).abs()

# ------------------------------------------------------------
# CANCER LEVEL
# ------------------------------------------------------------

cancer_summary = (
    sig.groupby("Cancer_Type")
    .agg(
        Significant_Gene_Count=("Gene", "nunique"),

        Enriched_Gene_Count=(
            "Direction",
            lambda x: (x == "Enriched").sum()
        ),

        Depleted_Gene_Count=(
            "Direction",
            lambda x: (x == "Depleted").sum()
        ),

        Best_FDR=("FDR", "min"),

        Max_OR=("Odds_Ratio", "max"),

        Strongest_Signal=("Signal_Score", "max")
    )
    .reset_index()
)

cancer_summary["Net_Enrichment"] = (
    cancer_summary["Enriched_Gene_Count"] -
    cancer_summary["Depleted_Gene_Count"]
)

# ------------------------------------------------------------
# TOP GENES
# ------------------------------------------------------------

print("\n========================================")
print("TOP RECURRENT GENES")
print("========================================")

print(
    gene_summary
    .sort_values(
        [
            "Enriched_Cancer_Count",
            "Significant_Cancer_Count",
            "Best_FDR"
        ],
        ascending=[False, False, True]
    )
    .to_string(index=False)
)

# ------------------------------------------------------------
# TOP CANCERS
# ------------------------------------------------------------

print("\n========================================")
print("TOP CANCER TYPES")
print("========================================")

print(
    cancer_summary
    .sort_values(
        [
            "Enriched_Gene_Count",
            "Significant_Gene_Count",
            "Best_FDR"
        ],
        ascending=[False, False, True]
    )
    .to_string(index=False)
)

# ------------------------------------------------------------
# TOP INDIVIDUAL ASSOCIATIONS
# ------------------------------------------------------------

print("\n========================================")
print("TOP INDIVIDUAL ASSOCIATIONS")
print("========================================")

cols = [
    "Cancer_Type",
    "Gene",
    "Cancer_Patients",
    "Gene_Mutated_Cancer_Patients",
    "Cancer_Prevalence_Percent",
    "Other_Prevalence_Percent",
    "Odds_Ratio",
    "Fisher_P",
    "FDR",
    "Direction",
    "Signal_Score"
]

print(
    sig
    .sort_values(
        "Signal_Score",
        ascending=False
    )[cols]
    .head(50)
    .to_string(index=False)
)

# ------------------------------------------------------------
# RECURRENT ENRICHED GENES
# ------------------------------------------------------------

print("\n========================================")
print("RECURRENT ENRICHED GENES")
print("========================================")

recurrent_enriched = gene_summary[
    gene_summary["Enriched_Cancer_Count"] >= 2
].copy()

print(
    recurrent_enriched
    .sort_values(
        [
            "Enriched_Cancer_Count",
            "Best_FDR"
        ],
        ascending=[False, True]
    )
    .to_string(index=False)
)

# ------------------------------------------------------------
# CANCERS WITH MULTIPLE ENRICHED CLOCK GENES
# ------------------------------------------------------------

print("\n========================================")
print("CANCERS WITH MULTIPLE ENRICHED CLOCK GENES")
print("========================================")

multi_gene_cancers = cancer_summary[
    cancer_summary["Enriched_Gene_Count"] >= 2
].copy()

print(
    multi_gene_cancers
    .sort_values(
        [
            "Enriched_Gene_Count",
            "Best_FDR"
        ],
        ascending=[False, True]
    )
    .to_string(index=False)
)

# ------------------------------------------------------------
# GENE × CANCER MATRIX
# ------------------------------------------------------------

matrix = sig.pivot_table(
    index="Cancer_Type",
    columns="Gene",
    values="Odds_Ratio",
    aggfunc="first"
)

matrix.to_csv(
    OUT_MATRIX,
    sep="\t"
)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

gene_summary.to_csv(
    OUT_GENE,
    sep="\t",
    index=False
)

cancer_summary.to_csv(
    OUT_CANCER,
    sep="\t",
    index=False
)

sig.sort_values(
    "Signal_Score",
    ascending=False
).to_csv(
    OUT_ASSOC,
    sep="\t",
    index=False
)

print("\n========================================")
print("SAVED")
print("========================================")

print(OUT_GENE)
print(OUT_CANCER)
print(OUT_ASSOC)
print(OUT_MATRIX)
