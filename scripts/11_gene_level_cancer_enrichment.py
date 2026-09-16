import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ============================================================
# PATHS
# ============================================================

primary_file = "results/denominator/primary_samples_with_patient_id.tsv"
phenotype_file = "data/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
mutation_file = "results/clock16/clock16_mutations.tsv.gz"

out_file = "results/clock16/REAL_CLOCK16_gene_cancer_enrichment.tsv"

# ============================================================
# LOAD
# ============================================================

primary = pd.read_csv(
    primary_file,
    sep="\t",
    low_memory=False
)

phenotype = pd.read_csv(
    phenotype_file,
    sep="\t",
    compression="gzip",
    low_memory=False
)

clock = pd.read_csv(
    mutation_file,
    sep="\t",
    compression="gzip",
    low_memory=False
)

print("\n========================================")
print("GENE × CANCER ENRICHMENT")
print("========================================")

print("Primary columns:", primary.columns.tolist())
print("Phenotype columns:", phenotype.columns.tolist())

# ============================================================
# CREATE PATIENT ID IN PRIMARY IF NECESSARY
# ============================================================

if "patient_id" not in primary.columns:

    primary["patient_id"] = (
        primary["sample"]
        .astype(str)
        .str.extract(r"^(TCGA-[^-]+-[^-]+)")[0]
    )

# ============================================================
# PHENOTYPE: SAMPLE → CANCER TYPE
# ============================================================

phenotype = phenotype[
    ["sample", "_primary_disease"]
].drop_duplicates()

phenotype = phenotype.rename(
    columns={"_primary_disease": "cancer_type"}
)

# ============================================================
# ATTACH CANCER TYPE TO PRIMARY SAMPLES
# ============================================================

primary = primary.merge(
    phenotype,
    on="sample",
    how="left",
    suffixes=("", "_pheno")
)

# If primary already had cancer_type, use phenotype value
if "cancer_type_pheno" in primary.columns:

    primary["cancer_type"] = primary["cancer_type_pheno"].combine_first(
        primary["cancer_type"]
    )

    primary = primary.drop(
        columns=["cancer_type_pheno"]
    )

# ============================================================
# QC
# ============================================================

print("\nPhenotype primary samples:", len(primary))

missing_cancer = primary["cancer_type"].isna().sum()

print(
    "Primary samples without cancer type:",
    missing_cancer
)

if missing_cancer > 0:
    print(
        "\nERROR: Some primary samples have no cancer type."
    )
    print(
        primary.loc[
            primary["cancer_type"].isna(),
            "sample"
        ].head(20)
    )
    raise SystemExit(1)

# One patient = one primary sample
primary_patients = (
    primary[
        ["patient_id", "cancer_type"]
    ]
    .drop_duplicates()
)

print(
    "Total primary patients with cancer type:",
    primary_patients["patient_id"].nunique()
)

# ============================================================
# CLOCK16 PATIENT ID
# ============================================================

clock["patient_id"] = (
    clock["sample"]
    .astype(str)
    .str.extract(r"^(TCGA-[^-]+-[^-]+)")[0]
)

# ============================================================
# CLOCK16 GENES
# ============================================================

clock_genes = [
    "TIMELESS",
    "PER2",
    "PER3",
    "PER1",
    "NPAS2",
    "RORB",
    "CLOCK",
    "RORC",
    "RORA",
    "ARNTL",
    "CSNK1E",
    "CRY1",
    "CRY2",
    "CSNK1D",
    "NR1D2",
    "NR1D1"
]

clock = clock[
    clock["gene"].isin(clock_genes)
].copy()

# ============================================================
# PRIMARY CLOCK16 ONLY
# ============================================================

# Use sample type code from TCGA barcode
clock["sample_type_code"] = (
    clock["sample"]
    .astype(str)
    .str[13:15]
)

clock_primary = clock[
    clock["sample_type_code"] == "01"
].copy()

print(
    "Primary CLOCK16 mutation events:",
    len(clock_primary)
)

print(
    "Primary CLOCK16 mutant patients:",
    clock_primary["patient_id"].nunique()
)

# ============================================================
# ATTACH CANCER TYPE
# ============================================================

clock_primary = clock_primary.merge(
    primary_patients,
    on="patient_id",
    how="inner"
)

print(
    "Cancer types:",
    clock_primary["cancer_type"].nunique()
)

# ============================================================
# UNIQUE PATIENT-GENE PAIRS
# ============================================================

patient_gene = (
    clock_primary[
        ["patient_id", "gene", "cancer_type"]
    ]
    .drop_duplicates()
)

# ============================================================
# CANCER TYPES
# ============================================================

cancer_types = sorted(
    primary_patients["cancer_type"]
    .dropna()
    .unique()
)

# ============================================================
# TESTS
# ============================================================

results = []

for cancer in cancer_types:

    cancer_patients = set(
        primary_patients.loc[
            primary_patients["cancer_type"] == cancer,
            "patient_id"
        ]
    )

    cancer_n = len(cancer_patients)

    for gene in clock_genes:

        gene_patients = set(
            patient_gene.loc[
                patient_gene["gene"] == gene,
                "patient_id"
            ]
        )

        # ----------------------------------------------------
        # 2x2 TABLE
        #
        #                Gene+   Gene-
        # Cancer           a       b
        # Other            c       d
        # ----------------------------------------------------

        a = len(cancer_patients & gene_patients)

        b = cancer_n - a

        other_patients = set(
            primary_patients.loc[
                primary_patients["cancer_type"] != cancer,
                "patient_id"
            ]
        )

        c = len(other_patients & gene_patients)

        d = len(other_patients) - c

        table = [
            [a, b],
            [c, d]
        ]

        odds_ratio, p_value = fisher_exact(
            table,
            alternative="two-sided"
        )

        cancer_prev = (
            100 * a / cancer_n
            if cancer_n > 0 else np.nan
        )

        other_n = len(other_patients)

        other_prev = (
            100 * c / other_n
            if other_n > 0 else np.nan
        )

        if odds_ratio > 1:
            direction = "Enriched"
        elif odds_ratio < 1:
            direction = "Depleted"
        else:
            direction = "Neutral"

        results.append({

            "Cancer_Type": cancer,

            "Gene": gene,

            "Cancer_Patients": cancer_n,

            "Gene_Mutated_Cancer_Patients": a,

            "Cancer_Prevalence_Percent": cancer_prev,

            "Other_Prevalence_Percent": other_prev,

            "Odds_Ratio": odds_ratio,

            "Fisher_P": p_value,

            "Direction": direction

        })

# ============================================================
# DATAFRAME
# ============================================================

results = pd.DataFrame(results)

# ============================================================
# FDR
# ============================================================

results["FDR"] = multipletests(
    results["Fisher_P"],
    method="fdr_bh"
)[1]

# ============================================================
# ORDER
# ============================================================

results = results.sort_values(
    ["FDR", "Fisher_P"]
).reset_index(drop=True)

# ============================================================
# SUMMARY
# ============================================================

print("\n========================================")
print("RESULT SUMMARY")
print("========================================")

print(
    "Total gene × cancer tests:",
    len(results)
)

print(
    "FDR < 0.05:",
    (results["FDR"] < 0.05).sum()
)

# ============================================================
# TOP SIGNIFICANT
# ============================================================

sig = results[
    results["FDR"] < 0.05
].copy()

print("\n=== TOP SIGNIFICANT ASSOCIATIONS ===")

if len(sig) == 0:

    print("No significant gene × cancer associations.")

else:

    print(
        sig[
            [
                "Cancer_Type",
                "Gene",
                "Cancer_Patients",
                "Gene_Mutated_Cancer_Patients",
                "Cancer_Prevalence_Percent",
                "Other_Prevalence_Percent",
                "Odds_Ratio",
                "Fisher_P",
                "FDR",
                "Direction"
            ]
        ]
        .head(50)
        .to_string(index=False)
    )

# ============================================================
# TOP ENRICHED
# ============================================================

print("\n=== TOP ENRICHED ===")

print(
    sig[
        sig["Direction"] == "Enriched"
    ][
        [
            "Cancer_Type",
            "Gene",
            "Cancer_Prevalence_Percent",
            "Other_Prevalence_Percent",
            "Odds_Ratio",
            "FDR"
        ]
    ]
    .sort_values(
        ["FDR", "Odds_Ratio"],
        ascending=[True, False]
    )
    .head(30)
    .to_string(index=False)
)

# ============================================================
# TOP DEPLETED
# ============================================================

print("\n=== TOP DEPLETED ===")

print(
    sig[
        sig["Direction"] == "Depleted"
    ][
        [
            "Cancer_Type",
            "Gene",
            "Cancer_Prevalence_Percent",
            "Other_Prevalence_Percent",
            "Odds_Ratio",
            "FDR"
        ]
    ]
    .sort_values(
        ["FDR", "Odds_Ratio"],
        ascending=[True, True]
    )
    .head(30)
    .to_string(index=False)
)

# ============================================================
# SAVE
# ============================================================

results.to_csv(
    out_file,
    sep="\t",
    index=False
)

print("\n========================================")
print("SAVED")
print("========================================")

print(out_file)

