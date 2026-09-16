#!/usr/bin/env python3

# ============================================================
# CLOCK16 GENE-LEVEL POOLED FIRTH LOGISTIC REGRESSION
#
# For each of the 16 CLOCK genes:
#
#   GENE_MUTATED ~ log1p(non-CLOCK16 burden) + cancer type
#
# Firth logistic regression
# Likelihood-ratio test for mutational burden
# BH-FDR across 16 genes
#
# ============================================================

import os
import warnings
import numpy as np
import pandas as pd

from scipy.optimize import minimize
from scipy.stats import chi2


# ============================================================
# INPUTS
# ============================================================

BURDEN_FILE = (
    "results/clock16/"
    "PAN_CANCER_patient_mutational_burden.tsv"
)

MUT_FILE = (
    "results/clock16/"
    "REAL_primary_CLOCK16_mutations.tsv.gz"
)

OUT = (
    "results/clock16/"
    "CLOCK16_gene_level_pooled_FIRTH_LR.tsv"
)

QC_OUT = (
    "results/clock16/"
    "CLOCK16_gene_level_pooled_FIRTH_LR_QC.tsv"
)


# ============================================================
# CLOCK16 GENES
# ============================================================

CLOCK_GENES = [
    "ARNTL",
    "CLOCK",
    "NPAS2",
    "PER1",
    "PER2",
    "PER3",
    "CRY1",
    "CRY2",
    "NR1D1",
    "NR1D2",
    "RORA",
    "RORB",
    "RORC",
    "CSNK1D",
    "CSNK1E",
    "TIMELESS"
]


MAX_ITER = 5000
TOL = 1e-9


# ============================================================
# FDR
# ============================================================

def bh_fdr(pvalues):

    p = np.asarray(pvalues, dtype=float)

    n = len(p)

    order = np.argsort(p)

    ranked = p[order]

    q = np.empty(n, dtype=float)

    prev = 1.0

    for i in range(n - 1, -1, -1):

        rank = i + 1

        value = (
            ranked[i] * n / rank
        )

        prev = min(prev, value)

        q[i] = prev

    out = np.empty(n, dtype=float)

    out[order] = q

    return np.clip(out, 0, 1)


# ============================================================
# FIRTH FUNCTIONS
# ============================================================

def sigmoid(z):

    z = np.clip(z, -700, 700)

    return 1.0 / (
        1.0 + np.exp(-z)
    )


def firth_loglik(beta, X, y):

    eta = X @ beta

    p = sigmoid(eta)

    eps = 1e-12

    ll = np.sum(
        y * np.log(
            np.clip(p, eps, 1 - eps)
        )
        +
        (1 - y)
        * np.log(
            np.clip(1 - p, eps, 1 - eps)
        )
    )

    W = p * (1 - p)

    XW = X * W[:, None]

    I = X.T @ XW

    sign, logdet = np.linalg.slogdet(I)

    if sign <= 0 or not np.isfinite(logdet):

        return -np.inf

    return (
        ll
        + 0.5 * logdet
    )


def objective(beta, X, y):

    value = firth_loglik(
        beta,
        X,
        y
    )

    if not np.isfinite(value):

        return 1e100

    return -value


def fit_firth(X, y):

    p = X.shape[1]

    start = np.zeros(p)

    result = minimize(
        objective,
        start,
        args=(X, y),
        method="BFGS",
        options={
            "maxiter": MAX_ITER,
            "gtol": TOL
        }
    )

    if (
        not result.success
        or not np.all(np.isfinite(result.x))
    ):

        result = minimize(
            objective,
            result.x,
            args=(X, y),
            method="L-BFGS-B",
            options={
                "maxiter": MAX_ITER,
                "ftol": 1e-12,
                "gtol": 1e-9
            }
        )

    beta = result.x

    ll = firth_loglik(
        beta,
        X,
        y
    )

    return {
        "beta": beta,
        "loglik": ll,
        "success": result.success,
        "message": str(result.message),
        "iterations": result.nit
    }


# ============================================================
# LOAD BURDEN DATA
# ============================================================

print("=" * 80)
print("CLOCK16 GENE-LEVEL POOLED FIRTH LR")
print("=" * 80)

print()

print("Loading burden data:")

df = pd.read_csv(
    BURDEN_FILE,
    sep="\t"
)

print(
    "Rows:",
    len(df)
)

print(
    "Unique patients:",
    df["patient_id"].nunique()
)


# ============================================================
# LOAD MUTATION DATA
# ============================================================

print()

print("Loading CLOCK16 mutation data:")

mut = pd.read_csv(
    MUT_FILE,
    sep="\t",
    compression="gzip"
)

print(
    "Rows:",
    len(mut)
)

print(
    "Columns:",
    list(mut.columns)
)


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

for col in [
    "patient_id",
    "gene"
]:

    if col not in mut.columns:

        raise ValueError(
            f"Missing mutation column: {col}"
        )


# ============================================================
# CHECK BURDEN DATA
# ============================================================

required_burden = [
    "patient_id",
    "cancer_type",
    "MC3_total_mutations",
    "MC3_CLOCK16_mutations"
]

for col in required_burden:

    if col not in df.columns:

        raise ValueError(
            f"Missing burden column: {col}"
        )


if df["patient_id"].duplicated().any():

    raise ValueError(
        "Duplicate patient IDs in burden data."
    )


# ============================================================
# NON-CLOCK16 BURDEN
# ============================================================

df["non_CLOCK16_mutations"] = (
    df["MC3_total_mutations"]
    -
    df["MC3_CLOCK16_mutations"]
)

if (
    df["non_CLOCK16_mutations"] < 0
).any():

    raise ValueError(
        "Negative non-CLOCK16 mutation burden."
    )


df["log1p_non_CLOCK16_burden"] = np.log1p(
    df["non_CLOCK16_mutations"]
)


# ============================================================
# MUTATION DATA
# ============================================================

mut["gene"] = (
    mut["gene"]
    .astype(str)
    .str.strip()
)

mut["patient_id"] = (
    mut["patient_id"]
    .astype(str)
    .str.strip()
)


# Keep only the 16 CLOCK genes

mut = mut[
    mut["gene"].isin(CLOCK_GENES)
].copy()


print()

print(
    "CLOCK16 mutation rows:",
    len(mut)
)

print(
    "CLOCK16 genes observed:",
    mut["gene"].nunique()
)

print()

print(
    mut["gene"]
    .value_counts()
    .reindex(CLOCK_GENES)
    .fillna(0)
    .astype(int)
    .to_string()
)


# ============================================================
# PATIENT × GENE MUTATION MATRIX
# ============================================================

gene_patient = (
    mut[
        ["patient_id", "gene"]
    ]
    .drop_duplicates()
)

gene_patient["mutated"] = 1

gene_matrix = (
    gene_patient
    .pivot_table(
        index="patient_id",
        columns="gene",
        values="mutated",
        aggfunc="max",
        fill_value=0
    )
)


for gene in CLOCK_GENES:

    if gene not in gene_matrix.columns:

        gene_matrix[gene] = 0


gene_matrix = gene_matrix[
    CLOCK_GENES
]


# ============================================================
# MERGE WITH PAN-CANCER DENOMINATOR
# ============================================================

analysis = df[
    [
        "patient_id",
        "cancer_type",
        "log1p_non_CLOCK16_burden"
    ]
].copy()


analysis = analysis.merge(
    gene_matrix,
    left_on="patient_id",
    right_index=True,
    how="left"
)


analysis[CLOCK_GENES] = (
    analysis[CLOCK_GENES]
    .fillna(0)
    .astype(int)
)


print()

print(
    "Analysis patients:",
    len(analysis)
)

print(
    "Analysis unique patients:",
    analysis["patient_id"].nunique()
)


# ============================================================
# CANCER TYPE DESIGN
# ============================================================

cancer_levels = sorted(
    analysis["cancer_type"].unique()
)

reference_cancer = cancer_levels[0]

print()

print(
    "Cancer reference:",
    reference_cancer
)

print(
    "Cancer types:",
    len(cancer_levels)
)


cancer_cat = pd.Categorical(
    analysis["cancer_type"],
    categories=cancer_levels
)

cancer_dummies = pd.get_dummies(
    cancer_cat,
    drop_first=True,
    dtype=float
)


# ============================================================
# BURDEN STANDARDIZATION
# ============================================================

burden = analysis[
    "log1p_non_CLOCK16_burden"
].to_numpy(dtype=float)

burden_mean = burden.mean()

burden_sd = burden.std(ddof=0)

if burden_sd == 0:

    raise ValueError(
        "Burden has zero variance."
    )

burden_z = (
    burden - burden_mean
) / burden_sd


# ============================================================
# BASE DESIGN MATRIX
#
# Intercept + burden + cancer type
# ============================================================

X_base_df = pd.DataFrame({

    "Intercept":
        np.ones(len(analysis)),

    "log1p_non_CLOCK16_burden_z":
        burden_z
})


for col in cancer_dummies.columns:

    X_base_df[col] = (
        cancer_dummies[col]
        .to_numpy()
    )


X_base = X_base_df.to_numpy(
    dtype=float
)


# ============================================================
# REDUCED MODEL
#
# Gene ~ cancer type
# ============================================================

X_reduced_df = pd.DataFrame({

    "Intercept":
        np.ones(len(analysis))
})


for col in cancer_dummies.columns:

    X_reduced_df[col] = (
        cancer_dummies[col]
        .to_numpy()
    )


X_reduced = X_reduced_df.to_numpy(
    dtype=float
)


print()

print(
    "Reduced model parameters:",
    X_reduced.shape[1]
)

print(
    "Base/full model parameters:",
    X_base.shape[1]
)


# ============================================================
# ANALYZE EACH GENE
# ============================================================

results = []

qc_rows = []


for gene in CLOCK_GENES:

    print()
    print("-" * 80)
    print("GENE:", gene)
    print("-" * 80)

    y = analysis[
        gene
    ].to_numpy(dtype=float)

    n_mut = int(
        y.sum()
    )

    n_wt = int(
        len(y) - n_mut
    )

    print(
        "Mutated:",
        n_mut,
        "WT:",
        n_wt
    )


    # --------------------------------------------------------
    # Minimum information check
    # --------------------------------------------------------

    if n_mut == 0:

        print(
            "SKIPPED: no mutated patients."
        )

        qc_rows.append({

            "Gene": gene,
            "N": len(y),
            "Mutated": n_mut,
            "Wild_Type": n_wt,
            "Status": "Skipped_no_mutants"

        })

        continue


    if n_wt == 0:

        print(
            "SKIPPED: no wild-type patients."
        )

        qc_rows.append({

            "Gene": gene,
            "N": len(y),
            "Mutated": n_mut,
            "Wild_Type": n_wt,
            "Status": "Skipped_no_wild_type"

        })

        continue


    # --------------------------------------------------------
    # FULL MODEL
    #
    # gene ~ burden + cancer type
    # --------------------------------------------------------

    X_full = X_base.copy()

    # Gene-specific reduced/full model is identical
    # in covariates, with gene as outcome.

    full = fit_firth(
        X_full,
        y
    )


    # --------------------------------------------------------
    # REDUCED MODEL
    #
    # gene ~ cancer type
    # --------------------------------------------------------

    reduced = fit_firth(
        X_reduced,
        y
    )


    # --------------------------------------------------------
    # LR TEST
    # --------------------------------------------------------

    LR = 2.0 * (
        full["loglik"]
        -
        reduced["loglik"]
    )

    if LR < 0:

        LR = 0.0


    p_lr = chi2.sf(
        LR,
        df=1
    )


    # --------------------------------------------------------
    # BURDEN EFFECT
    # --------------------------------------------------------

    beta_z = full["beta"][1]

    beta_original = (
        beta_z / burden_sd
    )

    OR = np.exp(
        np.clip(
            beta_original,
            -700,
            700
        )
    )


    # --------------------------------------------------------
    # CI
    # --------------------------------------------------------

    eta = X_full @ full["beta"]

    p_hat = sigmoid(eta)

    W = p_hat * (
        1 - p_hat
    )

    XW = X_full * W[:, None]

    I = X_full.T @ XW

    cov = np.linalg.pinv(I)

    variance_z = max(
        cov[1, 1],
        0
    )

    se_z = np.sqrt(
        variance_z
    )

    se_original = (
        se_z / burden_sd
    )


    CI_low = np.exp(
        np.clip(
            beta_original
            -
            1.96 * se_original,
            -700,
            700
        )
    )

    CI_high = np.exp(
        np.clip(
            beta_original
            +
            1.96 * se_original,
            -700,
            700
        )
    )


    # --------------------------------------------------------
    # BURDEN MEDIANS
    # --------------------------------------------------------

    gene_mut_burden = analysis.loc[
        analysis[gene] == 1,
        "log1p_non_CLOCK16_burden"
    ]

    gene_wt_burden = analysis.loc[
        analysis[gene] == 0,
        "log1p_non_CLOCK16_burden"
    ]


    # --------------------------------------------------------
    # RESULT
    # --------------------------------------------------------

    results.append({

        "Gene":
            gene,

        "N":
            len(y),

        "Mutated_Patients":
            n_mut,

        "Wild_Type_Patients":
            n_wt,

        "OR_per_log1p_non_CLOCK16_mutation":
            OR,

        "CI95_Low":
            CI_low,

        "CI95_High":
            CI_high,

        "Beta":
            beta_original,

        "SE":
            se_original,

        "Firth_LR_Statistic":
            LR,

        "Firth_LR_P":
            p_lr,

        "Direction":
            (
                "Positive"
                if beta_original > 0
                else "Negative"
            ),

        "Median_log1p_Burden_Mutated":
            gene_mut_burden.median(),

        "Median_log1p_Burden_WT":
            gene_wt_burden.median(),

        "Full_Model_Fitted":
            full["success"],

        "Reduced_Model_Fitted":
            reduced["success"],

        "Full_Iterations":
            full["iterations"],

        "Reduced_Iterations":
            reduced["iterations"]

    })


    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------

    qc_rows.append({

        "Gene":
            gene,

        "N":
            len(y),

        "Mutated":
            n_mut,

        "Wild_Type":
            n_wt,

        "Full_Model_Fitted":
            full["success"],

        "Reduced_Model_Fitted":
            reduced["success"],

        "Full_LogLik":
            full["loglik"],

        "Reduced_LogLik":
            reduced["loglik"],

        "LR":
            LR,

        "LR_P":
            p_lr,

        "Status":
            "Completed"

    })


# ============================================================
# RESULTS DATAFRAME
# ============================================================

results_df = pd.DataFrame(
    results
)


# ============================================================
# BH-FDR
# ============================================================

if len(results_df) > 0:

    results_df["FDR"] = bh_fdr(
        results_df["Firth_LR_P"]
        .to_numpy()
    )

    results_df = results_df.sort_values(
        "Firth_LR_P"
    ).reset_index(drop=True)


# ============================================================
# SIGNIFICANCE
# ============================================================

if len(results_df) > 0:

    results_df["Significant_FDR_0.05"] = (
        results_df["FDR"] < 0.05
    )


# ============================================================
# SAVE
# ============================================================

results_df.to_csv(
    OUT,
    sep="\t",
    index=False
)

qc_df = pd.DataFrame(
    qc_rows
)

qc_df.to_csv(
    QC_OUT,
    sep="\t",
    index=False
)


# ============================================================
# PRINT
# ============================================================

print()
print("=" * 80)
print("16-GENE POOLED FIRTH LR RESULTS")
print("=" * 80)

print()

if len(results_df) > 0:

    print(
        results_df[
            [
                "Gene",
                "N",
                "Mutated_Patients",
                "OR_per_log1p_non_CLOCK16_mutation",
                "CI95_Low",
                "CI95_High",
                "Firth_LR_P",
                "FDR",
                "Direction"
            ]
        ].to_string(index=False)
    )


print()

print(
    "Models attempted:",
    len(CLOCK_GENES)
)

print(
    "Models completed:",
    len(results_df)
)

if len(results_df) > 0:

    print(
        "FDR < 0.05:",
        int(
            (
                results_df["FDR"] < 0.05
            ).sum()
        )
    )


print()

print("Saved:")
print(OUT)
print(QC_OUT)

print()

print("=" * 80)
print("CLOCK16 GENE-LEVEL POOLED FIRTH LR COMPLETED")
print("=" * 80)
