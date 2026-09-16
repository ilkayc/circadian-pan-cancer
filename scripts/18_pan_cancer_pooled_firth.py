#!/usr/bin/env python3

# ============================================================
# CLOCK16 POOLED PAN-CANCER FIRTH LOGISTIC REGRESSION
#
# Model:
#   CLOCK16_mutated ~ log1p(non-CLOCK16 mutation burden)
#                     + cancer type
#
# Purpose:
#   Test whether CLOCK16 mutation status remains associated
#   with mutational burden after adjustment for cancer type.
#
# Firth logistic regression:
#   Jeffreys-prior penalized likelihood
#
# Likelihood-ratio test:
#   Full model vs reduced model
#
# ============================================================

import os
import warnings
import numpy as np
import pandas as pd

from scipy.optimize import minimize
from scipy.stats import chi2


# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT = "results/clock16/PAN_CANCER_patient_mutational_burden.tsv"

OUT = "results/clock16/PAN_CANCER_pooled_FIRTH_LR.tsv"

QC_OUT = "results/clock16/PAN_CANCER_pooled_FIRTH_LR_QC.tsv"


# ============================================================
# SETTINGS
# ============================================================

MAX_ITER = 5000
TOL = 1e-9


# ============================================================
# LOAD DATA
# ============================================================

print("=" * 80)
print("CLOCK16 POOLED PAN-CANCER FIRTH LOGISTIC REGRESSION")
print("=" * 80)

print()

df = pd.read_csv(INPUT, sep="\t")

print("Input:", INPUT)
print("Rows:", len(df))
print("Columns:", list(df.columns))

required = [
    "patient_id",
    "cancer_type",
    "CLOCK16_mutated",
    "MC3_total_mutations",
    "MC3_CLOCK16_mutations"
]

missing = [x for x in required if x not in df.columns]

if missing:
    raise ValueError(
        "Missing required columns: " + ", ".join(missing)
    )


# ============================================================
# DATA INTEGRITY
# ============================================================

if df["patient_id"].duplicated().any():
    raise ValueError("Duplicate patient_id detected.")

if df["cancer_type"].isna().any():
    raise ValueError("Missing cancer_type detected.")

if df["CLOCK16_mutated"].isna().any():
    raise ValueError("Missing CLOCK16_mutated detected.")

if df["MC3_total_mutations"].isna().any():
    raise ValueError("Missing MC3_total_mutations detected.")

if df["MC3_CLOCK16_mutations"].isna().any():
    raise ValueError("Missing MC3_CLOCK16_mutations detected.")


# ============================================================
# DEFINE NON-CLOCK16 BURDEN
# ============================================================

df["non_CLOCK16_mutations"] = (
    df["MC3_total_mutations"]
    - df["MC3_CLOCK16_mutations"]
)

if (df["non_CLOCK16_mutations"] < 0).any():
    raise ValueError(
        "Negative non-CLOCK16 mutation burden detected."
    )


df["log1p_non_CLOCK16_burden"] = np.log1p(
    df["non_CLOCK16_mutations"]
)


# ============================================================
# BINARY OUTCOME
# ============================================================

df["CLOCK16_status"] = (
    df["CLOCK16_mutated"]
    .astype(int)
)


# ============================================================
# BASIC SUMMARY
# ============================================================

print()
print("Global CLOCK16 status:")
print(
    df["CLOCK16_mutated"]
    .value_counts()
    .rename(
        index={
            False: "CLOCK16_wild_type",
            True: "CLOCK16_mutated"
        }
    )
)

print()

print("Cancer types:", df["cancer_type"].nunique())

print()

print("Non-CLOCK16 burden:")
print(
    df.groupby("CLOCK16_status")["non_CLOCK16_mutations"]
    .agg(
        N="size",
        Median="median",
        Mean="mean",
        Q1=lambda x: x.quantile(0.25),
        Q3=lambda x: x.quantile(0.75),
        Max="max"
    )
)


# ============================================================
# DESIGN MATRIX
# ============================================================

# Reference cancer type:
# first alphabetically after sorting

cancer_levels = sorted(
    df["cancer_type"].unique()
)

reference_cancer = cancer_levels[0]

print()
print("Cancer reference level:")
print(reference_cancer)

print()
print("Number of cancer types:", len(cancer_levels))


# Create dummy variables manually
cancer_cat = pd.Categorical(
    df["cancer_type"],
    categories=cancer_levels
)

cancer_dummies = pd.get_dummies(
    cancer_cat,
    drop_first=True,
    dtype=float
)


# ============================================================
# MODEL VARIABLES
# ============================================================

y = df["CLOCK16_status"].to_numpy(dtype=float)

burden = df[
    "log1p_non_CLOCK16_burden"
].to_numpy(dtype=float)


# Standardize burden for numerical stability ONLY.
#
# Interpretation will later be converted back to
# OR per 1-unit increase in log1p burden.

burden_mean = burden.mean()
burden_sd = burden.std(ddof=0)

if burden_sd == 0:
    raise ValueError("Mutation burden has zero variance.")

burden_z = (
    (burden - burden_mean)
    / burden_sd
)


# ============================================================
# FULL MODEL
#
# CLOCK16 ~ burden + cancer type
# ============================================================

X_full_df = pd.DataFrame({
    "Intercept": np.ones(len(df)),
    "log1p_non_CLOCK16_burden_z": burden_z
})

for col in cancer_dummies.columns:
    X_full_df[col] = cancer_dummies[col].to_numpy()

X_full = X_full_df.to_numpy(dtype=float)


# ============================================================
# REDUCED MODEL
#
# CLOCK16 ~ cancer type
#
# Used for LR test of burden.
# ============================================================

X_reduced_df = pd.DataFrame({
    "Intercept": np.ones(len(df))
})

for col in cancer_dummies.columns:
    X_reduced_df[col] = cancer_dummies[col].to_numpy()

X_reduced = X_reduced_df.to_numpy(dtype=float)


# ============================================================
# CHECK MATRIX RANK
# ============================================================

rank_full = np.linalg.matrix_rank(X_full)
rank_reduced = np.linalg.matrix_rank(X_reduced)

print()
print("Full model:")
print("  N =", X_full.shape[0])
print("  Parameters =", X_full.shape[1])
print("  Rank =", rank_full)

print()
print("Reduced model:")
print("  N =", X_reduced.shape[0])
print("  Parameters =", X_reduced.shape[1])
print("  Rank =", rank_reduced)


if rank_full != X_full.shape[1]:
    raise ValueError(
        "Full design matrix is rank deficient."
    )

if rank_reduced != X_reduced.shape[1]:
    raise ValueError(
        "Reduced design matrix is rank deficient."
    )


# ============================================================
# FIRTH LOGISTIC REGRESSION
#
# Penalized log-likelihood:
#
#   l(beta) + 0.5 log |I(beta)|
#
# where:
#
#   I(beta) = X' W X
#
# ============================================================

def sigmoid(z):

    z = np.clip(z, -700, 700)

    return 1.0 / (
        1.0 + np.exp(-z)
    )


def firth_loglik(beta, X, y):

    eta = X @ beta

    p = sigmoid(eta)

    # ordinary binomial log likelihood
    eps = 1e-12

    ll = np.sum(
        y * np.log(np.clip(p, eps, 1 - eps))
        +
        (1 - y)
        * np.log(np.clip(1 - p, eps, 1 - eps))
    )

    W = p * (1 - p)

    # Fisher information
    XW = X * W[:, None]

    I = X.T @ XW

    sign, logdet = np.linalg.slogdet(I)

    if sign <= 0 or not np.isfinite(logdet):
        return -np.inf

    return ll + 0.5 * logdet


def objective(beta, X, y):

    val = firth_loglik(beta, X, y)

    if not np.isfinite(val):
        return 1e100

    return -val


def fit_firth(X, y, start=None):

    p = X.shape[1]

    if start is None:
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

    # If BFGS reports precision issues but objective is finite,
    # perform a second optimization with L-BFGS-B.

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
        "message": result.message,
        "iterations": result.nit
    }


# ============================================================
# FIT REDUCED MODEL
# ============================================================

print()
print("Fitting reduced Firth model...")

reduced = fit_firth(
    X_reduced,
    y
)

print(
    "Reduced model success:",
    reduced["success"]
)

print(
    "Reduced model log-likelihood:",
    reduced["loglik"]
)


# ============================================================
# FIT FULL MODEL
# ============================================================

print()
print("Fitting full Firth model...")

full = fit_firth(
    X_full,
    y
)

print(
    "Full model success:",
    full["success"]
)

print(
    "Full model log-likelihood:",
    full["loglik"]
)


if not full["success"]:
    warnings.warn(
        "Full model optimizer did not report success."
    )


# ============================================================
# BURDEN EFFECT
# ============================================================

beta_z = full["beta"][1]


# Convert standardized beta back to original scale:
#
# beta_original = beta_z / SD

beta_original = (
    beta_z / burden_sd
)

OR = np.exp(beta_original)


# ============================================================
# APPROXIMATE COVARIANCE MATRIX
#
# Based on inverse Fisher information at fitted beta.
# ============================================================

eta = X_full @ full["beta"]

p_hat = sigmoid(eta)

W = p_hat * (1 - p_hat)

XW = X_full * W[:, None]

I = X_full.T @ XW

cov = np.linalg.pinv(I)

se_z = np.sqrt(
    max(cov[1, 1], 0)
)

se_original = (
    se_z / burden_sd
)


CI_low = np.exp(
    beta_original
    - 1.96 * se_original
)

CI_high = np.exp(
    beta_original
    + 1.96 * se_original
)


# ============================================================
# LIKELIHOOD-RATIO TEST
#
# H0:
#   burden coefficient = 0
#
# Full:
#   burden + cancer type
#
# Reduced:
#   cancer type only
#
# LR = 2 * (LL_full - LL_reduced)
#
# df = 1
# ============================================================

LR = 2.0 * (
    full["loglik"]
    - reduced["loglik"]
)

if LR < 0:
    LR = 0.0

LR_p = chi2.sf(
    LR,
    df=1
)


# ============================================================
# GLOBAL BURDEN DIFFERENCE
# ============================================================

from scipy.stats import mannwhitneyu

mutated_burden = df.loc[
    df["CLOCK16_mutated"],
    "non_CLOCK16_mutations"
]

wt_burden = df.loc[
    ~df["CLOCK16_mutated"],
    "non_CLOCK16_mutations"
]

mw = mannwhitneyu(
    mutated_burden,
    wt_burden,
    alternative="two-sided"
)


# ============================================================
# RESULT TABLE
# ============================================================

result = pd.DataFrame([{

    "Analysis":
        "Pooled pan-cancer Firth LR",

    "N":
        len(df),

    "Cancer_Types":
        df["cancer_type"].nunique(),

    "CLOCK16_Mutated":
        int(df["CLOCK16_mutated"].sum()),

    "CLOCK16_Wild_Type":
        int((~df["CLOCK16_mutated"]).sum()),

    "Reference_Cancer":
        reference_cancer,

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
        LR_p,

    "Direction":
        "Positive"
        if beta_original > 0
        else "Negative",

    "Reduced_LogLik":
        reduced["loglik"],

    "Full_LogLik":
        full["loglik"],

    "Full_Model_Iterations":
        full["iterations"],

    "Reduced_Model_Iterations":
        reduced["iterations"],

    "Global_MW_P":
        mw.pvalue,

    "Global_Mutated_Median_Burden":
        mutated_burden.median(),

    "Global_WT_Median_Burden":
        wt_burden.median()

}])


# ============================================================
# QC
# ============================================================

qc = pd.DataFrame({

    "Metric": [

        "Input rows",

        "Unique patients",

        "Cancer types",

        "Full model parameters",

        "Reduced model parameters",

        "Full model rank",

        "Reduced model rank",

        "Full model fitted",

        "Reduced model fitted",

        "Full log-likelihood",

        "Reduced log-likelihood",

        "LR statistic",

        "LR degrees_of_freedom",

        "LR P-value",

        "OR",

        "CI low",

        "CI high",

        "Mutation burden variable",

        "Burden standardized before fitting",

        "No negative non-CLOCK16 burden"

    ],

    "Value": [

        len(df),

        df["patient_id"].nunique(),

        df["cancer_type"].nunique(),

        X_full.shape[1],

        X_reduced.shape[1],

        rank_full,

        rank_reduced,

        full["success"],

        reduced["success"],

        full["loglik"],

        reduced["loglik"],

        LR,

        1,

        LR_p,

        OR,

        CI_low,

        CI_high,

        "log1p(non_CLOCK16_mutations)",

        True,

        bool(
            (df["non_CLOCK16_mutations"] >= 0).all()
        )

    ]

})


# ============================================================
# SAVE
# ============================================================

result.to_csv(
    OUT,
    sep="\t",
    index=False
)

qc.to_csv(
    QC_OUT,
    sep="\t",
    index=False
)


# ============================================================
# PRINT FINAL RESULT
# ============================================================

print()
print("=" * 80)
print("POOLED PAN-CANCER FIRTH LR RESULT")
print("=" * 80)

print()

print(
    result[
        [
            "N",
            "Cancer_Types",
            "CLOCK16_Mutated",
            "CLOCK16_Wild_Type",
            "OR_per_log1p_non_CLOCK16_mutation",
            "CI95_Low",
            "CI95_High",
            "Firth_LR_Statistic",
            "Firth_LR_P",
            "Direction"
        ]
    ].to_string(index=False)
)

print()

print("Interpretation:")
print(
    "OR > 1 indicates that higher non-CLOCK16 "
    "mutational burden is associated with increased "
    "odds of CLOCK16 mutation after adjustment "
    "for cancer type."
)

print()

print("Saved:")
print(OUT)
print(QC_OUT)

print()
print("=" * 80)
print("POOLED PAN-CANCER FIRTH LR ANALYSIS COMPLETED")
print("=" * 80)
