#!/usr/bin/env python3

import os
import warnings
import numpy as np
import pandas as pd

from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT = "results/clock16/PAN_CANCER_patient_mutational_burden.tsv"

OUT = "results/clock16/PAN_CANCER_burden_adjusted_FIRTH_LR.tsv"

QC_OUT = "results/clock16/PAN_CANCER_burden_adjusted_FIRTH_LR_QC.tsv"


# ============================================================
# NUMERICAL SETTINGS
# ============================================================

EPS = 1e-8
MAXITER = 1000


# ============================================================
# FIRTH PENALIZED LOG-LIKELIHOOD
#
# l*(beta) = l(beta) + 0.5 log |X' W X|
#
# This is the Jeffreys-prior penalized likelihood used
# for Firth logistic regression.
# ============================================================

def firth_penalized_loglik(beta, X, y):

    beta = np.asarray(beta, dtype=float)

    eta = X @ beta
    p = expit(eta)

    # Stable Bernoulli log-likelihood
    loglik = np.sum(
        y * np.log(np.clip(p, EPS, 1 - EPS))
        +
        (1 - y) * np.log(np.clip(1 - p, EPS, 1 - EPS))
    )

    W = p * (1.0 - p)

    XWX = X.T @ (W[:, None] * X)

    # Numerical stabilization
    XWX = XWX + np.eye(XWX.shape[0]) * EPS

    sign, logdet = np.linalg.slogdet(XWX)

    if sign <= 0 or not np.isfinite(logdet):
        return -np.inf

    return loglik + 0.5 * logdet


# ============================================================
# NEGATIVE OBJECTIVE
# ============================================================

def objective(beta, X, y):

    value = firth_penalized_loglik(beta, X, y)

    if not np.isfinite(value):
        return 1e100

    return -value


# ============================================================
# FIT FIRTH LOGISTIC REGRESSION
# ============================================================

def fit_firth(X, y):

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)

    # Starting values from ordinary logistic regression
    # are not required; zero initialization is robust.
    beta0 = np.zeros(X.shape[1], dtype=float)

    result = minimize(
        objective,
        beta0,
        args=(X, y),
        method="BFGS",
        options={
            "maxiter": MAXITER,
            "gtol": 1e-8
        }
    )

    # If BFGS fails, try L-BFGS-B
    if not result.success or not np.all(np.isfinite(result.x)):

        result = minimize(
            objective,
            beta0,
            args=(X, y),
            method="L-BFGS-B",
            options={
                "maxiter": MAXITER,
                "ftol": 1e-12,
                "gtol": 1e-8
            }
        )

    beta = np.asarray(result.x, dtype=float)

    pll = firth_penalized_loglik(beta, X, y)

    return {
        "beta": beta,
        "pll": pll,
        "success": bool(result.success),
        "message": str(result.message),
        "iterations": getattr(result, "nit", np.nan)
    }


# ============================================================
# PROFILE / APPROXIMATE 95% CI
#
# For the manuscript, the primary significance test is LR.
# CI is retained from the local curvature of the Firth
# penalized likelihood.
# ============================================================

def firth_se(beta, X, y):

    p = expit(X @ beta)
    W = p * (1.0 - p)

    XWX = X.T @ (W[:, None] * X)
    XWX = XWX + np.eye(XWX.shape[0]) * EPS

    try:
        cov = np.linalg.inv(XWX)
        se = np.sqrt(np.maximum(np.diag(cov), 0))
        return se
    except Exception:
        return np.repeat(np.nan, len(beta))


# ============================================================
# LOAD DATA
# ============================================================

print("=" * 80)
print("CLOCK16 PAN-CANCER BURDEN-ADJUSTED ASSOCIATION")
print("FIRTH PENALIZED LIKELIHOOD-RATIO TEST")
print("=" * 80)

df = pd.read_csv(INPUT, sep="\t")

print("\nInput:", INPUT)
print("Rows:", len(df))
print("Cancer types:", df["cancer_type"].nunique())

required = [
    "patient_id",
    "cancer_type",
    "CLOCK16_mutated",
    "MC3_total_mutations",
    "MC3_CLOCK16_mutations"
]

missing = [x for x in required if x not in df.columns]

if missing:
    raise ValueError(f"Missing columns: {missing}")


# ============================================================
# NON-CLOCK16 MUTATIONAL BURDEN
# ============================================================

df["non_CLOCK16_mutations"] = (
    df["MC3_total_mutations"]
    -
    df["MC3_CLOCK16_mutations"]
)

if (df["non_CLOCK16_mutations"] < 0).any():

    raise ValueError(
        "Negative non-CLOCK16 mutation burden detected."
    )

df["log1p_non_CLOCK16"] = np.log1p(
    df["non_CLOCK16_mutations"]
)

df["CLOCK16_mutated"] = (
    df["CLOCK16_mutated"].astype(int)
)


# ============================================================
# GLOBAL CHECK
# ============================================================

print("\nGlobal CLOCK16 status:")

print(
    df["CLOCK16_mutated"]
    .value_counts()
    .sort_index()
    .rename({
        0: "CLOCK16_wild_type",
        1: "CLOCK16_mutated"
    })
)

print("\nGlobal burden:")

for status, label in [
    (1, "CLOCK16-mutated"),
    (0, "CLOCK16-wild-type")
]:

    x = df.loc[
        df["CLOCK16_mutated"] == status,
        "non_CLOCK16_mutations"
    ]

    print(
        f"{label}: "
        f"N={len(x)}, "
        f"median={x.median():.1f}, "
        f"mean={x.mean():.2f}"
    )


# ============================================================
# PER-CANCER FIRTH LR ANALYSIS
# ============================================================

results = []
qc = []

cancers = sorted(
    df["cancer_type"].dropna().unique()
)

for cancer in cancers:

    sub = df[
        df["cancer_type"] == cancer
    ].copy()

    N = len(sub)

    y = sub["CLOCK16_mutated"].values.astype(float)

    n_mut = int(y.sum())
    n_wt = int(N - n_mut)

    x = sub["log1p_non_CLOCK16"].values.astype(float)

    # --------------------------------------------------------
    # Basic requirement
    # --------------------------------------------------------

    if n_mut == 0 or n_wt == 0:

        qc.append({
            "Cancer_Type": cancer,
            "N": N,
            "CLOCK16_Mutated": n_mut,
            "CLOCK16_Wild_Type": n_wt,
            "Full_Converged": False,
            "Null_Converged": False,
            "LR_Valid": False,
            "Reason": "Only one outcome class"
        })

        continue


    # --------------------------------------------------------
    # Full model
    #
    # CLOCK16 status ~ non-CLOCK16 burden
    # --------------------------------------------------------

    X_full = np.column_stack([
        np.ones(N),
        x
    ])

    # --------------------------------------------------------
    # Null model
    #
    # CLOCK16 status ~ 1
    # --------------------------------------------------------

    X_null = np.ones((N, 1))


    # --------------------------------------------------------
    # FIT
    # --------------------------------------------------------

    full = fit_firth(
        X_full,
        y
    )

    null = fit_firth(
        X_null,
        y
    )


    # --------------------------------------------------------
    # LR TEST
    #
    # Difference in penalized likelihood
    # --------------------------------------------------------

    lr_stat = np.nan
    p_lr = np.nan

    if (
        np.isfinite(full["pll"])
        and
        np.isfinite(null["pll"])
    ):

        lr_stat = 2.0 * (
            full["pll"] -
            null["pll"]
        )

        # Numerical noise can produce tiny negatives
        lr_stat = max(0.0, lr_stat)

        p_lr = chi2.sf(
            lr_stat,
            df=1
        )


    # --------------------------------------------------------
    # EFFECT ESTIMATE
    # --------------------------------------------------------

    beta = full["beta"][1]

    se = firth_se(
        full["beta"],
        X_full,
        y
    )[1]

    OR = np.exp(beta)

    if np.isfinite(se):

        CI_low = np.exp(
            beta - 1.96 * se
        )

        CI_high = np.exp(
            beta + 1.96 * se
        )

    else:

        CI_low = np.nan
        CI_high = np.nan


    # --------------------------------------------------------
    # DIRECTION
    # --------------------------------------------------------

    direction = (
        "Positive"
        if beta > 0
        else "Negative"
    )


    # --------------------------------------------------------
    # SAVE RESULT
    # --------------------------------------------------------

    results.append({

        "Cancer_Type": cancer,

        "N": N,

        "CLOCK16_Mutated": n_mut,

        "CLOCK16_Wild_Type": n_wt,

        "OR_per_log1p_non_CLOCK16_mutation": OR,

        "CI95_Low": CI_low,

        "CI95_High": CI_high,

        "Beta": beta,

        "SE": se,

        "Firth_LR_Statistic": lr_stat,

        "Firth_LR_P": p_lr,

        "Direction": direction,

        "Full_Iterations": full["iterations"],

        "Null_Iterations": null["iterations"],

        "Median_non_CLOCK16_mutations":
            sub["non_CLOCK16_mutations"].median()

    })


    qc.append({

        "Cancer_Type": cancer,

        "N": N,

        "CLOCK16_Mutated": n_mut,

        "CLOCK16_Wild_Type": n_wt,

        "Full_Converged": full["success"],

        "Null_Converged": null["success"],

        "Full_PLL": full["pll"],

        "Null_PLL": null["pll"],

        "LR_Statistic": lr_stat,

        "LR_Valid": np.isfinite(p_lr),

        "Full_Message": full["message"],

        "Null_Message": null["message"]

    })


# ============================================================
# RESULTS DATAFRAME
# ============================================================

res = pd.DataFrame(results)

if len(res) == 0:

    raise RuntimeError(
        "No valid cancer-specific models were fitted."
    )


# ============================================================
# MULTIPLE TESTING
# ============================================================

valid = res["Firth_LR_P"].notna()

res["FDR"] = np.nan

if valid.sum() > 0:

    res.loc[valid, "FDR"] = multipletests(
        res.loc[valid, "Firth_LR_P"],
        method="fdr_bh"
    )[1]


# ============================================================
# SORT
# ============================================================

res = res.sort_values(
    ["FDR", "Firth_LR_P"],
    na_position="last"
).reset_index(drop=True)


# ============================================================
# SAVE
# ============================================================

res.to_csv(
    OUT,
    sep="\t",
    index=False
)

qc_df = pd.DataFrame(qc)

qc_df.to_csv(
    QC_OUT,
    sep="\t",
    index=False
)


# ============================================================
# PRINT RESULTS
# ============================================================

print("\nResults:")

print(
    res[
        [
            "Cancer_Type",
            "N",
            "CLOCK16_Mutated",
            "OR_per_log1p_non_CLOCK16_mutation",
            "CI95_Low",
            "CI95_High",
            "Firth_LR_Statistic",
            "Firth_LR_P",
            "FDR",
            "Direction"
        ]
    ]
    .to_string(index=False)
)


print("\nSignificant Firth LR associations:")
print(
    int(
        (res["FDR"] < 0.05)
        .sum()
    )
)

print("\nFDR < 0.05:")

sig = res[
    res["FDR"] < 0.05
]

if len(sig) > 0:

    print(
        sig[
            [
                "Cancer_Type",
                "N",
                "CLOCK16_Mutated",
                "OR_per_log1p_non_CLOCK16_mutation",
                "CI95_Low",
                "CI95_High",
                "Firth_LR_P",
                "FDR"
            ]
        ]
        .to_string(index=False)
    )

else:

    print("None")


# ============================================================
# FINAL QC
# ============================================================

print("\nQC:")
print(
    f"Models attempted: {len(cancers)}"
)

print(
    f"Models fitted: {len(res)}"
)

print(
    f"Valid LR tests: "
    f"{res['Firth_LR_P'].notna().sum()}"
)

print(
    f"FDR < 0.05: "
    f"{(res['FDR'] < 0.05).sum()}"
)

print("\nSaved:")
print(OUT)
print(QC_OUT)

print("\n" + "=" * 80)
print("FIRTH LR BURDEN-ADJUSTED ANALYSIS COMPLETED")
print("=" * 80)

