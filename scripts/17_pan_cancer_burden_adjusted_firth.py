#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

# ============================================================
# CLOCK16 PAN-CANCER BURDEN-ADJUSTED ASSOCIATION
# FIRTH LOGISTIC REGRESSION
#
# No external firth package required.
#
# Outcome:
#   CLOCK16_mutated
#
# Predictor:
#   log1p(non-CLOCK16 mutation burden)
#
# Analysis:
#   Cancer-type specific
#
# ============================================================

INPUT = "results/clock16/PAN_CANCER_patient_mutational_burden.tsv"

OUT = "results/clock16/PAN_CANCER_burden_adjusted_FIRTH.tsv"
QC_OUT = "results/clock16/PAN_CANCER_burden_adjusted_FIRTH_QC.tsv"


print("=" * 80)
print("CLOCK16 PAN-CANCER BURDEN-ADJUSTED ASSOCIATION")
print("FIRTH LOGISTIC REGRESSION")
print("=" * 80)

# ============================================================
# 1. LOAD DATA
# ============================================================

df = pd.read_csv(INPUT, sep="\t")

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

print("\nInput:", INPUT)
print("Rows:", len(df))
print("Cancer types:", df["cancer_type"].nunique())

# ============================================================
# 2. CONSTRUCT NON-CLOCK16 BURDEN
# ============================================================

df["non_CLOCK16_mutations"] = (
    df["MC3_total_mutations"] -
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
# 3. FIRTH LOGISTIC REGRESSION
# ============================================================
#
# Penalized log-likelihood:
#
# l_F = l + 0.5 * log |I(beta)|
#
# where I(beta) = X'WX
#
# Optimization is performed by Newton-Raphson.
#
# ============================================================

def firth_logistic(X, y, max_iter=200, tol=1e-8):

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)

    n, p = X.shape

    beta = np.zeros(p)

    def penalized_loglik(beta):

        eta = X @ beta

        # stable logistic log-likelihood
        loglik = np.sum(
            y * eta -
            np.logaddexp(0, eta)
        )

        mu = expit(eta)

        W = mu * (1.0 - mu)

        # prevent singular matrix
        W = np.maximum(W, 1e-12)

        XWX = X.T @ (W[:, None] * X)

        sign, logdet = np.linalg.slogdet(XWX)

        if sign <= 0 or not np.isfinite(logdet):
            return -np.inf

        return loglik + 0.5 * logdet

    for iteration in range(max_iter):

        eta = X @ beta
        mu = expit(eta)

        W = mu * (1.0 - mu)
        W = np.maximum(W, 1e-12)

        XWX = X.T @ (W[:, None] * X)

        try:
            XWX_inv = np.linalg.inv(XWX)
        except np.linalg.LinAlgError:
            XWX_inv = np.linalg.pinv(XWX)

        # Hat diagonal
        WX = W[:, None] * X

        H_diag = np.sum(
            (X @ XWX_inv) * X,
            axis=1
        ) * W

        # Firth adjusted score
        adjusted_y = y - mu + H_diag * (0.5 - mu)

        score = X.T @ adjusted_y

        try:
            step = np.linalg.solve(XWX, score)
        except np.linalg.LinAlgError:
            step = XWX_inv @ score

        # ----------------------------------------------------
        # Step-halving
        # ----------------------------------------------------

        old_ll = penalized_loglik(beta)

        step_factor = 1.0
        accepted = False

        for _ in range(30):

            beta_new = beta + step_factor * step

            new_ll = penalized_loglik(beta_new)

            if np.isfinite(new_ll) and new_ll >= old_ll:
                accepted = True
                break

            step_factor *= 0.5

        if not accepted:
            break

        delta = np.max(np.abs(beta_new - beta))

        beta = beta_new

        if delta < tol:
            break

    # final covariance
    eta = X @ beta
    mu = expit(eta)

    W = mu * (1.0 - mu)
    W = np.maximum(W, 1e-12)

    XWX = X.T @ (W[:, None] * X)

    try:
        cov = np.linalg.inv(XWX)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(XWX)

    se = np.sqrt(
        np.maximum(np.diag(cov), 0)
    )

    return beta, se, iteration + 1


# ============================================================
# 4. GLOBAL BURDEN CHECK
# ============================================================

mut = df.loc[
    df["CLOCK16_mutated"] == 1,
    "non_CLOCK16_mutations"
]

wt = df.loc[
    df["CLOCK16_mutated"] == 0,
    "non_CLOCK16_mutations"
]

from scipy.stats import mannwhitneyu

mw_p = mannwhitneyu(
    mut,
    wt,
    alternative="two-sided"
).pvalue

print("\nGlobal burden:")
print("CLOCK16-mutated:", len(mut))
print("CLOCK16-wild-type:", len(wt))
print("Mutated median:", np.median(mut))
print("Wild-type median:", np.median(wt))
print("Mann-Whitney P:", mw_p)

# ============================================================
# 5. CANCER-SPECIFIC FIRTH MODELS
# ============================================================

results = []
qc = []

for cancer in sorted(df["cancer_type"].dropna().unique()):

    sub = df[
        df["cancer_type"] == cancer
    ].copy()

    n = len(sub)

    y = sub["CLOCK16_mutated"].astype(int).values

    n_mut = int(y.sum())
    n_wt = int(n - n_mut)

    x = sub["log1p_non_CLOCK16"].astype(float).values

    # Need both outcome classes
    if n_mut == 0 or n_wt == 0:

        qc.append({
            "Cancer_Type": cancer,
            "N": n,
            "CLOCK16_Mutated": n_mut,
            "CLOCK16_Wild_Type": n_wt,
            "Status": "SKIPPED_SINGLE_OUTCOME"
        })

        continue

    # --------------------------------------------------------
    # Model
    # intercept + burden
    # --------------------------------------------------------

    X = np.column_stack([
        np.ones(n),
        x
    ])

    try:

        beta, se, iterations = firth_logistic(
            X,
            y
        )

        b = beta[1]
        s = se[1]

        OR = np.exp(b)

        CI_low = np.exp(
            b - 1.96 * s
        )

        CI_high = np.exp(
            b + 1.96 * s
        )

        # Wald statistic for burden coefficient
        if s > 0 and np.isfinite(s):

            z = b / s

            p = 2.0 * (
                1.0 -
                __import__("scipy").stats.norm.cdf(
                    abs(z)
                )
            )

        else:

            p = np.nan

        direction = (
            "Positive"
            if b > 0
            else "Negative"
        )

        results.append({

            "Cancer_Type": cancer,
            "N": n,
            "CLOCK16_Mutated": n_mut,
            "CLOCK16_Wild_Type": n_wt,

            "OR_per_log1p_non_CLOCK16_mutation":
                OR,

            "CI95_Low":
                CI_low,

            "CI95_High":
                CI_high,

            "Beta":
                b,

            "SE":
                s,

            "P":
                p,

            "Direction":
                direction,

            "Iterations":
                iterations,

            "Median_non_CLOCK16_mutations":
                np.median(
                    sub["non_CLOCK16_mutations"]
                )

        })

        qc.append({
            "Cancer_Type": cancer,
            "N": n,
            "CLOCK16_Mutated": n_mut,
            "CLOCK16_Wild_Type": n_wt,
            "Status": "OK",
            "Iterations": iterations
        })

    except Exception as e:

        qc.append({
            "Cancer_Type": cancer,
            "N": n,
            "CLOCK16_Mutated": n_mut,
            "CLOCK16_Wild_Type": n_wt,
            "Status": "ERROR",
            "Error": str(e)
        })

# ============================================================
# 6. MULTIPLE TESTING
# ============================================================

res = pd.DataFrame(results)

if len(res) > 0:

    valid = res["P"].notna()

    res["FDR"] = np.nan

    if valid.sum() > 0:

        res.loc[valid, "FDR"] = multipletests(
            res.loc[valid, "P"],
            method="fdr_bh"
        )[1]

    res = res.sort_values(
        ["FDR", "P"],
        na_position="last"
    )

# ============================================================
# 7. SAVE
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

print("\nResults:")
print(res.to_string(index=False))

print("\nSaved:")
print(OUT)

print(QC_OUT)

print("\n" + "=" * 80)
print("FIRTH BURDEN-ADJUSTED ANALYSIS COMPLETED")
print("=" * 80)

