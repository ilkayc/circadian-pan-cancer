#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
from scipy.optimize import minimize
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

# ============================================================
# INPUT / OUTPUT
# ============================================================

BURDEN = "results/clock16/PAN_CANCER_patient_mutational_burden.tsv"
MUT = "results/clock16/REAL_primary_CLOCK16_mutations.tsv.gz"

OUT = "results/clock16/CLOCK16_functional_CORRECTED_pooled_FIRTH_LR.tsv"
QC_OUT = "results/clock16/CLOCK16_functional_CORRECTED_pooled_FIRTH_LR_QC.tsv"

# ============================================================
# FIRTH LOGISTIC REGRESSION
# ============================================================

def firth_logistic(X, y, max_iter=1000):

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)

    p = X.shape[1]

    beta0 = np.zeros(p)

    def objective(beta):

        eta = X @ beta

        # stable log-likelihood
        loglik = np.sum(
            y * eta - np.logaddexp(0, eta)
        )

        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -700, 700)))

        W = mu * (1.0 - mu)

        XWX = X.T @ (W[:, None] * X)

        sign, logdet = np.linalg.slogdet(
            XWX + np.eye(p) * 1e-10
        )

        if sign <= 0:
            return 1e100

        penalty = 0.5 * logdet

        return -(loglik + penalty)

    result = minimize(
        objective,
        beta0,
        method="BFGS",
        options={
            "maxiter": max_iter,
            "gtol": 1e-8
        }
    )

    beta = result.x

    # numerical Hessian approximation through weighted information
    eta = X @ beta

    mu = 1.0 / (
        1.0 + np.exp(-np.clip(eta, -700, 700))
    )

    W = mu * (1.0 - mu)

    XWX = X.T @ (W[:, None] * X)

    cov = np.linalg.pinv(
        XWX + np.eye(p) * 1e-10
    )

    se = np.sqrt(
        np.maximum(np.diag(cov), 0)
    )

    return result, beta, se


# ============================================================
# PARSE SIFT / POLYPHEN
# ============================================================

def parse_numeric(x):

    if pd.isna(x):
        return np.nan

    s = str(x).strip()

    nums = re.findall(
        r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?",
        s
    )

    if not nums:
        return np.nan

    try:
        return float(nums[0])
    except:
        return np.nan


# ============================================================
# LOAD BURDEN
# ============================================================

print("=" * 80)
print("CLOCK16 CORRECTED FUNCTIONAL MUTATION")
print("POOLED PAN-CANCER FIRTH LR")
print("=" * 80)

burden = pd.read_csv(
    BURDEN,
    sep="\t"
)

print("\nLoading burden:")
print("Rows:", len(burden))
print("Patients:", burden["patient_id"].nunique())

# ------------------------------------------------------------
# CLOCK16 status
# ------------------------------------------------------------

if burden["CLOCK16_mutated"].dtype != bool:

    burden["CLOCK16_mutated"] = (
        burden["CLOCK16_mutated"]
        .astype(str)
        .str.lower()
        .isin(["true", "1", "yes"])
    )

# non-CLOCK16 burden
burden["non_CLOCK16_mutations"] = (
    burden["MC3_total_mutations"]
    - burden["MC3_CLOCK16_mutations"]
)

if (burden["non_CLOCK16_mutations"] < 0).any():
    raise ValueError(
        "Negative non-CLOCK16 mutation burden detected."
    )

burden["log1p_non_CLOCK16"] = np.log1p(
    burden["non_CLOCK16_mutations"]
)

# ============================================================
# LOAD MUTATIONS
# ============================================================

print("\nLoading CLOCK16 mutations:")

mut = pd.read_csv(
    MUT,
    sep="\t",
    compression="gzip",
    low_memory=False
)

print("Rows:", len(mut))
print("Genes:", mut["gene"].nunique())
print("Patients:", mut["patient_id"].nunique())

# ============================================================
# PARSE ANNOTATIONS
# ============================================================

mut["SIFT_numeric"] = mut["SIFT"].apply(
    parse_numeric
)

mut["PolyPhen_numeric"] = mut["PolyPhen"].apply(
    parse_numeric
)

# ------------------------------------------------------------
# Nonsynonymous
# ------------------------------------------------------------

nonsyn = mut[
    mut["effect_class"]
    .astype(str)
    .str.lower()
    .eq("nonsynonymous")
].copy()

# ============================================================
# FUNCTIONAL DEFINITIONS
# ============================================================

# SIFT:
# <= 0.05 = deleterious
#
# PolyPhen:
# >= 0.446 = damaging
#
# Damaging:
# either predictor positive
#
# Strong_Damaging:
# both predictors positive

nonsyn["SIFT_deleterious"] = (
    nonsyn["SIFT_numeric"].notna()
    & (nonsyn["SIFT_numeric"] <= 0.05)
)

nonsyn["PolyPhen_damaging"] = (
    nonsyn["PolyPhen_numeric"].notna()
    & (nonsyn["PolyPhen_numeric"] >= 0.446)
)

nonsyn["Damaging"] = (
    nonsyn["SIFT_deleterious"]
    | nonsyn["PolyPhen_damaging"]
)

nonsyn["Strong_Damaging"] = (
    nonsyn["SIFT_deleterious"]
    & nonsyn["PolyPhen_damaging"]
)

print("\nFunctional mutation rows:")
print(
    "Nonsynonymous:",
    len(nonsyn)
)

print(
    "Damaging:",
    nonsyn["Damaging"].sum()
)

print(
    "Strong_Damaging:",
    nonsyn["Strong_Damaging"].sum()
)

print("\nFunctional patient counts:")

for endpoint in [
    "Damaging",
    "Strong_Damaging"
]:

    print(
        endpoint,
        "=",
        nonsyn.loc[
            nonsyn[endpoint],
            "patient_id"
        ].nunique()
    )

# ============================================================
# PATIENT LEVEL ENDPOINTS
# ============================================================

patients = burden[
    [
        "patient_id",
        "cancer_type",
        "log1p_non_CLOCK16"
    ]
].copy()

# ------------------------------------------------------------
# Nonsynonymous
# ------------------------------------------------------------

nonsyn_patients = set(
    nonsyn["patient_id"].unique()
)

patients["Nonsynonymous"] = (
    patients["patient_id"]
    .isin(nonsyn_patients)
)

# ------------------------------------------------------------
# Damaging
# ------------------------------------------------------------

damaging_patients = set(
    nonsyn.loc[
        nonsyn["Damaging"],
        "patient_id"
    ].unique()
)

patients["Damaging"] = (
    patients["patient_id"]
    .isin(damaging_patients)
)

# ------------------------------------------------------------
# Strong damaging
# ------------------------------------------------------------

strong_patients = set(
    nonsyn.loc[
        nonsyn["Strong_Damaging"],
        "patient_id"
    ].unique()
)

patients["Strong_Damaging"] = (
    patients["patient_id"]
    .isin(strong_patients)
)

# ============================================================
# CANCER DESIGN MATRIX
# ============================================================

cancer = pd.Categorical(
    patients["cancer_type"]
)

reference = sorted(
    patients["cancer_type"].unique()
)[0]

print("\nCancer reference:")
print(reference)

cancer_dummies = pd.get_dummies(
    patients["cancer_type"],
    drop_first=True,
    dtype=float
)

X_base = pd.DataFrame(
    {
        "Intercept": 1.0,
        "log1p_non_CLOCK16": patients[
            "log1p_non_CLOCK16"
        ].values
    }
)

X_base = pd.concat(
    [
        X_base,
        cancer_dummies.reset_index(drop=True)
    ],
    axis=1
)

X_base = X_base.astype(float)

X_reduced = X_base.drop(
    columns=["log1p_non_CLOCK16"]
)

# ============================================================
# FIRTH LR TEST
# ============================================================

def fit_endpoint(y):

    y = np.asarray(y, dtype=float)

    # Full model
    full_res, beta, se = firth_logistic(
        X_base.values,
        y
    )

    # Reduced model
    red_res, beta_red, se_red = firth_logistic(
        X_reduced.values,
        y
    )

    if not full_res.success:
        print(
            "WARNING: full model:",
            full_res.message
        )

    if not red_res.success:
        print(
            "WARNING: reduced model:",
            red_res.message
        )

    # Penalized objective returns negative penalized log-likelihood
    ll_full = -full_res.fun
    ll_red = -red_res.fun

    LR = max(
        0.0,
        2.0 * (ll_full - ll_red)
    )

    P = chi2.sf(
        LR,
        df=1
    )

    # log1p burden coefficient is column 1
    beta_burden = beta[1]
    se_burden = se[1]

    OR = np.exp(
        np.clip(beta_burden, -700, 700)
    )

    CI_low = np.exp(
        np.clip(
            beta_burden - 1.96 * se_burden,
            -700,
            700
        )
    )

    CI_high = np.exp(
        np.clip(
            beta_burden + 1.96 * se_burden,
            -700,
            700
        )
    )

    return {
        "OR": OR,
        "CI_low": CI_low,
        "CI_high": CI_high,
        "LR": LR,
        "P": P,
        "full_success": full_res.success,
        "reduced_success": red_res.success
    }


# ============================================================
# RUN ENDPOINTS
# ============================================================

endpoints = [
    "Nonsynonymous",
    "Damaging",
    "Strong_Damaging"
]

results = []

for endpoint in endpoints:

    y = patients[endpoint].astype(int).values

    cases = int(y.sum())
    controls = int(len(y) - cases)

    print("\n" + "-" * 80)
    print("ENDPOINT:", endpoint)
    print("Cases:", cases)
    print("Controls:", controls)
    print("-" * 80)

    # Need both outcome classes
    if cases == 0 or controls == 0:

        print(
            "SKIPPED: endpoint has only one outcome class."
        )

        results.append(
            {
                "Endpoint": endpoint,
                "N": len(y),
                "Cases": cases,
                "Controls": controls,
                "OR_per_log1p_non_CLOCK16_mutation": np.nan,
                "CI95_Low": np.nan,
                "CI95_High": np.nan,
                "Firth_LR_Statistic": np.nan,
                "Firth_LR_P": np.nan,
                "Status": "No outcome variation"
            }
        )

        continue

    fit = fit_endpoint(y)

    direction = (
        "Positive"
        if fit["OR"] > 1
        else "Negative"
    )

    results.append(
        {
            "Endpoint": endpoint,
            "N": len(y),
            "Cases": cases,
            "Controls": controls,
            "OR_per_log1p_non_CLOCK16_mutation": fit["OR"],
            "CI95_Low": fit["CI_low"],
            "CI95_High": fit["CI_high"],
            "Firth_LR_Statistic": fit["LR"],
            "Firth_LR_P": fit["P"],
            "Direction": direction,
            "Full_Model_Success": fit["full_success"],
            "Reduced_Model_Success": fit["reduced_success"]
        }
    )

# ============================================================
# FDR
# ============================================================

res = pd.DataFrame(results)

valid = res["Firth_LR_P"].notna()

res["FDR"] = np.nan

if valid.sum() > 0:

    res.loc[
        valid,
        "FDR"
    ] = multipletests(
        res.loc[
            valid,
            "Firth_LR_P"
        ],
        method="fdr_bh"
    )[1]

# ============================================================
# OUTPUT
# ============================================================

print("\n" + "=" * 80)
print("RESULTS")
print("=" * 80)

print(
    res[
        [
            "Endpoint",
            "N",
            "Cases",
            "Controls",
            "OR_per_log1p_non_CLOCK16_mutation",
            "CI95_Low",
            "CI95_High",
            "Firth_LR_P",
            "FDR",
            "Direction"
        ]
    ].to_string(index=False)
)

res.to_csv(
    OUT,
    sep="\t",
    index=False
)

qc = pd.DataFrame(
    {
        "Metric": [
            "Total patients",
            "Total mutation rows",
            "Nonsynonymous rows",
            "Nonsynonymous patients",
            "Damaging mutation rows",
            "Damaging patients",
            "Strong damaging mutation rows",
            "Strong damaging patients",
            "Models attempted",
            "Models completed",
            "FDR < 0.05"
        ],
        "Value": [
            len(patients),
            len(mut),
            len(nonsyn),
            nonsyn_patients.__len__(),
            int(nonsyn["Damaging"].sum()),
            len(damaging_patients),
            int(nonsyn["Strong_Damaging"].sum()),
            len(strong_patients),
            len(endpoints),
            int(
                res["Firth_LR_P"].notna().sum()
            ),
            int(
                (
                    res["FDR"] < 0.05
                ).sum()
            )
        ]
    }
)

qc.to_csv(
    QC_OUT,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUT)
print(QC_OUT)

print("\n" + "=" * 80)
print("CORRECTED FUNCTIONAL FIRTH LR COMPLETED")
print("=" * 80)

