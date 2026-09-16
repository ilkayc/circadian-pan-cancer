
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

# ============================================================
# CLOCK16 FUNCTIONAL MUTATION — POOLED PAN-CANCER FIRTH LR
#
# Endpoints:
#   1. Nonsynonymous
#   2. Damaging
#   3. Strong damaging
#
# Model:
#   Functional mutation status
#       ~ log1p(non-CLOCK16 mutation burden)
#       + cancer type
#
# Inference:
#   Firth logistic regression
#   Likelihood-ratio test
#   FDR correction
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
    "CLOCK16_functional_pooled_FIRTH_LR.tsv"
)

QC_OUT = (
    "results/clock16/"
    "CLOCK16_functional_pooled_FIRTH_LR_QC.tsv"
)


print("=" * 80)
print("CLOCK16 FUNCTIONAL MUTATION — POOLED PAN-CANCER FIRTH LR")
print("=" * 80)


# ============================================================
# Firth logistic regression
# ============================================================

def firth_logistic(X, y, max_iter=1000, tol=1e-9):

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)

    n, p = X.shape

    beta = np.zeros(p)

    def loglik(beta):

        eta = X @ beta

        # stable log(1 + exp(eta))
        log_term = np.logaddexp(0, eta)

        return np.sum(y * eta - log_term)

    converged = False

    for iteration in range(max_iter):

        eta = X @ beta

        # stable logistic
        mu = np.empty_like(eta)

        positive = eta >= 0
        mu[positive] = 1.0 / (1.0 + np.exp(-eta[positive]))

        exp_eta = np.exp(eta[~positive])
        mu[~positive] = exp_eta / (1.0 + exp_eta)

        W = mu * (1.0 - mu)

        # avoid zero weights
        W = np.maximum(W, 1e-12)

        XtWX = X.T @ (W[:, None] * X)

        try:
            XtWX_inv = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            XtWX_inv = np.linalg.pinv(XtWX)

        # Hat diagonal
        H = np.sum(
            (X @ XtWX_inv) * X,
            axis=1
        ) * W

        # Firth adjustment
        adjustment = H * (0.5 - mu)

        score = X.T @ (y - mu + adjustment)

        try:
            step = np.linalg.solve(XtWX, score)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(XtWX) @ score

        # line search
        current_ll = loglik(beta)

        step_factor = 1.0
        accepted = False

        while step_factor > 1e-8:

            beta_new = beta + step_factor * step

            new_ll = loglik(beta_new)

            if np.isfinite(new_ll) and new_ll >= current_ll:
                accepted = True
                break

            step_factor *= 0.5

        if not accepted:
            break

        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            converged = True
            break

        beta = beta_new

    # final information matrix
    eta = X @ beta

    mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -700, 700)))

    W = np.maximum(mu * (1.0 - mu), 1e-12)

    XtWX = X.T @ (W[:, None] * X)

    try:
        cov = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(XtWX)

    se = np.sqrt(np.maximum(np.diag(cov), 0))

    return {
        "beta": beta,
        "se": se,
        "loglik": loglik(beta),
        "converged": converged,
        "iterations": iteration + 1
    }


# ============================================================
# Load burden data
# ============================================================

print("\nLoading burden data:")

df = pd.read_csv(
    BURDEN_FILE,
    sep="\t"
)

print("Rows:", len(df))
print("Unique patients:", df["patient_id"].nunique())


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
        f"Missing burden columns: {missing}"
    )


# non-CLOCK16 burden
df["non_CLOCK16_mutations"] = (
    df["MC3_total_mutations"]
    - df["MC3_CLOCK16_mutations"]
)

if (df["non_CLOCK16_mutations"] < 0).any():
    raise ValueError(
        "Negative non-CLOCK16 mutation burden detected."
    )

df["log1p_non_CLOCK16"] = np.log1p(
    df["non_CLOCK16_mutations"]
)


# ============================================================
# Load CLOCK16 mutation data
# ============================================================

print("\nLoading CLOCK16 mutation data:")

mut = pd.read_csv(
    MUT_FILE,
    sep="\t"
)

print("Rows:", len(mut))
print("Unique patients:", mut["patient_id"].nunique())
print("Genes:", mut["gene"].nunique())


required_mut = [
    "patient_id",
    "gene",
    "effect_class"
]

missing_mut = [
    x for x in required_mut
    if x not in mut.columns
]

if missing_mut:
    raise ValueError(
        f"Missing mutation columns: {missing_mut}"
    )


# ============================================================
# Inspect effect classes
# ============================================================

print("\nEffect classes:")

print(
    mut["effect_class"]
    .value_counts(dropna=False)
    .to_string()
)


# ============================================================
# Define functional categories
# ============================================================

effect = (
    mut["effect_class"]
    .astype(str)
    .str.lower()
    .str.strip()
)

# Broad nonsynonymous definition
nonsynonymous_mask = effect.str.contains(
    "nonsyn",
    regex=False
)

# Damaging definition
damaging_mask = effect.str.contains(
    "damaging",
    regex=False
)

# Strong damaging definition
strong_damaging_mask = effect.str.contains(
    "strong",
    regex=False
) & effect.str.contains(
    "damaging",
    regex=False
)


# ============================================================
# Patient-level functional status
# ============================================================

functional = pd.DataFrame({
    "patient_id": df["patient_id"]
}).drop_duplicates()


def patient_status(mask):

    patients = set(
        mut.loc[mask, "patient_id"]
        .dropna()
        .astype(str)
    )

    return functional["patient_id"].astype(str).isin(
        patients
    ).astype(int)


functional["Nonsynonymous"] = patient_status(
    nonsynonymous_mask
)

functional["Damaging"] = patient_status(
    damaging_mask
)

functional["Strong_Damaging"] = patient_status(
    strong_damaging_mask
)


# merge
df["patient_id"] = df["patient_id"].astype(str)

functional["patient_id"] = (
    functional["patient_id"].astype(str)
)

df = df.merge(
    functional,
    on="patient_id",
    how="left"
)


# ============================================================
# Integrity checks
# ============================================================

print("\nFunctional patient counts:")

for endpoint in [
    "Nonsynonymous",
    "Damaging",
    "Strong_Damaging"
]:

    print(
        endpoint,
        "=",
        int(df[endpoint].sum())
    )


print("\nPatients outside denominator:")

mutation_patients = set(
    mut["patient_id"]
    .dropna()
    .astype(str)
)

denominator_patients = set(
    df["patient_id"].astype(str)
)

outside = mutation_patients - denominator_patients

print(len(outside))

if len(outside) > 0:
    raise ValueError(
        "Mutation patients outside denominator."
    )


# ============================================================
# Cancer design matrix
# ============================================================

cancer = pd.get_dummies(
    df["cancer_type"],
    drop_first=False,
    dtype=float
)

cancer_levels = sorted(
    cancer.columns.tolist()
)

reference = cancer_levels[0]

print("\nCancer reference:")
print(reference)

# Remove reference level
cancer = cancer.drop(
    columns=[reference]
)

X_cov = np.column_stack([
    np.ones(len(df)),
    df["log1p_non_CLOCK16"].values,
    cancer.values
])

print(
    "\nCovariate matrix:",
    X_cov.shape
)


# ============================================================
# Analysis
# ============================================================

endpoints = [
    "Nonsynonymous",
    "Damaging",
    "Strong_Damaging"
]

results = []
qc = []


for endpoint in endpoints:

    print("\n" + "-" * 80)
    print("ENDPOINT:", endpoint)
    print("-" * 80)

    y = df[endpoint].astype(float).values

    n_cases = int(y.sum())
    n_controls = int(len(y) - n_cases)

    print(
        "Cases:",
        n_cases,
        "Controls:",
        n_controls
    )

    # Reduced model:
    # functional status ~ cancer
    X_reduced = np.column_stack([
        np.ones(len(df)),
        cancer.values
    ])

    # Full model:
    # functional status ~ burden + cancer
    X_full = X_cov

    print(
        "Reduced parameters:",
        X_reduced.shape[1]
    )

    print(
        "Full parameters:",
        X_full.shape[1]
    )

    # fit reduced
    reduced = firth_logistic(
        X_reduced,
        y
    )

    # fit full
    full = firth_logistic(
        X_full,
        y
    )

    ll_reduced = reduced["loglik"]
    ll_full = full["loglik"]

    LR = 2.0 * (
        ll_full - ll_reduced
    )

    if LR < 0:
        LR = 0.0

    # one tested parameter: burden
    p = chi2.sf(
        LR,
        df=1
    )

    beta = full["beta"][1]
    se = full["se"][1]

    OR = np.exp(beta)

    CI_low = np.exp(
        beta - 1.96 * se
    )

    CI_high = np.exp(
        beta + 1.96 * se
    )

    results.append({
        "Endpoint": endpoint,
        "N": len(df),
        "Cases": n_cases,
        "Controls": n_controls,
        "OR_per_log1p_non_CLOCK16_mutation": OR,
        "CI95_Low": CI_low,
        "CI95_High": CI_high,
        "Firth_LR_Statistic": LR,
        "Firth_LR_P": p,
        "Direction": (
            "Positive"
            if beta > 0
            else "Negative"
        ),
        "Reduced_LogLik": ll_reduced,
        "Full_LogLik": ll_full,
        "Iterations_Reduced": reduced["iterations"],
        "Iterations_Full": full["iterations"]
    })

    qc.append({
        "Endpoint": endpoint,
        "N": len(df),
        "Cases": n_cases,
        "Controls": n_controls,
        "Reduced_Converged": reduced["converged"],
        "Full_Converged": full["converged"],
        "LR_Valid": bool(
            np.isfinite(LR)
            and np.isfinite(p)
        )
    })


# ============================================================
# FDR
# ============================================================

res = pd.DataFrame(results)

res["FDR"] = multipletests(
    res["Firth_LR_P"],
    method="fdr_bh"
)[1]

res = res.sort_values(
    "Firth_LR_P"
).reset_index(drop=True)


# ============================================================
# Output
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
            "OR_per_log1p_non_CLOCK16_mutation",
            "CI95_Low",
            "CI95_High",
            "Firth_LR_P",
            "FDR",
            "Direction"
        ]
    ].to_string(index=False)
)


print("\nFDR < 0.05:")

sig = res[
    res["FDR"] < 0.05
]

print(
    sig[
        [
            "Endpoint",
            "Cases",
            "OR_per_log1p_non_CLOCK16_mutation",
            "CI95_Low",
            "CI95_High",
            "Firth_LR_P",
            "FDR"
        ]
    ].to_string(index=False)
)


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


print("\nSaved:")
print(OUT)
print(QC_OUT)

print("\n" + "=" * 80)
print("FUNCTIONAL POOLED FIRTH LR COMPLETED")
print("=" * 80)

