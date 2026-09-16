import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

INPUT = "results/clock16/PAN_CANCER_patient_mutational_burden.tsv"

OUT_ASSOC = "results/clock16/PAN_CANCER_burden_adjusted_association.tsv"
OUT_SUMMARY = "results/clock16/PAN_CANCER_burden_adjusted_summary.tsv"
OUT_QC = "results/clock16/PAN_CANCER_burden_adjusted_QC.tsv"

print("="*80)
print("CLOCK16 MUTATION ASSOCIATION — MUTATIONAL BURDEN ADJUSTMENT")
print("="*80)

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

# ------------------------------------------------------------
# 1. CREATE NON-CLOCK16 MUTATION BURDEN
# ------------------------------------------------------------

df["MC3_non_CLOCK16_mutations"] = (
    df["MC3_total_mutations"] -
    df["MC3_CLOCK16_mutations"]
)

if (df["MC3_non_CLOCK16_mutations"] < 0).any():
    raise ValueError(
        "ERROR: negative non-CLOCK16 mutation burden detected."
    )

# log1p transformation because burden is highly right-skewed
df["log_MC3_total"] = np.log1p(df["MC3_total_mutations"])
df["log_MC3_non_CLOCK16"] = np.log1p(
    df["MC3_non_CLOCK16_mutations"]
)

# ------------------------------------------------------------
# 2. BASIC GLOBAL COMPARISON
# ------------------------------------------------------------

mut = df[df["CLOCK16_mutated"] == True]
wt = df[df["CLOCK16_mutated"] == False]

print("\nCLOCK16-mutated:", len(mut))
print("CLOCK16-wild-type:", len(wt))

print("\nGlobal mutation burden:")

print(
    "Mutated median:",
    mut["MC3_total_mutations"].median()
)

print(
    "Wild-type median:",
    wt["MC3_total_mutations"].median()
)

u, p = mannwhitneyu(
    mut["MC3_total_mutations"],
    wt["MC3_total_mutations"],
    alternative="two-sided"
)

print("Mann-Whitney P:", p)

# ------------------------------------------------------------
# 3. CANCER-SPECIFIC LOGISTIC REGRESSION
#
# CLOCK16 status ~ log(non-CLOCK16 burden)
#
# This asks whether CLOCK16 mutation status is associated
# with mutation burden within each cancer type.
# ------------------------------------------------------------

results = []

for cancer, sub in df.groupby("cancer_type"):

    n = len(sub)
    n_mut = int(sub["CLOCK16_mutated"].sum())
    n_wt = n - n_mut

    if n_mut == 0 or n_wt == 0:
        continue

    x = sub["log_MC3_non_CLOCK16"].astype(float)
    y = sub["CLOCK16_mutated"].astype(int)

    X = pd.DataFrame({
        "log_non_CLOCK16_burden": x
    })

    X = sm.add_constant(X)

    try:

        model = sm.Logit(y, X).fit(
            disp=False,
            maxiter=200
        )

        beta = model.params["log_non_CLOCK16_burden"]
        se = model.bse["log_non_CLOCK16_burden"]

        OR = np.exp(beta)

        CI_low = np.exp(beta - 1.96 * se)
        CI_high = np.exp(beta + 1.96 * se)

        pval = model.pvalues["log_non_CLOCK16_burden"]

        results.append({
            "Cancer_Type": cancer,
            "N": n,
            "CLOCK16_Mutated": n_mut,
            "CLOCK16_WildType": n_wt,
            "OR_per_log1p_non_CLOCK16_mutation": OR,
            "CI95_Low": CI_low,
            "CI95_High": CI_high,
            "P": pval
        })

    except Exception as e:

        print(
            f"WARNING: model failed for {cancer}: {e}"
        )

# ------------------------------------------------------------
# 4. FDR
# ------------------------------------------------------------

res = pd.DataFrame(results)

if len(res) > 0:

    res["FDR"] = multipletests(
        res["P"],
        method="fdr_bh"
    )[1]

    res["Direction"] = np.where(
        res["OR_per_log1p_non_CLOCK16_mutation"] > 1,
        "Positive",
        "Negative"
    )

    res = res.sort_values(
        ["FDR", "P"]
    )

# ------------------------------------------------------------
# 5. SUMMARY
# ------------------------------------------------------------

summary = pd.DataFrame([
    ["Patients", len(df)],
    ["Cancer_types", df["cancer_type"].nunique()],
    ["CLOCK16_mutated_patients", len(mut)],
    ["CLOCK16_wildtype_patients", len(wt)],
    ["Global_MannWhitney_P", p],
    ["Cancer_types_tested", len(res)],
    [
        "FDR<0.05",
        int((res["FDR"] < 0.05).sum())
        if len(res) else 0
    ],
    [
        "FDR<0.10",
        int((res["FDR"] < 0.10).sum())
        if len(res) else 0
    ]
], columns=["Metric", "Value"])

# ------------------------------------------------------------
# 6. QC
# ------------------------------------------------------------

qc = pd.DataFrame([
    ["Input_rows", len(df)],
    ["Unique_patients", df["patient_id"].nunique()],
    ["Cancer_types", df["cancer_type"].nunique()],
    ["Missing_values_total", int(df.isna().sum().sum())],
    [
        "Negative_non_CLOCK16_burden",
        int((df["MC3_non_CLOCK16_mutations"] < 0).sum())
    ],
    [
        "CLOCK16_status_true",
        int(df["CLOCK16_mutated"].sum())
    ],
    [
        "CLOCK16_status_false",
        int((~df["CLOCK16_mutated"]).sum())
    ]
], columns=["Metric", "Value"])

# ------------------------------------------------------------
# 7. SAVE
# ------------------------------------------------------------

res.to_csv(
    OUT_ASSOC,
    sep="\t",
    index=False
)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False
)

qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False
)

print("\n" + "="*80)
print("RESULTS")
print("="*80)

if len(res):

    print(
        res[
            [
                "Cancer_Type",
                "N",
                "CLOCK16_Mutated",
                "OR_per_log1p_non_CLOCK16_mutation",
                "CI95_Low",
                "CI95_High",
                "P",
                "FDR",
                "Direction"
            ]
        ].to_string(index=False)
    )

print("\nSaved:")
print(OUT_ASSOC)
print(OUT_SUMMARY)
print(OUT_QC)

print("\n" + "="*80)
print("BURDEN-ADJUSTED ANALYSIS COMPLETED")
print("="*80)
