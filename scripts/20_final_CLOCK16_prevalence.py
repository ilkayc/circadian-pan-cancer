import pandas as pd
import numpy as np
import os

INPUT = "results/clock16/PAN_CANCER_patient_mutational_burden.tsv"
OUT = "results/clock16/CLOCK16_final_cancer_prevalence.tsv"
QC_OUT = "results/clock16/CLOCK16_final_prevalence_QC.tsv"

print("=" * 80)
print("FINAL CLOCK16 PATIENT-LEVEL PREVALENCE")
print("=" * 80)

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------

df = pd.read_csv(INPUT, sep="\t")

print("\nInput:", INPUT)
print("Rows:", len(df))
print("Unique patients:", df["patient_id"].nunique())
print("Cancer types:", df["cancer_type"].nunique())

# ------------------------------------------------------------
# Basic integrity
# ------------------------------------------------------------

assert df["patient_id"].nunique() == len(df), \
    "ERROR: duplicate patient IDs detected"

assert df["CLOCK16_mutated"].isin([True, False]).all(), \
    "ERROR: unexpected CLOCK16 status"

# ------------------------------------------------------------
# Cancer-specific prevalence
# ------------------------------------------------------------

rows = []

for cancer, sub in df.groupby("cancer_type", sort=False):

    n = len(sub)

    mutated = int(sub["CLOCK16_mutated"].sum())
    wt = n - mutated

    prevalence = mutated / n
    prevalence_pct = prevalence * 100

    # Wilson 95% CI
    z = 1.96
    phat = prevalence

    denom = 1 + z**2 / n

    center = (
        phat + z**2 / (2 * n)
    ) / denom

    half = (
        z * np.sqrt(
            (phat * (1 - phat) / n)
            + z**2 / (4 * n**2)
        )
    ) / denom

    ci_low = max(0, center - half)
    ci_high = min(1, center + half)

    rows.append({
        "Cancer_Type": cancer,
        "N": n,
        "CLOCK16_Mutated": mutated,
        "CLOCK16_Wild_Type": wt,
        "Prevalence": prevalence,
        "Prevalence_Percent": prevalence_pct,
        "CI95_Low": ci_low,
        "CI95_High": ci_high,
        "CI95_Low_Percent": ci_low * 100,
        "CI95_High_Percent": ci_high * 100
    })

prev = pd.DataFrame(rows)

# ------------------------------------------------------------
# Sort by prevalence
# ------------------------------------------------------------

prev = prev.sort_values(
    ["Prevalence", "Cancer_Type"],
    ascending=[False, True]
).reset_index(drop=True)

# ------------------------------------------------------------
# Global totals
# ------------------------------------------------------------

total_n = len(df)
total_mutated = int(df["CLOCK16_mutated"].sum())
total_wt = total_n - total_mutated

total_prevalence = total_mutated / total_n

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

qc = pd.DataFrame([
    ["Total_Patients", total_n],
    ["Unique_Patients", df["patient_id"].nunique()],
    ["Cancer_Types", df["cancer_type"].nunique()],
    ["Total_CLOCK16_Mutated", total_mutated],
    ["Total_CLOCK16_Wild_Type", total_wt],
    ["Sum_Cancer_Patients", int(prev["N"].sum())],
    ["Sum_Cancer_Mutated", int(prev["CLOCK16_Mutated"].sum())],
    ["Sum_Cancer_Wild_Type", int(prev["CLOCK16_Wild_Type"].sum())],
    ["Overall_Prevalence", total_prevalence],
    ["Overall_Prevalence_Percent", total_prevalence * 100],
    ["Invalid_Prevalence", int(
        ((prev["Prevalence"] < 0) | (prev["Prevalence"] > 1)).sum()
    )],
    ["Invalid_CI", int(
        (
            (prev["CI95_Low"] < 0)
            | (prev["CI95_High"] > 1)
            | (prev["CI95_Low"] > prev["CI95_High"])
        ).sum()
    )]
], columns=["QC_Metric", "Value"])

# ------------------------------------------------------------
# Final assertions
# ------------------------------------------------------------

assert total_n == 10593, \
    f"ERROR: expected 10593 patients, found {total_n}"

assert total_mutated == 1165, \
    f"ERROR: expected 1165 CLOCK16-mutated patients, found {total_mutated}"

assert total_wt == 9428, \
    f"ERROR: expected 9428 CLOCK16-wild-type patients, found {total_wt}"

assert prev["N"].sum() == total_n
assert prev["CLOCK16_Mutated"].sum() == total_mutated
assert prev["CLOCK16_Wild_Type"].sum() == total_wt

assert qc.loc[
    qc["QC_Metric"] == "Invalid_Prevalence", "Value"
].iloc[0] == 0

assert qc.loc[
    qc["QC_Metric"] == "Invalid_CI", "Value"
].iloc[0] == 0

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

prev.to_csv(OUT, sep="\t", index=False)
qc.to_csv(QC_OUT, sep="\t", index=False)

# ------------------------------------------------------------
# Print final QC
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("FINAL QC")
print("=" * 80)

for _, row in qc.iterrows():
    print(f"{row['QC_Metric']}: {row['Value']}")

print("\n" + "=" * 80)
print("CANCER-SPECIFIC PREVALENCE")
print("=" * 80)

print(
    prev[
        [
            "Cancer_Type",
            "N",
            "CLOCK16_Mutated",
            "CLOCK16_Wild_Type",
            "Prevalence_Percent",
            "CI95_Low_Percent",
            "CI95_High_Percent"
        ]
    ].to_string(index=False)
)

print("\nSaved:")
print(OUT)
print(QC_OUT)

print("\nSTATUS: FINAL PREVALENCE QC PASSED")
