#!/usr/bin/env python3

import pandas as pd
import numpy as np

INPUT = "results/clock16/PAN_CANCER_patient_mutational_burden.tsv"
OUT = "results/clock16/PAN_CANCER_burden_integrity_check.tsv"

print("=" * 80)
print("CLOCK16 NON-CLOCK16 MUTATION BURDEN INTEGRITY CHECK")
print("=" * 80)

df = pd.read_csv(INPUT, sep="\t")

print("\nColumns:")
print(df.columns.tolist())

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

# ------------------------------------------------------------
# 1. NON-CLOCK16 BURDEN
# ------------------------------------------------------------

df["non_CLOCK16_mutations"] = (
    df["MC3_total_mutations"]
    - df["MC3_CLOCK16_mutations"]
)

# ------------------------------------------------------------
# 2. INTEGRITY
# ------------------------------------------------------------

print("\nNegative non-CLOCK16 burden:")
print((df["non_CLOCK16_mutations"] < 0).sum())

if (df["non_CLOCK16_mutations"] < 0).any():
    raise ValueError("Negative mutation burden detected.")

print("\nCLOCK16 mutation count summary:")
print(df["MC3_CLOCK16_mutations"].describe())

print("\nNon-CLOCK16 mutation burden summary:")
print(df["non_CLOCK16_mutations"].describe())

# ------------------------------------------------------------
# 3. MUTATED VS WT
# ------------------------------------------------------------

print("\nMutation burden by CLOCK16 status:")

summary = (
    df
    .groupby("CLOCK16_mutated")["non_CLOCK16_mutations"]
    .agg(
        N="size",
        Median="median",
        Mean="mean",
        Q1=lambda x: x.quantile(0.25),
        Q3=lambda x: x.quantile(0.75),
        Max="max"
    )
)

print(summary.to_string())

# ------------------------------------------------------------
# 4. LOG1P BURDEN
# ------------------------------------------------------------

df["log1p_non_CLOCK16_burden"] = np.log1p(
    df["non_CLOCK16_mutations"]
)

print("\nlog1p(non-CLOCK16 burden):")
print(
    df.groupby("CLOCK16_mutated")["log1p_non_CLOCK16_burden"]
      .describe()
      .to_string()
)

# ------------------------------------------------------------
# 5. CANCER-SPECIFIC COUNTS
# ------------------------------------------------------------

cancer_summary = (
    df.groupby("cancer_type")
      .agg(
          N=("patient_id", "size"),
          CLOCK16_Mutated=("CLOCK16_mutated", "sum"),
          Median_NonCLOCK16_Burden=(
              "non_CLOCK16_mutations",
              "median"
          ),
          Mean_NonCLOCK16_Burden=(
              "non_CLOCK16_mutations",
              "mean"
          )
      )
      .reset_index()
)

cancer_summary["CLOCK16_Mutation_Rate"] = (
    cancer_summary["CLOCK16_Mutated"]
    / cancer_summary["N"]
)

print("\nCancer-specific summary:")
print(cancer_summary.to_string(index=False))

# ------------------------------------------------------------
# 6. SAVE
# ------------------------------------------------------------

df.to_csv(
    OUT,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUT)

print("\n" + "=" * 80)
print("INTEGRITY CHECK COMPLETED")
print("=" * 80)
