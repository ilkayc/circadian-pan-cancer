#!/usr/bin/env python3

import pandas as pd
import numpy as np

INPUT = "results/clock16/REAL_primary_CLOCK16_mutations.tsv.gz"
OUT = "results/clock16/CLOCK16_functional_annotation_QC.tsv"

print("=" * 80)
print("CLOCK16 FUNCTIONAL ANNOTATION QC")
print("=" * 80)

print("\nLoading:")
print(INPUT)

df = pd.read_csv(INPUT, sep="\t", compression="gzip", low_memory=False)

print("\nRows:", len(df))
print("Columns:")
print(df.columns.tolist())

# ------------------------------------------------------------
# REQUIRED COLUMNS
# ------------------------------------------------------------

required = [
    "patient_id",
    "gene",
    "effect",
    "effect_class",
    "SIFT",
    "PolyPhen"
]

missing = [x for x in required if x not in df.columns]

if missing:
    raise ValueError(f"Missing columns: {missing}")

# ------------------------------------------------------------
# EFFECT CLASS
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("EFFECT CLASS")
print("-" * 80)

print(df["effect_class"].value_counts(dropna=False).to_string())

# ------------------------------------------------------------
# RAW EFFECT
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("RAW EFFECT")
print("-" * 80)

print(
    df["effect"]
    .value_counts(dropna=False)
    .head(50)
    .to_string()
)

# ------------------------------------------------------------
# SIFT RAW VALUES
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("SIFT RAW VALUES")
print("-" * 80)

print("Missing:", df["SIFT"].isna().sum())
print("Non-missing:", df["SIFT"].notna().sum())

print(
    df["SIFT"]
    .astype(str)
    .value_counts(dropna=False)
    .head(50)
    .to_string()
)

# ------------------------------------------------------------
# POLYPHEN RAW VALUES
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("POLYPHEN RAW VALUES")
print("-" * 80)

print("Missing:", df["PolyPhen"].isna().sum())
print("Non-missing:", df["PolyPhen"].notna().sum())

print(
    df["PolyPhen"]
    .astype(str)
    .value_counts(dropna=False)
    .head(50)
    .to_string()
)

# ------------------------------------------------------------
# NUMERIC CONVERSION
# ------------------------------------------------------------

def numeric_annotation(x):

    if pd.isna(x):
        return np.nan

    s = str(x).strip()

    # Common forms:
    # "0.01"
    # "0.01(Damaging)"
    # "Damaging"
    # "deleterious(0.01)"
    # "benign(0.01)"

    import re

    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)

    if len(nums) == 0:
        return np.nan

    try:
        return float(nums[0])
    except:
        return np.nan


df["SIFT_numeric"] = df["SIFT"].apply(numeric_annotation)
df["PolyPhen_numeric"] = df["PolyPhen"].apply(numeric_annotation)

# ------------------------------------------------------------
# NONSYNONYMOUS ONLY
# ------------------------------------------------------------

nonsyn = df[df["effect_class"].astype(str).str.lower() == "nonsynonymous"].copy()

print("\n" + "-" * 80)
print("NONSYNONYMOUS MUTATIONS")
print("-" * 80)

print("Rows:", len(nonsyn))
print("Unique patients:", nonsyn["patient_id"].nunique())

# ------------------------------------------------------------
# SIFT / POLYPHEN NUMERIC SUMMARY
# ------------------------------------------------------------

for col in ["SIFT_numeric", "PolyPhen_numeric"]:

    print("\n", col)

    x = nonsyn[col].dropna()

    print("Non-missing:", len(x))

    if len(x) > 0:
        print(x.describe().to_string())

# ------------------------------------------------------------
# COMMON DAMAGE THRESHOLDS
# ------------------------------------------------------------

# SIFT:
# <= 0.05 conventionally considered deleterious
#
# PolyPhen:
# > 0.446 often considered possibly/probably damaging
# depending on annotation system.
#
# IMPORTANT:
# These are QC indicators here, NOT final classification.

nonsyn["SIFT_deleterious_QC"] = (
    nonsyn["SIFT_numeric"].notna()
    & (nonsyn["SIFT_numeric"] <= 0.05)
)

nonsyn["PolyPhen_damaging_QC"] = (
    nonsyn["PolyPhen_numeric"].notna()
    & (nonsyn["PolyPhen_numeric"] >= 0.446)
)

nonsyn["Either_damaging_QC"] = (
    nonsyn["SIFT_deleterious_QC"]
    | nonsyn["PolyPhen_damaging_QC"]
)

nonsyn["Both_damaging_QC"] = (
    nonsyn["SIFT_deleterious_QC"]
    & nonsyn["PolyPhen_damaging_QC"]
)

print("\n" + "-" * 80)
print("DAMAGE QC")
print("-" * 80)

print(
    "\nSIFT deleterious:",
    nonsyn["SIFT_deleterious_QC"].sum()
)

print(
    "PolyPhen damaging:",
    nonsyn["PolyPhen_damaging_QC"].sum()
)

print(
    "Either:",
    nonsyn["Either_damaging_QC"].sum()
)

print(
    "Both:",
    nonsyn["Both_damaging_QC"].sum()
)

# ------------------------------------------------------------
# UNIQUE PATIENT COUNTS
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("PATIENT-LEVEL DAMAGE COUNTS")
print("-" * 80)

for col in [
    "SIFT_deleterious_QC",
    "PolyPhen_damaging_QC",
    "Either_damaging_QC",
    "Both_damaging_QC"
]:

    patients = nonsyn.loc[nonsyn[col], "patient_id"].nunique()

    print(f"{col}: {patients} patients")

# ------------------------------------------------------------
# EFFECT CLASS × SIFT / POLYPHEN
# ------------------------------------------------------------

summary = (
    df.groupby("effect_class", dropna=False)
      .agg(
          Mutation_Rows=("patient_id", "size"),
          Unique_Patients=("patient_id", "nunique"),
          SIFT_NonMissing=("SIFT_numeric", lambda x: x.notna().sum()),
          PolyPhen_NonMissing=("PolyPhen_numeric", lambda x: x.notna().sum())
      )
      .reset_index()
)

print("\n" + "-" * 80)
print("ANNOTATION COMPLETENESS BY EFFECT CLASS")
print("-" * 80)

print(summary.to_string(index=False))

# ------------------------------------------------------------
# GENE LEVEL DAMAGE
# ------------------------------------------------------------

gene_summary = (
    nonsyn.groupby("gene")
    .agg(
        Nonsynonymous_Rows=("patient_id", "size"),
        Nonsynonymous_Patients=("patient_id", "nunique"),
        SIFT_Deleterious_Rows=("SIFT_deleterious_QC", "sum"),
        PolyPhen_Damaging_Rows=("PolyPhen_damaging_QC", "sum"),
        Either_Damaging_Rows=("Either_damaging_QC", "sum"),
        Both_Damaging_Rows=("Both_damaging_QC", "sum"),
        Either_Damaging_Patients=(
            "patient_id",
            lambda x: x[nonsyn.loc[x.index, "Either_damaging_QC"]].nunique()
        )
    )
    .reset_index()
    .sort_values("Either_Damaging_Rows", ascending=False)
)

print("\n" + "-" * 80)
print("GENE-LEVEL DAMAGE QC")
print("-" * 80)

print(gene_summary.to_string(index=False))

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

out_summary = summary.copy()

out_summary.to_csv(
    OUT,
    sep="\t",
    index=False
)

gene_summary.to_csv(
    "results/clock16/CLOCK16_functional_annotation_gene_QC.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUT)

print(
    "results/clock16/CLOCK16_functional_annotation_gene_QC.tsv"
)

print("\n" + "=" * 80)
print("FUNCTIONAL ANNOTATION QC COMPLETED")
print("=" * 80)
