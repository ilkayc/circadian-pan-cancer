#!/usr/bin/env python3

import os
import itertools
import numpy as np
import pandas as pd

# ============================================================
# CLOCK16 PAN-CANCER PAIRWISE CO-MUTATION
# PATIENT-LEVEL + CANCER MAPPING + FIXED-MARGIN PERMUTATION
# ============================================================

PROJECT = os.getcwd()
RESULTS = os.path.join(PROJECT, "results", "clock16")
DOWNLOADS = os.path.expanduser("~/Downloads/CLOCK16_FINAL")
OUTDIR = os.path.join(DOWNLOADS, "Pairwise")

os.makedirs(RESULTS, exist_ok=True)
os.makedirs(OUTDIR, exist_ok=True)

CLOCK16 = [
    "TIMELESS", "NPAS2", "CRY1", "CRY2",
    "ARNTL", "CSNK1D", "CSNK1E",
    "PER1", "PER2", "PER3",
    "RORA", "RORB", "RORC",
    "CLOCK", "NR1D1", "NR1D2"
]

EVENT_FILE = os.path.join(
    RESULTS,
    "REAL_CLOCK16_multiple_same_gene_events.tsv"
)

print("=" * 90)
print("CLOCK16 PAN-CANCER PAIRWISE CO-MUTATION")
print("=" * 90)

# ============================================================
# 1. LOAD CLOCK16 PATIENT-GENE EVENTS
# ============================================================

if not os.path.exists(EVENT_FILE):
    raise FileNotFoundError(EVENT_FILE)

events = pd.read_csv(
    EVENT_FILE,
    sep="\t",
    low_memory=False
)

events.columns = [
    str(c).strip().lower().replace(" ", "_")
    for c in events.columns
]

required = {"patient_id", "gene"}

if not required.issubset(events.columns):
    raise ValueError(
        f"Required columns missing. Available: {events.columns.tolist()}"
    )

events["patient_id"] = events["patient_id"].astype(str).str.strip()
events["gene"] = (
    events["gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

events = events[
    events["gene"].isin(CLOCK16)
].copy()

events = events[
    ["patient_id", "gene"]
].drop_duplicates()

print("\nCLOCK16 event input:")
print(EVENT_FILE)
print("Shape:", events.shape)
print("Patients:", events["patient_id"].nunique())
print("CLOCK16 genes:", events["gene"].nunique())

# ============================================================
# 2. FIND PATIENT -> CANCER TYPE METADATA
# ============================================================

print("\nSearching for patient-level cancer metadata...")

metadata_candidates = [
    os.path.join(
        RESULTS,
        "PAN_CANCER_patient_mutational_burden.tsv"
    ),
    os.path.join(
        RESULTS,
        "PAN_CANCER_burden_integrity_check.tsv"
    ),
]

# Search additional files if necessary
for fn in sorted(os.listdir(RESULTS)):
    if fn.endswith(".tsv"):
        path = os.path.join(RESULTS, fn)
        if path not in metadata_candidates:
            metadata_candidates.append(path)

metadata = None
metadata_path = None

for path in metadata_candidates:

    if not os.path.exists(path):
        continue

    try:
        tmp = pd.read_csv(
            path,
            sep="\t",
            low_memory=False
        )
    except Exception:
        continue

    tmp.columns = [
        str(c).strip().lower().replace(" ", "_")
        for c in tmp.columns
    ]

    if "patient_id" not in tmp.columns:
        continue

    cancer_col = None

    for c in [
        "cancer_type",
        "cancer",
        "project",
        "tumor_type"
    ]:
        if c in tmp.columns:
            cancer_col = c
            break

    if cancer_col is None:
        continue

    tmp["patient_id"] = (
        tmp["patient_id"]
        .astype(str)
        .str.strip()
    )

    tmp[cancer_col] = (
        tmp[cancer_col]
        .astype(str)
        .str.strip()
    )

    overlap = (
        set(events["patient_id"])
        & set(tmp["patient_id"])
    )

    if len(overlap) == 0:
        continue

    metadata = tmp[
        ["patient_id", cancer_col]
    ].drop_duplicates()

    metadata = metadata.rename(
        columns={cancer_col: "cancer_type"}
    )

    metadata_path = path

    print("Selected metadata:")
    print(metadata_path)
    print("Patient overlap:", len(overlap))
    print("Cancer types:", metadata["cancer_type"].nunique())

    break

if metadata is None:
    raise ValueError(
        "Could not find patient_id -> cancer_type metadata."
    )

# ============================================================
# 3. MERGE PATIENT-GENE EVENTS WITH CANCER TYPE
# ============================================================

events = events.merge(
    metadata,
    on="patient_id",
    how="left"
)

missing_cancer = events["cancer_type"].isna().sum()

print("\nPatient-to-cancer mapping:")
print("Events without cancer type:", missing_cancer)

if missing_cancer > 0:
    print(
        events[
            events["cancer_type"].isna()
        ]["patient_id"]
        .drop_duplicates()
        .head(20)
        .to_string(index=False)
    )

events = events.dropna(
    subset=["cancer_type"]
).copy()

print(
    "Mapped CLOCK16 events:",
    len(events)
)

print(
    "Mapped patients:",
    events["patient_id"].nunique()
)

print(
    "Mapped cancer types:",
    events["cancer_type"].nunique()
)

# ============================================================
# 4. PATIENT × CANCER × GENE BINARY MATRIX
# ============================================================

patient_cancer = (
    events[
        ["patient_id", "cancer_type"]
    ]
    .drop_duplicates()
)

# Check patient -> cancer uniqueness
patient_n_cancer = (
    patient_cancer
    .groupby("patient_id")["cancer_type"]
    .nunique()
)

if (patient_n_cancer > 1).any():

    bad = int(
        (patient_n_cancer > 1).sum()
    )

    raise ValueError(
        f"{bad} patients map to multiple cancer types."
    )

patient_cancer = (
    patient_cancer
    .set_index("patient_id")
)

binary = (
    events
    .assign(value=1)
    .pivot_table(
        index="patient_id",
        columns="gene",
        values="value",
        aggfunc="max",
        fill_value=0
    )
)

for gene in CLOCK16:
    if gene not in binary.columns:
        binary[gene] = 0

binary = binary[CLOCK16].astype(int)

binary = binary.join(
    patient_cancer,
    how="inner"
)

print("\nFinal analysis matrix:")
print("Patients:", len(binary))
print(
    "Cancer types:",
    binary["cancer_type"].nunique()
)

print(
    "CLOCK16-mutated patients:",
    int(
        (
            binary[CLOCK16].sum(axis=1) > 0
        ).sum()
    )
)

# ============================================================
# 5. PAN-CANCER FIXED-MARGIN PERMUTATION
# ============================================================

N_PERM = 10000
rng = np.random.default_rng(20260813)

results = []

cancer_types = sorted(
    binary["cancer_type"].unique()
)

print("\n" + "=" * 90)
print("RUNNING PAN-CANCER PAIRWISE ANALYSIS")
print("=" * 90)

print(
    "Cancer types:",
    len(cancer_types)
)

print(
    "Gene pairs per cancer:",
    len(list(itertools.combinations(CLOCK16, 2)))
)

print(
    "Permutations per pair:",
    N_PERM
)

for cancer in cancer_types:

    sub = binary[
        binary["cancer_type"] == cancer
    ].copy()

    N = len(sub)

    if N < 20:
        print(
            f"SKIP {cancer}: N={N} (<20)"
        )
        continue

    X = sub[CLOCK16].to_numpy(
        dtype=np.int8
    )

    counts = X.sum(axis=0)

    for i, j in itertools.combinations(
        range(len(CLOCK16)), 2
    ):

        g1 = CLOCK16[i]
        g2 = CLOCK16[j]

        n1 = int(counts[i])
        n2 = int(counts[j])

        observed = int(
            (
                X[:, i]
                & X[:, j]
            ).sum()
        )

        expected_independence = (
            n1 * n2 / N
        )

        # ----------------------------------------------------
        # Fixed-margin null
        # ----------------------------------------------------

        null = np.empty(
            N_PERM,
            dtype=np.int16
        )

        population = np.arange(N)

        for k in range(N_PERM):

            p1 = rng.choice(
                population,
                size=n1,
                replace=False
            )

            p2 = rng.choice(
                population,
                size=n2,
                replace=False
            )

            null[k] = np.intersect1d(
                p1,
                p2,
                assume_unique=True
            ).size

        null_mean = float(
            null.mean()
        )

        null_sd = float(
            null.std(ddof=1)
        )

        oe = (
            observed / null_mean
            if null_mean > 0
            else np.nan
        )

        z = (
            (
                observed - null_mean
            ) / null_sd
            if null_sd > 0
            else np.nan
        )

        p = (
            (
                np.sum(
                    np.abs(
                        null - null_mean
                    )
                    >= abs(
                        observed - null_mean
                    )
                )
                + 1
            )
            / (N_PERM + 1)
        )

        results.append({
            "Cancer_Type": cancer,
            "N": N,
            "Gene1": g1,
            "Gene2": g2,
            "Gene1_Mutated": n1,
            "Gene2_Mutated": n2,
            "Observed_CoMutation": observed,
            "Expected_Independence": expected_independence,
            "Permutation_Expected": null_mean,
            "Null_SD": null_sd,
            "Observed_Expected_Ratio": oe,
            "Permutation_Z": z,
            "Permutation_P": p
        })

    print(
        f"{cancer:45s} N={N:5d}"
    )

results = pd.DataFrame(results)

# ============================================================
# 6. BENJAMINI-HOCHBERG FDR
# ============================================================

def bh_fdr(values):

    p = np.asarray(
        values,
        dtype=float
    )

    n = len(p)

    order = np.argsort(p)

    ranked = p[order]

    q = (
        ranked
        * n
        / np.arange(1, n + 1)
    )

    q = np.minimum.accumulate(
        q[::-1]
    )[::-1]

    output = np.empty(n)

    output[order] = q

    return np.clip(
        output,
        0,
        1
    )

results["FDR"] = bh_fdr(
    results["Permutation_P"]
)

results["Direction"] = np.where(
    results["Observed_Expected_Ratio"] > 1,
    "Enriched",
    np.where(
        results["Observed_Expected_Ratio"] < 1,
        "Depleted",
        "Neutral"
    )
)

# ============================================================
# 7. SAVE MAIN RESULTS
# ============================================================

main_project = os.path.join(
    RESULTS,
    "PAN_CANCER_CLOCK16_pairwise_burden_preserving.tsv"
)

main_download = os.path.join(
    OUTDIR,
    "PAN_CANCER_CLOCK16_pairwise_burden_preserving.tsv"
)

results.to_csv(
    main_project,
    sep="\t",
    index=False
)

results.to_csv(
    main_download,
    sep="\t",
    index=False
)

# ============================================================
# 8. SIGNIFICANT RESULTS
# ============================================================

sig = results[
    results["FDR"] < 0.05
].sort_values(
    ["FDR", "Permutation_P"]
)

sig_project = os.path.join(
    RESULTS,
    "PAN_CANCER_CLOCK16_pairwise_SIGNIFICANT.tsv"
)

sig_download = os.path.join(
    OUTDIR,
    "PAN_CANCER_CLOCK16_pairwise_SIGNIFICANT.tsv"
)

sig.to_csv(
    sig_project,
    sep="\t",
    index=False
)

sig.to_csv(
    sig_download,
    sep="\t",
    index=False
)

# ============================================================
# 9. CANCER-LEVEL SUMMARY
# ============================================================

summary = (
    results
    .groupby("Cancer_Type")
    .agg(
        N=("N", "first"),
        Tested_Pairs=("Gene1", "size"),
        Significant_Pairs=(
            "FDR",
            lambda x:
            int((x < 0.05).sum())
        ),
        Enriched_Pairs=(
            "Direction",
            lambda x:
            int((x == "Enriched").sum())
        ),
        Depleted_Pairs=(
            "Direction",
            lambda x:
            int((x == "Depleted").sum())
        ),
        Best_FDR=("FDR", "min"),
        Max_OE_Ratio=(
            "Observed_Expected_Ratio",
            "max"
        ),
        Min_OE_Ratio=(
            "Observed_Expected_Ratio",
            "min"
        )
    )
    .reset_index()
)

summary = summary.sort_values(
    [
        "Significant_Pairs",
        "Best_FDR"
    ],
    ascending=[
        False,
        True
    ]
)

summary_project = os.path.join(
    RESULTS,
    "PAN_CANCER_CLOCK16_pairwise_summary_by_cancer.tsv"
)

summary_download = os.path.join(
    OUTDIR,
    "PAN_CANCER_CLOCK16_pairwise_summary_by_cancer.tsv"
)

summary.to_csv(
    summary_project,
    sep="\t",
    index=False
)

summary.to_csv(
    summary_download,
    sep="\t",
    index=False
)

# ============================================================
# 10. GENE-PAIR REPRODUCIBILITY ACROSS CANCERS
# ============================================================

pair_summary = (
    results
    .groupby(
        ["Gene1", "Gene2"]
    )
    .agg(
        Cancer_Types=(
            "Cancer_Type",
            "nunique"
        ),
        Significant_Cancer_Count=(
            "FDR",
            lambda x:
            int((x < 0.05).sum())
        ),
        Enriched_Cancer_Count=(
            "Direction",
            lambda x:
            int((x == "Enriched").sum())
        ),
        Depleted_Cancer_Count=(
            "Direction",
            lambda x:
            int((x == "Depleted").sum())
        ),
        Best_FDR=(
            "FDR",
            "min"
        ),
        Median_OE_Ratio=(
            "Observed_Expected_Ratio",
            "median"
        ),
        Max_OE_Ratio=(
            "Observed_Expected_Ratio",
            "max"
        )
    )
    .reset_index()
)

pair_summary = pair_summary.sort_values(
    [
        "Significant_Cancer_Count",
        "Best_FDR"
    ],
    ascending=[
        False,
        True
    ]
)

pair_project = os.path.join(
    RESULTS,
    "PAN_CANCER_CLOCK16_pair_summary.tsv"
)

pair_download = os.path.join(
    OUTDIR,
    "PAN_CANCER_CLOCK16_pair_summary.tsv"
)

pair_summary.to_csv(
    pair_project,
    sep="\t",
    index=False
)

pair_summary.to_csv(
    pair_download,
    sep="\t",
    index=False
)

# ============================================================
# 11. FINAL QC
# ============================================================

print("\n" + "=" * 90)
print("FINAL PAN-CANCER CLOCK16 PAIRWISE QC")
print("=" * 90)

print(
    "Patients:",
    binary.index.nunique()
)

print(
    "Cancer types:",
    binary["cancer_type"].nunique()
)

print(
    "Total tested pairs:",
    len(results)
)

print(
    "Expected maximum if all 32 cancers:",
    32 * 120
)

print(
    "FDR < 0.05:",
    int(
        (results["FDR"] < 0.05).sum()
    )
)

print(
    "Enriched:",
    int(
        (results["Direction"] == "Enriched").sum()
    )
)

print(
    "Depleted:",
    int(
        (results["Direction"] == "Depleted").sum()
    )
)

print(
    "FDR range:",
    f"{results['FDR'].min():.6g}"
    f" - "
    f"{results['FDR'].max():.6g}"
)

print("\nTop 20 results:")

print(
    results
    .sort_values(
        ["FDR", "Permutation_P"]
    )
    .head(20)
    [
        [
            "Cancer_Type",
            "Gene1",
            "Gene2",
            "N",
            "Observed_CoMutation",
            "Permutation_Expected",
            "Observed_Expected_Ratio",
            "Permutation_Z",
            "Permutation_P",
            "FDR",
            "Direction"
        ]
    ]
    .to_string(index=False)
)

print("\nFiles saved:")
print(main_project)
print(main_download)
print(sig_project)
print(sig_download)
print(summary_project)
print(summary_download)
print(pair_project)
print(pair_download)

print("\n" + "=" * 90)
print("PAN-CANCER CLOCK16 PAIRWISE ANALYSIS COMPLETE")
print("=" * 90)

