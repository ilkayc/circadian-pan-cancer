import os
import itertools
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# ============================================================
# CLOCK16 FINAL PAN-CANCER PAIRWISE CO-MUTATION ANALYSIS
# ============================================================

BASE = os.path.abspath(".")

CLOCK16_GENES = [
    "ARNTL", "CLOCK", "CRY1", "CRY2",
    "CSNK1D", "CSNK1E", "NPAS2",
    "NR1D1", "NR1D2",
    "PER1", "PER2", "PER3",
    "RORA", "RORB", "RORC",
    "TIMELESS"
]

MIN_N = 20
N_PERM = 10000
SEED = 20260813

rng = np.random.default_rng(SEED)

EVENT_FILE = (
    "results/clock16/"
    "REAL_CLOCK16_patient_gene_event_counts.tsv"
)

META_FILE = (
    "results/clock16/"
    "PAN_CANCER_patient_mutational_burden.tsv"
)

OUTDIR = "results/clock16"
FIGDIR = os.path.expanduser(
    "~/Downloads/CLOCK16_FINAL/Pairwise"
)

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 1. LOAD
# ============================================================

print("=" * 80)
print("CLOCK16 FINAL PAN-CANCER PAIRWISE CO-MUTATION")
print("=" * 80)

events = pd.read_csv(EVENT_FILE, sep="\t")
meta = pd.read_csv(META_FILE, sep="\t")

print("\nEVENT INPUT")
print(EVENT_FILE)
print("Shape:", events.shape)
print("Columns:", list(events.columns))

print("\nMETADATA INPUT")
print(META_FILE)
print("Shape:", meta.shape)
print("Columns:", list(meta.columns))

# ============================================================
# 2. VALIDATE
# ============================================================

required_events = {"patient_id", "gene", "Mutation_Events"}
required_meta = {
    "patient_id",
    "cancer_type",
    "CLOCK16_mutated",
    "MC3_total_mutations",
    "MC3_CLOCK16_mutations"
}

missing_e = required_events - set(events.columns)
missing_m = required_meta - set(meta.columns)

if missing_e:
    raise ValueError(f"Missing event columns: {missing_e}")

if missing_m:
    raise ValueError(f"Missing metadata columns: {missing_m}")

events["patient_id"] = events["patient_id"].astype(str)
events["gene"] = events["gene"].astype(str)

meta["patient_id"] = meta["patient_id"].astype(str)
meta["cancer_type"] = meta["cancer_type"].astype(str)

events = events[events["gene"].isin(CLOCK16_GENES)].copy()

print("\nCLOCK16 genes detected:")
print(sorted(events["gene"].unique()))

print("CLOCK16 genes:", events["gene"].nunique())

# ============================================================
# 3. PATIENT-LEVEL CLOCK16 MATRIX
# ============================================================

# Any event in a gene = gene mutated in that patient.
gene_binary = (
    events.assign(mutated=1)
    .drop_duplicates(["patient_id", "gene"])
    .pivot_table(
        index="patient_id",
        columns="gene",
        values="mutated",
        fill_value=0
    )
)

gene_binary = gene_binary.reindex(
    columns=CLOCK16_GENES,
    fill_value=0
)

gene_binary = gene_binary.astype(int)

print("\nPATIENT-LEVEL CLOCK16 MATRIX")
print("Patients:", len(gene_binary))
print("Genes:", gene_binary.shape[1])

clock16_patients = gene_binary.index

# ============================================================
# 4. MAP CANCER TYPE
# ============================================================

meta_sub = meta[
    meta["patient_id"].isin(clock16_patients)
].copy()

# Keep one row per patient.
if meta_sub["patient_id"].duplicated().any():
    meta_sub = (
        meta_sub
        .sort_values("patient_id")
        .drop_duplicates("patient_id", keep="first")
    )

merged = gene_binary.reset_index().merge(
    meta_sub[
        [
            "patient_id",
            "cancer_type",
            "MC3_total_mutations",
            "MC3_CLOCK16_mutations"
        ]
    ],
    on="patient_id",
    how="left",
    validate="one_to_one"
)

print("\nPATIENT METADATA MAPPING")
print("CLOCK16 patients:", len(clock16_patients))
print(
    "Patients with cancer type:",
    merged["cancer_type"].notna().sum()
)
print(
    "Patients without cancer type:",
    merged["cancer_type"].isna().sum()
)

if merged["cancer_type"].isna().any():
    missing = merged.loc[
        merged["cancer_type"].isna(),
        "patient_id"
    ].tolist()

    print("WARNING: unmapped patients:")
    print(missing[:20])

merged = merged.dropna(subset=["cancer_type"]).copy()

# ============================================================
# 5. CHECK AGAINST MASTER COHORT
# ============================================================

print("\nCOHORT INTEGRITY")

master_clock16 = meta.loc[
    meta["CLOCK16_mutated"].astype(bool),
    ["patient_id", "cancer_type"]
].copy()

master_clock16["patient_id"] = (
    master_clock16["patient_id"].astype(str)
)

print(
    "Master CLOCK16-mutated patients:",
    master_clock16["patient_id"].nunique()
)

print(
    "Mapped CLOCK16 patients:",
    merged["patient_id"].nunique()
)

if (
    master_clock16["patient_id"].nunique()
    != merged["patient_id"].nunique()
):
    print(
        "WARNING: master and event-derived CLOCK16 cohorts differ."
    )

# ============================================================
# 6. CANCER DISTRIBUTION
# ============================================================

cancer_counts = (
    merged["cancer_type"]
    .value_counts()
    .sort_values(ascending=False)
)

print("\nCLOCK16 PATIENTS BY CANCER")

for cancer, n in cancer_counts.items():
    print(f"{cancer:45s} N={n:5d}")

# ============================================================
# 7. PAIRWISE PERMUTATION FUNCTION
# ============================================================

def pair_test(X, gene1, gene2, burden, n_perm=10000):
    """
    Permutation test preserving each patient's CLOCK16 gene burden.

    Null:
    For each patient, the number of mutated CLOCK16 genes is preserved.
    Gene labels are randomly reassigned within each patient.

    This tests whether a specific pair co-occurs more/less often
    than expected from patient-level CLOCK16 mutation burden.
    """

    a = X[gene1].to_numpy(dtype=np.int8)
    b = X[gene2].to_numpy(dtype=np.int8)

    observed = int(np.sum(a & b))

    k = burden.to_numpy(dtype=np.int16)

    # Expected co-mutation under the burden-preserving null:
    # probability that both selected genes are among k mutated genes.
    g = X.shape[1]

    expected_per_patient = np.where(
        k >= 2,
        (k * (k - 1)) / (g * (g - 1)),
        0
    )

    expected = float(expected_per_patient.sum())

    # Observed / expected.
    ratio = (
        observed / expected
        if expected > 0
        else np.nan
    )

    # --------------------------------------------------------
    # Permutation
    # --------------------------------------------------------

    perm_counts = np.empty(n_perm, dtype=np.int16)

    for i in range(n_perm):

        count = 0

        for kk in np.unique(k):

            idx = np.where(k == kk)[0]

            if kk < 2 or len(idx) == 0:
                continue

            # Number of mutated genes per patient is preserved.
            #
            # For each patient, randomly choose kk genes.
            # We only need to track whether gene1 and gene2
            # are simultaneously selected.
            #
            # Vectorized across patients in this burden stratum.

            n = len(idx)

            # Random scores for the two genes.
            # Approximate exact within-patient randomization
            # using random ranks among all 16 genes.

            # Generate n x G random matrix.
            r = rng.random((n, g))

            # Indices of selected genes.
            selected = np.argpartition(
                r,
                kth=kk - 1,
                axis=1
            )[:, :kk]

            has1 = np.any(
                selected == list(X.columns).index(gene1),
                axis=1
            )

            has2 = np.any(
                selected == list(X.columns).index(gene2),
                axis=1
            )

            count += int(np.sum(has1 & has2))

        perm_counts[i] = count

    # Two-sided empirical permutation P.
    p = (
        (np.sum(
            np.abs(perm_counts - expected)
            >= np.abs(observed - expected)
        ) + 1)
        / (n_perm + 1)
    )

    sd = float(np.std(perm_counts, ddof=1))

    z = (
        (observed - expected) / sd
        if sd > 0
        else np.nan
    )

    if observed > expected:
        direction = "Enriched"
    elif observed < expected:
        direction = "Depleted"
    else:
        direction = "Neutral"

    return {
        "Observed_CoMutation": observed,
        "Permutation_Expected": expected,
        "Observed_Expected_Ratio": ratio,
        "Permutation_SD": sd,
        "Permutation_Z": z,
        "Permutation_P": p,
        "Direction": direction
    }

# ============================================================
# 8. RUN ALL CANCERS
# ============================================================

pairs = list(itertools.combinations(CLOCK16_GENES, 2))

print("\n" + "=" * 80)
print("PAIRWISE ANALYSIS")
print("=" * 80)

print("Genes:", len(CLOCK16_GENES))
print("Pairs per cancer:", len(pairs))
print("Minimum N:", MIN_N)
print("Permutations per pair:", N_PERM)

results = []

for cancer, sub in merged.groupby("cancer_type"):

    sub = sub.copy()
    n = len(sub)

    if n < MIN_N:
        print(
            f"SKIP {cancer}: N={n} (<{MIN_N})"
        )
        continue

    print(
        f"RUN  {cancer}: N={n}"
    )

    X = sub[CLOCK16_GENES].astype(int)

    burden = X.sum(axis=1)

    # Need variation in burden.
    if burden.sum() == 0:
        print("  SKIP: no CLOCK16 burden")
        continue

    for gene1, gene2 in pairs:

        res = pair_test(
            X,
            gene1,
            gene2,
            burden,
            n_perm=N_PERM
        )

        res.update({
            "Cancer_Type": cancer,
            "Gene1": gene1,
            "Gene2": gene2,
            "N": n,
            "CLOCK16_Mutated_Patients": int(len(sub)),
            "Mean_CLOCK16_Genes_Per_Patient":
                float(burden.mean())
        })

        results.append(res)

# ============================================================
# 9. RESULTS
# ============================================================

df = pd.DataFrame(results)

if df.empty:
    raise RuntimeError(
        "No cancer type passed the minimum N threshold."
    )

df = df[
    [
        "Cancer_Type",
        "Gene1",
        "Gene2",
        "N",
        "CLOCK16_Mutated_Patients",
        "Mean_CLOCK16_Genes_Per_Patient",
        "Observed_CoMutation",
        "Permutation_Expected",
        "Observed_Expected_Ratio",
        "Permutation_SD",
        "Permutation_Z",
        "Permutation_P",
        "Direction"
    ]
].copy()

# ============================================================
# 10. FDR
# ============================================================

# Cancer-specific FDR:
df["FDR_within_cancer"] = np.nan

for cancer, idx in df.groupby("Cancer_Type").groups.items():

    p = df.loc[idx, "Permutation_P"].to_numpy()

    df.loc[idx, "FDR_within_cancer"] = (
        multipletests(
            p,
            method="fdr_bh"
        )[1]
    )

# Global pan-cancer FDR:
df["FDR_global"] = (
    multipletests(
        df["Permutation_P"].to_numpy(),
        method="fdr_bh"
    )[1]
)

# ============================================================
# 11. ROBUST CLASSIFICATION
# ============================================================

df["Robust_Significant"] = (
    (df["FDR_global"] < 0.05)
    &
    (df["FDR_within_cancer"] < 0.05)
)

# ============================================================
# 12. SAVE RAW RESULTS
# ============================================================

raw_file = (
    f"{OUTDIR}/"
    "PAN_CANCER_CLOCK16_pairwise_FINAL.tsv"
)

df.to_csv(
    raw_file,
    sep="\t",
    index=False
)

raw_fig_file = (
    f"{FIGDIR}/"
    "PAN_CANCER_CLOCK16_pairwise_FINAL.tsv"
)

df.to_csv(
    raw_fig_file,
    sep="\t",
    index=False
)

# ============================================================
# 13. SIGNIFICANT RESULTS
# ============================================================

sig = df[
    df["Robust_Significant"]
].copy()

sig = sig.sort_values(
    [
        "FDR_global",
        "FDR_within_cancer",
        "Permutation_P"
    ]
)

sig_file = (
    f"{OUTDIR}/"
    "PAN_CANCER_CLOCK16_pairwise_ROBUST.tsv"
)

sig.to_csv(
    sig_file,
    sep="\t",
    index=False
)

# ============================================================
# 14. PAIR-LEVEL RECURRENCE ACROSS CANCERS

# ============================================================

print("\n" + "=" * 80)
print("PAIR-LEVEL RECURRENCE ACROSS CANCERS")
print("=" * 80)

pair_summary = (
    df.groupby(["Gene1", "Gene2"], as_index=False)
      .agg(
          Cancers_Tested=("Cancer_Type", "nunique"),
          FDR05_Cancers=(
              "FDR_within_cancer",
              lambda x: int((x < 0.05).sum())
          ),
          Global_FDR05_Cancers=(
              "FDR_global",
              lambda x: int((x < 0.05).sum())
          ),
          Robust_Cancers=(
              "Robust_Significant",
              "sum"
          ),
          Enriched_Cancers=(
              "Direction",
              lambda x: int((x == "Enriched").sum())
          ),
          Depleted_Cancers=(
              "Direction",
              lambda x: int((x == "Depleted").sum())
          )
      )
)

# Recurrent associations across at least 2 cancers

pair_summary["Recurrent_FDR05"] = (
    pair_summary["FDR05_Cancers"] >= 2
)

pair_summary["Recurrent_Global_FDR05"] = (
    pair_summary["Global_FDR05_Cancers"] >= 2
)

pair_summary["Recurrent_Robust"] = (
    pair_summary["Robust_Cancers"] >= 2
)

# Sort strongest recurrent pairs first

pair_summary = pair_summary.sort_values(
    [
        "Recurrent_Robust",
        "Recurrent_FDR05",
        "Robust_Cancers",
        "FDR05_Cancers",
        "Enriched_Cancers"
    ],
    ascending=[False, False, False, False, False]
).reset_index(drop=True)

# Save pair-level summary

pair_summary_file = os.path.join(
    OUT_DIR,
    "PAN_CANCER_CLOCK16_pair_summary.tsv"
)

pair_summary.to_csv(
    pair_summary_file,
    sep="\t",
    index=False
)

print("\nTop recurrent gene pairs:")
print(
    pair_summary.head(20).to_string(index=False)
)

print("\nSaved:")
print(pair_summary_file)

# ============================================================
# ROBUST RECURRENT PAIRS
# ============================================================

robust_pairs = pair_summary[
    pair_summary["Recurrent_Robust"]
].copy()

robust_pair_file = os.path.join(
    OUT_DIR,
    "PAN_CANCER_CLOCK16_pair_summary_ROBUST.tsv"
)

robust_pairs.to_csv(
    robust_pair_file,
    sep="\t",
    index=False
)

print("\nRobust recurrent gene pairs:",
      len(robust_pairs))

print("Saved:")
print(robust_pair_file)

# 15. CANCER SUMMARY
# ============================================================

cancer_summary = (
    df.groupby("Cancer_Type")
    .agg(
        N=("N", "first"),
        Pairs_Tested=("Gene1", "size"),
        FDR05_Pairs_Within_Cancer=(
            "FDR_within_cancer",
            lambda x: int((x < 0.05).sum())
        ),
        Global_FDR05_Pairs=(
            "FDR_global",
            lambda x: int((x < 0.05).sum())
        ),
        Robust_Pairs=(
            "Robust_Significant",
            "sum"
        ),
        Enriched_Pairs=(
            "Direction",
            lambda x: int((x == "Enriched").sum())
        ),
        Depleted_Pairs=(
            "Direction",
            lambda x: int((x == "Depleted").sum())
        )
    )
    .reset_index()
)

cancer_summary = cancer_summary.sort_values(
    "N",
    ascending=False
)

cancer_file = (
    f"{OUTDIR}/"
    "PAN_CANCER_CLOCK16_pairwise_cancer_summary_FINAL.tsv"
)

cancer_summary.to_csv(
    cancer_file,
    sep="\t",
    index=False
)

# ============================================================
# 16. QC
# ============================================================

qc = []

qc.append(
    ("Master_CLOCK16_patients", master_clock16["patient_id"].nunique())
)

qc.append(
    ("Event_CLOCK16_patients", len(clock16_patients))
)

qc.append(
    ("Mapped_CLOCK16_patients", merged["patient_id"].nunique())
)

qc.append(
    ("Cancer_types_in_master_metadata",
     meta["cancer_type"].nunique())
)

qc.append(
    ("Cancer_types_with_CLOCK16",
     merged["cancer_type"].nunique())
)

qc.append(
    ("Cancer_types_tested",
     df["Cancer_Type"].nunique())
)

qc.append(
    ("Pairs_per_cancer",
     len(pairs))
)

qc.append(
    ("Total_pairwise_tests",
     len(df))
)

qc.append(
    ("Global_FDR05",
     int((df["FDR_global"] < 0.05).sum()))
)

qc.append(
    ("Within_cancer_FDR05",
     int((df["FDR_within_cancer"] < 0.05).sum()))
)

qc.append(
    ("Robust_significant",
     int(df["Robust_Significant"].sum()))
)

qc_df = pd.DataFrame(
    qc,
    columns=["Metric", "Value"]
)

qc_file = (
    f"{OUTDIR}/"
    "PAN_CANCER_CLOCK16_pairwise_FINAL_QC.tsv"
)

qc_df.to_csv(
    qc_file,
    sep="\t",
    index=False
)

# ============================================================
# 17. FINAL REPORT
# ============================================================

print("\n" + "=" * 80)
print("FINAL PAN-CANCER PAIRWISE RESULT")
print("=" * 80)

print(
    "Master CLOCK16 patients:",
    master_clock16["patient_id"].nunique()
)

print(
    "Event-derived CLOCK16 patients:",
    len(clock16_patients)
)

print(
    "Mapped CLOCK16 patients:",
    merged["patient_id"].nunique()
)

print(
    "Cancer types in master metadata:",
    meta["cancer_type"].nunique()
)

print(
    "Cancer types containing CLOCK16:",
    merged["cancer_type"].nunique()
)

print(
    "Cancer types tested:",
    df["Cancer_Type"].nunique()
)

print(
    "Pairs per cancer:",
    len(pairs)
)

print(
    "Total pairwise tests:",
    len(df)
)

print(
    "Within-cancer FDR < 0.05:",
    int((df["FDR_within_cancer"] < 0.05).sum())
)

print(
    "Global FDR < 0.05:",
    int((df["FDR_global"] < 0.05).sum())
)

print(
    "Robust significant:",
    int(df["Robust_Significant"].sum())
)

print("\nTOP ROBUST RESULTS")

if len(sig) > 0:
    print(
        sig[
            [
                "Cancer_Type",
                "Gene1",
                "Gene2",
                "N",
                "Observed_CoMutation",
                "Permutation_Expected",
                "Observed_Expected_Ratio",
                "Permutation_P",
                "FDR_within_cancer",
                "FDR_global",
                "Direction"
            ]
        ]
        .head(30)
        .to_string(index=False)
    )
else:
    print("No robust pairwise associations survived both FDR criteria.")

print("\nTOP RECURRENT PAIRS")

print(
    pair_summary.head(30).to_string(index=False)
)

print("\nFILES SAVED:")
print(raw_file)
print(raw_fig_file)
print(sig_file)
print(pair_file)
print(cancer_file)
print(qc_file)

print("\nANALYSIS COMPLETE.")
