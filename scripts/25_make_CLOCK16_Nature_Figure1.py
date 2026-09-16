import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ============================================================
# CLOCK16 — NATURE-STYLE FIGURE 1
# Panels A-D generated directly from source data
# ============================================================

PROJECT = "results/clock16"
OUT = os.path.expanduser("~/Downloads/CLOCK16_FINAL/Figures")
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------
# SOURCE FILES
# ------------------------------------------------------------

FILE_A = f"{PROJECT}/REAL_clock16_cancer_prevalence.tsv"
FILE_B = f"{PROJECT}/REAL_CLOCK16_cancer_association_any_v2.tsv"
FILE_C = f"{PROJECT}/CLOCK16_gene_level_pooled_FIRTH_LR.tsv"
FILE_D = f"{PROJECT}/CLOCK16_functional_CORRECTED_pooled_FIRTH_LR.tsv"

# ------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------

PNG = f"{OUT}/Figure1_CLOCK16_Nature.png"
PDF = f"{OUT}/Figure1_CLOCK16_Nature.pdf"
SVG = f"{OUT}/Figure1_CLOCK16_Nature.svg"

# ------------------------------------------------------------
# GLOBAL STYLE
# ------------------------------------------------------------

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

# ============================================================
# LOAD
# ============================================================

A = pd.read_csv(FILE_A, sep="\t")
B = pd.read_csv(FILE_B, sep="\t")
C = pd.read_csv(FILE_C, sep="\t")
D = pd.read_csv(FILE_D, sep="\t")

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------

required_A = [
    "cancer_type",
    "Total_Primary_Patients",
    "CLOCK16_Mutated_Patients",
    "Prevalence_Percent",
    "Lower_95CI",
    "Upper_95CI"
]

required_B = [
    "phenotype",
    "cancer_type",
    "total_primary_patients",
    "mutated_patients",
    "prevalence_percent",
    "odds_ratio",
    "lower_95CI",
    "upper_95CI",
    "fisher_p",
    "FDR"
]

required_C = [
    "Gene",
    "OR_per_log1p_non_CLOCK16_mutation",
    "CI95_Low",
    "CI95_High",
    "FDR"
]

required_D = [
    "Endpoint",
    "Cases",
    "Controls",
    "OR_per_log1p_non_CLOCK16_mutation",
    "CI95_Low",
    "CI95_High",
    "Firth_LR_P",
    "FDR"
]

for name, df, req in [
    ("A", A, required_A),
    ("B", B, required_B),
    ("C", C, required_C),
    ("D", D, required_D),
]:
    missing = [x for x in req if x not in df.columns]
    if missing:
        raise ValueError(
            f"Figure 1{name} missing columns: {missing}\n"
            f"Available: {df.columns.tolist()}"
        )

# ============================================================
# CLEAN
# ============================================================

# A
A = A.copy()
A["cancer_type"] = A["cancer_type"].astype(str)
A = A.sort_values("Prevalence_Percent", ascending=True)

# B — IMPORTANT: use the Any_CLOCK16 rows only
B = B[B["phenotype"].astype(str) == "Any_CLOCK16"].copy()

B["cancer_type"] = B["cancer_type"].astype(str)
B["odds_ratio"] = pd.to_numeric(B["odds_ratio"])
B["lower_95CI"] = pd.to_numeric(B["lower_95CI"])
B["upper_95CI"] = pd.to_numeric(B["upper_95CI"])
B["FDR"] = pd.to_numeric(B["FDR"])

B = B.sort_values("odds_ratio", ascending=True)

# C
C = C.copy()
C["OR_per_log1p_non_CLOCK16_mutation"] = pd.to_numeric(
    C["OR_per_log1p_non_CLOCK16_mutation"]
)
C["CI95_Low"] = pd.to_numeric(C["CI95_Low"])
C["CI95_High"] = pd.to_numeric(C["CI95_High"])
C["FDR"] = pd.to_numeric(C["FDR"])

C = C.sort_values(
    "OR_per_log1p_non_CLOCK16_mutation",
    ascending=True
)

# D
D = D.copy()

D["OR_per_log1p_non_CLOCK16_mutation"] = pd.to_numeric(
    D["OR_per_log1p_non_CLOCK16_mutation"]
)
D["CI95_Low"] = pd.to_numeric(D["CI95_Low"])
D["CI95_High"] = pd.to_numeric(D["CI95_High"])

# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(figsize=(12.0, 9.2))

gs = fig.add_gridspec(
    2, 2,
    width_ratios=[1.15, 1.0],
    height_ratios=[1.0, 0.90],
    wspace=0.32,
    hspace=0.34
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# ============================================================
# PANEL A — PREVALENCE
# ============================================================

y = np.arange(len(A))

x = A["Prevalence_Percent"].values
lo = A["Lower_95CI"].values
hi = A["Upper_95CI"].values

axA.errorbar(
    x,
    y,
    xerr=[x - lo, hi - x],
    fmt="o",
    markersize=3.6,
    linewidth=0.8,
    capsize=2.2
)

axA.set_yticks(y)
axA.set_yticklabels(A["cancer_type"])
axA.set_xlabel("CLOCK16-mutated patients (%)")
axA.set_title("Pan-cancer prevalence", loc="left", fontweight="bold")

axA.grid(axis="x", linewidth=0.35, alpha=0.25)
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

# ============================================================
# PANEL B — CANCER ASSOCIATION
# ============================================================

y = np.arange(len(B))

OR = B["odds_ratio"].values
lo = B["lower_95CI"].values
hi = B["upper_95CI"].values

sig = B["FDR"].values < 0.05

axB.axvline(
    1,
    linewidth=0.8,
    linestyle="--"
)

for i in range(len(B)):
    axB.plot(
        [lo[i], hi[i]],
        [y[i], y[i]],
        linewidth=1.0
    )

    axB.plot(
        OR[i],
        y[i],
        "o",
        markersize=4.0
    )

axB.set_yticks(y)
axB.set_yticklabels(B["cancer_type"])
axB.set_xscale("log")
axB.set_xlabel("Odds ratio")
axB.set_title("Cancer-type association", loc="left", fontweight="bold")

axB.grid(axis="x", linewidth=0.35, alpha=0.25)
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# ============================================================
# PANEL C — GENE-LEVEL FIRTH
# ============================================================

y = np.arange(len(C))

OR = C["OR_per_log1p_non_CLOCK16_mutation"].values
lo = C["CI95_Low"].values
hi = C["CI95_High"].values

axC.axvline(
    1,
    linewidth=0.8,
    linestyle="--"
)

for i in range(len(C)):
    axC.plot(
        [lo[i], hi[i]],
        [y[i], y[i]],
        linewidth=1.0
    )
    axC.plot(
        OR[i],
        y[i],
        "o",
        markersize=4.0
    )

axC.set_yticks(y)
axC.set_yticklabels(C["Gene"])
axC.set_xscale("log")
axC.set_xlabel(
    "Odds ratio per log1p non-CLOCK16 mutation burden"
)
axC.set_title(
    "Clock-gene associations",
    loc="left",
    fontweight="bold"
)

axC.grid(axis="x", linewidth=0.35, alpha=0.25)
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)

# ============================================================
# PANEL D — FUNCTIONAL CONSEQUENCE
# ============================================================

y = np.arange(len(D))

OR = D["OR_per_log1p_non_CLOCK16_mutation"].values
lo = D["CI95_Low"].values
hi = D["CI95_High"].values

axD.axvline(
    1,
    linewidth=0.8,
    linestyle="--"
)

for i in range(len(D)):
    axD.plot(
        [lo[i], hi[i]],
        [y[i], y[i]],
        linewidth=1.2
    )
    axD.plot(
        OR[i],
        y[i],
        "o",
        markersize=5
    )

axD.set_yticks(y)
axD.set_yticklabels(D["Endpoint"])
axD.set_xscale("log")
axD.set_xlabel(
    "Odds ratio per log1p non-CLOCK16 mutation burden"
)
axD.set_title(
    "Functional mutation classes",
    loc="left",
    fontweight="bold"
)

axD.grid(axis="x", linewidth=0.35, alpha=0.25)
axD.spines["top"].set_visible(False)
axD.spines["right"].set_visible(False)

# ============================================================
# PANEL LABELS
# ============================================================

for ax, label in zip(
    [axA, axB, axC, axD],
    ["a", "b", "c", "d"]
):
    ax.text(
        -0.12,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="right"
    )

# ============================================================
# FINAL LAYOUT
# ============================================================

fig.subplots_adjust(
    left=0.24,
    right=0.97,
    bottom=0.09,
    top=0.96
)

# High-resolution raster
fig.savefig(
    PNG,
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)

# Vector PDF
fig.savefig(
    PDF,
    bbox_inches="tight",
    facecolor="white"
)

# Vector SVG
fig.savefig(
    SVG,
    bbox_inches="tight",
    facecolor="white"
)

plt.close(fig)

# ============================================================
# QC
# ============================================================

print("=" * 80)
print("CLOCK16 NATURE-STYLE FIGURE 1")
print("=" * 80)

print(f"Panel A: {len(A)} cancer types")
print(f"Panel B: {len(B)} cancer types")
print(f"Panel C: {len(C)} genes")
print(f"Panel D: {len(D)} functional endpoints")

print("\nOutputs:")
print(PNG)
print(PDF)
print(SVG)

print("\nFigure 1 created successfully.")
print("=" * 80)
