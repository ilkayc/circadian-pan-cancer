import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

BASE = Path("results/clock16")

OUT = (
    Path.home()
    / "Downloads"
    / "CLOCK16_FINAL"
    / "Figures"
)

OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# NATURE-STYLE GLOBAL SETTINGS
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04
})

# ============================================================
# INPUT FILES
# ============================================================

PREV = BASE / "CLOCK16_final_cancer_prevalence.tsv"

CANCER = BASE / "REAL_CLOCK16_cancer_association_any_v2.tsv"

GENE = BASE / "CLOCK16_gene_level_pooled_FIRTH_LR.tsv"

FUNCTIONAL = (
    BASE /
    "CLOCK16_functional_CORRECTED_pooled_FIRTH_LR.tsv"
)

# ============================================================
# LOAD DATA
# ============================================================

prev = pd.read_csv(PREV, sep="\t")
cancer = pd.read_csv(CANCER, sep="\t")
gene = pd.read_csv(GENE, sep="\t")
functional = pd.read_csv(FUNCTIONAL, sep="\t")

# ============================================================
# QC
# ============================================================

assert len(prev) == 32
assert len(cancer) == 32
assert len(gene) == 16
assert len(functional) == 3

assert (
    prev["Prevalence"].between(0, 1).all()
)

assert (
    cancer["odds_ratio"] > 0
).all()

assert (
    gene["OR_per_log1p_non_CLOCK16_mutation"] > 0
).all()

assert (
    functional["OR_per_log1p_non_CLOCK16_mutation"] > 0
).all()

# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(
    figsize=(11.8, 10.5),
    facecolor="white"
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.08,
    right=0.98,
    bottom=0.07,
    top=0.97,
    wspace=0.30,
    hspace=0.38
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        direction="out",
        length=3,
        width=0.7
    )

def panel_label(ax, label):
    ax.text(
        -0.14,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="bottom",
        ha="left"
    )

def add_null_line(ax):
    ax.axvline(
        1,
        linestyle="--",
        linewidth=0.8,
        alpha=0.65
    )

# ============================================================
# PANEL A — PREVALENCE
# ============================================================

a = prev.sort_values(
    "Prevalence",
    ascending=True
).reset_index(drop=True)

y = np.arange(len(a))

axA.hlines(
    y,
    a["CI95_Low"],
    a["CI95_High"],
    linewidth=1.2
)

axA.scatter(
    a["Prevalence"],
    y,
    s=24,
    zorder=3
)

axA.set_yticks(y)
axA.set_yticklabels(
    a["Cancer_Type"],
    fontsize=7
)

axA.set_xlabel(
    "CLOCK16-mutant prevalence"
)

axA.set_xlim(
    0,
    max(a["CI95_High"]) * 1.08
)

axA.set_xticks(
    [0, 0.1, 0.2, 0.3]
)

axA.set_xticklabels(
    ["0%", "10%", "20%", "30%"]
)

axA.set_title(
    "Pan-cancer prevalence",
    loc="left",
    fontweight="bold"
)

panel_label(axA, "A")
clean_axis(axA)

# ============================================================
# PANEL B — CANCER ASSOCIATION
# ============================================================

b = cancer.sort_values(
    "odds_ratio",
    ascending=True
).reset_index(drop=True)

y = np.arange(len(b))

OR = b["odds_ratio"].to_numpy()
LOW = b["lower_95CI"].to_numpy()
HIGH = b["upper_95CI"].to_numpy()

xerr = np.vstack([
    OR - LOW,
    HIGH - OR
])

sig = b["FDR"] < 0.05

axB.errorbar(
    OR,
    y,
    xerr=xerr,
    fmt="o",
    markersize=4.5,
    capsize=2.5,
    linewidth=1.0,
    zorder=3
)

axB.scatter(
    OR[sig],
    y[sig],
    s=22,
    zorder=4
)

add_null_line(axB)

axB.set_xscale("log")

axB.set_yticks(y)
axB.set_yticklabels(
    b["cancer_type"],
    fontsize=6.8
)

axB.set_xlabel(
    "Odds ratio (95% CI)"
)

axB.set_title(
    "Cancer-type association",
    loc="left",
    fontweight="bold"
)

panel_label(axB, "B")
clean_axis(axB)

# ============================================================
# PANEL C — GENE ASSOCIATION
# ============================================================

c = gene.sort_values(
    "OR_per_log1p_non_CLOCK16_mutation",
    ascending=True
).reset_index(drop=True)

y = np.arange(len(c))

OR = c[
    "OR_per_log1p_non_CLOCK16_mutation"
].to_numpy()

LOW = c["CI95_Low"].to_numpy()
HIGH = c["CI95_High"].to_numpy()

xerr = np.vstack([
    OR - LOW,
    HIGH - OR
])

axC.errorbar(
    OR,
    y,
    xerr=xerr,
    fmt="o",
    markersize=5,
    capsize=2.5,
    linewidth=1.0
)

axC.axvline(
    1,
    linestyle="--",
    linewidth=0.8,
    alpha=0.65
)

axC.set_yticks(y)
axC.set_yticklabels(
    c["Gene"],
    fontsize=8
)

axC.set_xlabel(
    "Odds ratio per log1p non-CLOCK16 mutation burden (95% CI)"
)

axC.set_title(
    "Gene-level association",
    loc="left",
    fontweight="bold"
)

panel_label(axC, "C")
clean_axis(axC)

# ============================================================
# PANEL D — FUNCTIONAL ASSOCIATION
# ============================================================

d_order = [
    "Nonsynonymous",
    "Damaging",
    "Strong_Damaging"
]

d = functional.copy()

d["Endpoint"] = pd.Categorical(
    d["Endpoint"],
    categories=d_order,
    ordered=True
)

d = d.sort_values(
    "Endpoint"
).reset_index(drop=True)

y = np.arange(len(d))

OR = d[
    "OR_per_log1p_non_CLOCK16_mutation"
].to_numpy()

LOW = d["CI95_Low"].to_numpy()
HIGH = d["CI95_High"].to_numpy()

xerr = np.vstack([
    OR - LOW,
    HIGH - OR
])

axD.errorbar(
    OR,
    y,
    xerr=xerr,
    fmt="o",
    markersize=5.5,
    capsize=3,
    linewidth=1.1
)

axD.axvline(
    1,
    linestyle="--",
    linewidth=0.8,
    alpha=0.65
)

axD.set_yticks(y)

axD.set_yticklabels(
    [
        "Nonsynonymous",
        "Damaging",
        "Strong damaging"
    ],
    fontsize=8.5
)

axD.set_xlabel(
    "Odds ratio per log1p non-CLOCK16 mutation burden (95% CI)"
)

axD.set_title(
    "Functional mutation classes",
    loc="left",
    fontweight="bold"
)

panel_label(axD, "D")
clean_axis(axD)

# ============================================================
# CONSISTENT X AXES WHERE APPROPRIATE
# ============================================================

axC.set_xlim(2.0, 4.1)
axD.set_xlim(2.5, 4.0)

# ============================================================
# SAVE COMBINED FIGURE
# ============================================================

png = OUT / "Figure1_Nature_CLOCK16_pan_cancer_ABCD.png"
pdf = OUT / "Figure1_Nature_CLOCK16_pan_cancer_ABCD.pdf"

fig.savefig(
    png,
    dpi=600,
    facecolor="white"
)

fig.savefig(
    pdf,
    facecolor="white"
)

plt.close(fig)

# ============================================================
# ALSO SAVE INDIVIDUAL HIGH-QUALITY PANELS
# ============================================================

# Re-open individual panels from the same data is unnecessary;
# the combined figure is the primary publication figure.
# Existing individual figures are retained.

print("=" * 80)
print("NATURE-STYLE FIGURE 1 — FINAL")
print("=" * 80)

print("\nPanel A:")
print(f"32 cancer types; prevalence range "
      f"{a['Prevalence'].min()*100:.2f}%–"
      f"{a['Prevalence'].max()*100:.2f}%")

print("\nPanel B:")
print(f"32 cancer types")
print(f"FDR < 0.05: {(b['FDR'] < 0.05).sum()}")

print("\nPanel C:")
print(f"16 genes")
print(f"FDR < 0.05: {(c['FDR'] < 0.05).sum()}")

print("\nPanel D:")
print(f"3 functional endpoints")
print(f"FDR < 0.05: {(d['FDR'] < 0.05).sum()}")

print("\nOutput PNG:")
print(png)

print("\nOutput PDF:")
print(pdf)

print("\n" + "=" * 80)
print("NATURE-STYLE FIGURE 1 CREATED")
print("=" * 80)
