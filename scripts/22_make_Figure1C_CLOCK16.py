import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

INPUT = Path(
    "results/clock16/CLOCK16_gene_level_pooled_FIRTH_LR.tsv"
)

PROJECT_PNG = Path(
    "results/clock16/Figure1C_CLOCK16_gene_Firth.png"
)

PROJECT_PDF = Path(
    "results/clock16/Figure1C_CLOCK16_gene_Firth.pdf"
)

DOWNLOADS_DIR = (
    Path.home()
    / "Downloads"
    / "CLOCK16_FINAL"
    / "Figures"
)

DOWNLOADS_PNG = (
    DOWNLOADS_DIR
    / "Figure1C_CLOCK16_gene_Firth.png"
)

DOWNLOADS_PDF = (
    DOWNLOADS_DIR
    / "Figure1C_CLOCK16_gene_Firth.pdf"
)

# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INPUT, sep="\t")

# ============================================================
# QC
# ============================================================

assert len(df) == 16, f"Expected 16 genes, found {len(df)}"
assert df["Gene"].nunique() == 16

required = [
    "Gene",
    "OR_per_log1p_non_CLOCK16_mutation",
    "CI95_Low",
    "CI95_High",
    "Firth_LR_P",
    "FDR"
]

for col in required:
    assert col in df.columns, f"Missing column: {col}"

assert (df["OR_per_log1p_non_CLOCK16_mutation"] > 0).all()
assert (df["CI95_Low"] > 0).all()
assert (df["CI95_High"] > 0).all()
assert (
    df["CI95_Low"] <= df["CI95_High"]
).all()

# ============================================================
# ORDER
# Strongest association first
# ============================================================

df = df.sort_values(
    "OR_per_log1p_non_CLOCK16_mutation",
    ascending=False
).reset_index(drop=True)

OR = df[
    "OR_per_log1p_non_CLOCK16_mutation"
].to_numpy()

LOW = df["CI95_Low"].to_numpy()
HIGH = df["CI95_High"].to_numpy()

FDR = df["FDR"].to_numpy()

significant = FDR < 0.05

y = np.arange(len(df))

xerr = np.vstack([
    OR - LOW,
    HIGH - OR
])

# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(8.5, 7.5)
)

# CI + point
ax.errorbar(
    OR,
    y,
    xerr=xerr,
    fmt="o",
    markersize=5,
    capsize=3,
    linewidth=1.2
)

# Highlight FDR-significant genes
ax.scatter(
    OR[significant],
    y[significant],
    s=34
)

# Null
ax.axvline(
    1,
    linestyle="--",
    linewidth=1
)

# Log scale
ax.set_xscale("log")

# Gene labels
ax.set_yticks(y)
ax.set_yticklabels(
    df["Gene"],
    fontsize=10
)

ax.invert_yaxis()

ax.set_xlabel(
    "Odds ratio per log1p non-CLOCK16 mutation burden (95% CI)",
    fontsize=10.5
)

ax.set_ylabel(
    "CLOCK gene",
    fontsize=11
)

ax.set_title(
    "Gene-level burden-adjusted Firth associations",
    fontsize=13,
    fontweight="bold",
    pad=12
)

# Grid
ax.grid(
    axis="x",
    linestyle=":",
    linewidth=0.6,
    alpha=0.5
)

# Clean frame
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Panel label
ax.text(
    -0.12,
    1.02,
    "C",
    transform=ax.transAxes,
    fontsize=16,
    fontweight="bold",
    va="bottom"
)

# ============================================================
# SAVE
# ============================================================

DOWNLOADS_DIR.mkdir(
    parents=True,
    exist_ok=True
)

fig.tight_layout()

fig.savefig(
    PROJECT_PNG,
    dpi=600,
    bbox_inches="tight"
)

fig.savefig(
    PROJECT_PDF,
    bbox_inches="tight"
)

fig.savefig(
    DOWNLOADS_PNG,
    dpi=600,
    bbox_inches="tight"
)

fig.savefig(
    DOWNLOADS_PDF,
    bbox_inches="tight"
)

plt.close(fig)

# ============================================================
# QC
# ============================================================

print("=" * 80)
print("FIGURE 1C — FINAL")
print("=" * 80)

print(f"Genes:              {len(df)}")
print(f"FDR < 0.05:         {significant.sum()}")
print(f"FDR >= 0.05:        {(~significant).sum()}")

print("\nOR range:")
print(
    f"{OR.min():.6f} - {OR.max():.6f}"
)

print("\nGene order:")
print(
    df[
        [
            "Gene",
            "OR_per_log1p_non_CLOCK16_mutation",
            "CI95_Low",
            "CI95_High",
            "FDR"
        ]
    ].to_string(index=False)
)

print("\nProject:")
print(PROJECT_PNG)
print(PROJECT_PDF)

print("\nDownloads:")
print(DOWNLOADS_PNG)
print(DOWNLOADS_PDF)

print("\n" + "=" * 80)
print("FIGURE 1C CREATED SUCCESSFULLY")
print("=" * 80)
