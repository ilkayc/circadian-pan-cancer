import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

INPUT = Path(
    "results/clock16/CLOCK16_functional_CORRECTED_pooled_FIRTH_LR.tsv"
)

PROJECT_PNG = Path(
    "results/clock16/Figure1D_CLOCK16_functional_Firth.png"
)

PROJECT_PDF = Path(
    "results/clock16/Figure1D_CLOCK16_functional_Firth.pdf"
)

DOWNLOADS_DIR = (
    Path.home()
    / "Downloads"
    / "CLOCK16_FINAL"
    / "Figures"
)

DOWNLOADS_PNG = (
    DOWNLOADS_DIR
    / "Figure1D_CLOCK16_functional_Firth.png"
)

DOWNLOADS_PDF = (
    DOWNLOADS_DIR
    / "Figure1D_CLOCK16_functional_Firth.pdf"
)

# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INPUT, sep="\t")

# ============================================================
# QC
# ============================================================

assert len(df) == 3, f"Expected 3 endpoints, found {len(df)}"

required = [
    "Endpoint",
    "N",
    "Cases",
    "Controls",
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

# Required endpoint names
expected = {
    "Nonsynonymous",
    "Damaging",
    "Strong_Damaging"
}

assert set(df["Endpoint"]) == expected

# ============================================================
# ORDER
# ============================================================

order = [
    "Nonsynonymous",
    "Damaging",
    "Strong_Damaging"
]

df["Endpoint"] = pd.Categorical(
    df["Endpoint"],
    categories=order,
    ordered=True
)

df = df.sort_values(
    "Endpoint"
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
    figsize=(8.5, 4.8)
)

ax.errorbar(
    OR,
    y,
    xerr=xerr,
    fmt="o",
    markersize=6,
    capsize=3,
    linewidth=1.3
)

# Highlight significant endpoints
ax.scatter(
    OR[significant],
    y[significant],
    s=40
)

# Null
ax.axvline(
    1,
    linestyle="--",
    linewidth=1
)

# Log scale
ax.set_xscale("log")

# Labels
ax.set_yticks(y)
ax.set_yticklabels(
    [
        "Nonsynonymous",
        "Damaging",
        "Strong damaging"
    ],
    fontsize=10
)

ax.invert_yaxis()

ax.set_xlabel(
    "Odds ratio per log1p non-CLOCK16 mutation burden (95% CI)",
    fontsize=10.5
)

ax.set_ylabel(
    "Functional mutation class",
    fontsize=11
)

ax.set_title(
    "Functional CLOCK16 mutation associations",
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
    1.04,
    "D",
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
print("FIGURE 1D — FINAL")
print("=" * 80)

print(f"Endpoints:          {len(df)}")
print(f"FDR < 0.05:         {significant.sum()}")
print(f"FDR >= 0.05:        {(~significant).sum()}")

print("\nFinal results:")
print(
    df[
        [
            "Endpoint",
            "Cases",
            "Controls",
            "OR_per_log1p_non_CLOCK16_mutation",
            "CI95_Low",
            "CI95_High",
            "Firth_LR_P",
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
print("FIGURE 1D CREATED SUCCESSFULLY")
print("=" * 80)
