import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

INPUT = Path(
    "results/clock16/REAL_CLOCK16_cancer_association_any_v2.tsv"
)

PROJECT_PNG = Path(
    "results/clock16/Figure1B_CLOCK16_cancer_association.png"
)

PROJECT_PDF = Path(
    "results/clock16/Figure1B_CLOCK16_cancer_association.pdf"
)

DOWNLOADS_DIR = Path.home() / "Downloads" / "CLOCK16_FINAL" / "Figures"

DOWNLOADS_PNG = DOWNLOADS_DIR / \
    "Figure1B_CLOCK16_cancer_association.png"

DOWNLOADS_PDF = DOWNLOADS_DIR / \
    "Figure1B_CLOCK16_cancer_association.pdf"


# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INPUT, sep="\t")


# ============================================================
# QC
# ============================================================

assert len(df) == 32, f"Expected 32 cancer types, found {len(df)}"
assert df["cancer_type"].nunique() == 32
assert df["phenotype"].nunique() == 1
assert df["phenotype"].iloc[0] == "Any_CLOCK16"

assert (df["odds_ratio"] > 0).all()
assert (df["lower_95CI"] > 0).all()
assert (df["upper_95CI"] > 0).all()
assert (df["lower_95CI"] <= df["upper_95CI"]).all()

# Sort strongest enrichment -> strongest depletion
df = df.sort_values(
    "odds_ratio",
    ascending=False
).reset_index(drop=True)


# ============================================================
# LABELS
# ============================================================

label_map = {
    "skin cutaneous melanoma": "Skin melanoma",
    "uterine corpus endometrioid carcinoma": "Uterine endometrioid",
    "bladder urothelial carcinoma": "Bladder urothelial",
    "stomach adenocarcinoma": "Stomach adenocarcinoma",
    "lung squamous cell carcinoma": "Lung squamous",
    "cervical & endocervical cancer": "Cervical/endocervical",
    "lung adenocarcinoma": "Lung adenocarcinoma",
    "esophageal carcinoma": "Esophageal carcinoma",
    "colon adenocarcinoma": "Colon adenocarcinoma",
    "liver hepatocellular carcinoma": "Liver hepatocellular",
    "head & neck squamous cell carcinoma": "Head & neck squamous",
    "kidney papillary cell carcinoma": "Kidney papillary",
    "uterine carcinosarcoma": "Uterine carcinosarcoma",
    "adrenocortical cancer": "Adrenocortical",
    "sarcoma": "Sarcoma",
    "diffuse large B-cell lymphoma": "DLBCL",
    "rectum adenocarcinoma": "Rectal adenocarcinoma",
    "breast invasive carcinoma": "Breast carcinoma",
    "kidney clear cell carcinoma": "Kidney clear cell",
    "prostate adenocarcinoma": "Prostate adenocarcinoma",
    "kidney chromophobe": "Kidney chromophobe",
    "thymoma": "Thymoma",
    "testicular germ cell tumor": "Testicular germ cell",
    "pancreatic adenocarcinoma": "Pancreatic adenocarcinoma",
    "glioblastoma multiforme": "Glioblastoma",
    "mesothelioma": "Mesothelioma",
    "cholangiocarcinoma": "Cholangiocarcinoma",
    "pheochromocytoma & paraganglioma": "Pheochromocytoma/paraganglioma",
    "uveal melanoma": "Uveal melanoma",
    "brain lower grade glioma": "Lower-grade glioma",
    "thyroid carcinoma": "Thyroid carcinoma",
    "ovarian serous cystadenocarcinoma": "Ovarian serous"
}

df["Label"] = df["cancer_type"].map(label_map)


# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(figsize=(9, 11))

y = np.arange(len(df))

OR = df["odds_ratio"].to_numpy()
LOW = df["lower_95CI"].to_numpy()
HIGH = df["upper_95CI"].to_numpy()

xerr = np.vstack([
    OR - LOW,
    HIGH - OR
])

significant = df["FDR"] < 0.05


# All estimates
ax.errorbar(
    OR,
    y,
    xerr=xerr,
    fmt="o",
    markersize=4.5,
    capsize=2.5,
    linewidth=1.1
)

# Significant estimates
ax.scatter(
    OR[significant],
    y[significant],
    s=28
)

# Null
ax.axvline(
    1,
    linestyle="--",
    linewidth=1
)


# ============================================================
# AXES
# ============================================================

ax.set_xscale("log")

ax.set_yticks(y)
ax.set_yticklabels(df["Label"], fontsize=9)

ax.invert_yaxis()

ax.set_xlabel(
    "Odds ratio (95% CI)",
    fontsize=11
)

ax.set_title(
    "Cancer-specific CLOCK16 mutation association",
    fontsize=13,
    fontweight="bold",
    pad=12
)


# ============================================================
# GRID / CLEAN STYLE
# ============================================================

ax.grid(
    axis="x",
    linestyle=":",
    linewidth=0.6,
    alpha=0.5
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


# ============================================================
# PANEL LABEL
# ============================================================

ax.text(
    -0.12,
    1.02,
    "B",
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
# FINAL QC
# ============================================================

print("=" * 80)
print("FIGURE 1B — FINAL")
print("=" * 80)

print(f"Cancer types:       {len(df)}")
print(f"FDR < 0.05:         {significant.sum()}")
print(f"Enriched OR > 1:    {(OR > 1).sum()}")
print(f"Depleted OR < 1:    {(OR < 1).sum()}")

print("\nOR range:")
print(
    f"{OR.min():.6f} - {OR.max():.6f}"
)

print("\nProject files:")
print(PROJECT_PNG)
print(PROJECT_PDF)

print("\nDownloads files:")
print(DOWNLOADS_PNG)
print(DOWNLOADS_PDF)

print("\n" + "=" * 80)
print("FIGURE 1B CREATED SUCCESSFULLY")
print("=" * 80)
