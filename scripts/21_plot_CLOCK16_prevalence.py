import pandas as pd
import matplotlib.pyplot as plt
import os

INPUT = "results/clock16/CLOCK16_final_cancer_prevalence.tsv"
OUT = "results/clock16/Figure1A_CLOCK16_prevalence.pdf"
PNG = "results/clock16/Figure1A_CLOCK16_prevalence.png"

df = pd.read_csv(INPUT, sep="\t")

# ------------------------------------------------------------
# Prepare
# ------------------------------------------------------------

df = df.sort_values("Prevalence_Percent", ascending=True).reset_index(drop=True)

# Shorter labels for figure readability
label_map = {
    "skin cutaneous melanoma": "Skin cutaneous melanoma",
    "bladder urothelial carcinoma": "Bladder urothelial carcinoma",
    "uterine corpus endometrioid carcinoma": "Uterine corpus endometrioid carcinoma",
    "stomach adenocarcinoma": "Stomach adenocarcinoma",
    "lung squamous cell carcinoma": "Lung squamous cell carcinoma",
    "cervical & endocervical cancer": "Cervical & endocervical cancer",
    "lung adenocarcinoma": "Lung adenocarcinoma",
    "esophageal carcinoma": "Esophageal carcinoma",
    "colon adenocarcinoma": "Colon adenocarcinoma",
    "liver hepatocellular carcinoma": "Liver hepatocellular carcinoma",
    "head & neck squamous cell carcinoma": "Head & neck squamous cell carcinoma",
    "kidney papillary cell carcinoma": "Kidney papillary cell carcinoma",
    "uterine carcinosarcoma": "Uterine carcinosarcoma",
    "adrenocortical cancer": "Adrenocortical cancer",
    "sarcoma": "Sarcoma",
    "diffuse large B-cell lymphoma": "Diffuse large B-cell lymphoma",
    "rectum adenocarcinoma": "Rectum adenocarcinoma",
    "breast invasive carcinoma": "Breast invasive carcinoma",
    "kidney clear cell carcinoma": "Kidney clear cell carcinoma",
    "prostate adenocarcinoma": "Prostate adenocarcinoma",
    "kidney chromophobe": "Kidney chromophobe",
    "thymoma": "Thymoma",
    "testicular germ cell tumor": "Testicular germ cell tumor",
    "pancreatic adenocarcinoma": "Pancreatic adenocarcinoma",
    "glioblastoma multiforme": "Glioblastoma multiforme",
    "mesothelioma": "Mesothelioma",
    "pheochromocytoma & paraganglioma": "Pheochromocytoma & paraganglioma",
    "cholangiocarcinoma": "Cholangiocarcinoma",
    "brain lower grade glioma": "Brain lower grade glioma",
    "uveal melanoma": "Uveal melanoma",
    "thyroid carcinoma": "Thyroid carcinoma",
    "ovarian serous cystadenocarcinoma": "Ovarian serous cystadenocarcinoma"
}

df["Label"] = df["Cancer_Type"].map(label_map).fillna(df["Cancer_Type"])

# ------------------------------------------------------------
# Overall prevalence
# ------------------------------------------------------------

overall = 1165 / 10593 * 100

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig_height = max(9, len(df) * 0.34)

fig, ax = plt.subplots(
    figsize=(8.5, fig_height)
)

y = range(len(df))

x = df["Prevalence_Percent"].values

xerr_low = (
    df["Prevalence_Percent"]
    - df["CI95_Low_Percent"]
).values

xerr_high = (
    df["CI95_High_Percent"]
    - df["Prevalence_Percent"]
).values

ax.errorbar(
    x,
    y,
    xerr=[xerr_low, xerr_high],
    fmt="o",
    markersize=5,
    capsize=3,
    linewidth=1.2,
    elinewidth=1.0
)

# Overall prevalence reference
ax.axvline(
    overall,
    linestyle="--",
    linewidth=1
)

ax.text(
    overall + 0.25,
    len(df) - 0.5,
    f"Overall = {overall:.2f}%",
    fontsize=9,
    va="top"
)

# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------

ax.set_yticks(list(y))
ax.set_yticklabels(df["Label"])

ax.set_xlabel(
    "CLOCK16 mutation prevalence (%)",
    fontsize=11
)

ax.set_ylabel(
    "Cancer type",
    fontsize=11
)

ax.set_title(
    "CLOCK16 Mutation Prevalence Across TCGA Cancer Types",
    fontsize=13,
    pad=12
)

# Keep plot clean
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(
    axis="y",
    labelsize=8
)

ax.tick_params(
    axis="x",
    labelsize=9
)

ax.grid(
    axis="x",
    linestyle=":",
    linewidth=0.6,
    alpha=0.5
)

plt.tight_layout()

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

os.makedirs(
    os.path.dirname(OUT),
    exist_ok=True
)

fig.savefig(
    OUT,
    bbox_inches="tight"
)

fig.savefig(
    PNG,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print("=" * 80)
print("CLOCK16 PREVALENCE FIGURE")
print("=" * 80)

print(f"Input: {INPUT}")
print(f"Overall prevalence: {overall:.4f}%")
print(f"Cancer types plotted: {len(df)}")

print("\nHighest prevalence:")
print(
    df.nlargest(5, "Prevalence_Percent")[
        [
            "Cancer_Type",
            "N",
            "CLOCK16_Mutated",
            "Prevalence_Percent"
        ]
    ].to_string(index=False)
)

print("\nSaved:")
print(OUT)
print(PNG)

print("=" * 80)
