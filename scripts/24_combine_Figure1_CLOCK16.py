import matplotlib.pyplot as plt
from matplotlib.image import imread
from pathlib import Path

# ============================================================
# INPUTS
# ============================================================

fig_dir = (
    Path.home()
    / "Downloads"
    / "CLOCK16_FINAL"
    / "Figures"
)

files = [
    fig_dir / "Figure1A_CLOCK16_prevalence.png",
    fig_dir / "Figure1B_CLOCK16_cancer_association.png",
    fig_dir / "Figure1C_CLOCK16_gene_Firth.png",
    fig_dir / "Figure1D_CLOCK16_functional_Firth.png",
]

for f in files:
    assert f.exists(), f"Missing figure: {f}"

# ============================================================
# LOAD IMAGES
# ============================================================

images = [imread(str(f)) for f in files]

# ============================================================
# COMBINED FIGURE
# ============================================================

fig, axes = plt.subplots(
    2,
    2,
    figsize=(16, 13)
)

labels = ["A", "B", "C", "D"]

for ax, img, label in zip(
    axes.flatten(),
    images,
    labels
):
    ax.imshow(img)
    ax.axis("off")

    ax.text(
        0.01,
        0.99,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="left"
    )

# ============================================================
# LAYOUT
# ============================================================

plt.subplots_adjust(
    left=0.02,
    right=0.98,
    bottom=0.02,
    top=0.98,
    wspace=0.03,
    hspace=0.04
)

# ============================================================
# OUTPUT
# ============================================================

png_out = fig_dir / "Figure1_CLOCK16_pan_cancer_ABCD.png"
pdf_out = fig_dir / "Figure1_CLOCK16_pan_cancer_ABCD.pdf"

fig.savefig(
    png_out,
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.05
)

fig.savefig(
    pdf_out,
    bbox_inches="tight",
    pad_inches=0.05
)

plt.close(fig)

# ============================================================
# QC
# ============================================================

print("=" * 80)
print("COMBINED FIGURE 1 — FINAL")
print("=" * 80)

for f in files:
    print("Input:", f)

print("\nCombined PNG:")
print(png_out)

print("\nCombined PDF:")
print(pdf_out)

print("\n" + "=" * 80)
print("FIGURE 1 A–D COMBINATION COMPLETED")
print("=" * 80)
