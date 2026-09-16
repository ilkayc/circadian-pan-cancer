import fitz
from pathlib import Path

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

project = Path("results/clock16")
outdir = Path.home() / "Downloads" / "CLOCK16_FINAL" / "Figures"
outdir.mkdir(parents=True, exist_ok=True)

panels = [
    ("A", project / "Figure1A_CLOCK16_prevalence.pdf"),
    ("B", project / "Figure1B_CLOCK16_cancer_association.pdf"),
    ("C", project / "Figure1C_CLOCK16_gene_Firth.pdf"),
    ("D", project / "Figure1D_CLOCK16_functional_Firth.pdf"),
]

for label, path in panels:
    if not path.exists():
        raise FileNotFoundError(f"Missing panel {label}: {path}")

# ------------------------------------------------------------
# READ PDF PAGES
# ------------------------------------------------------------

docs = []
for label, path in panels:
    doc = fitz.open(path)
    if len(doc) != 1:
        raise ValueError(f"{path} has {len(doc)} pages; expected 1.")
    docs.append((label, doc))

# ------------------------------------------------------------
# TARGET PAGE
# Approx. journal-style 2-column figure
# ------------------------------------------------------------

CELL_W = 360
CELL_H = 300

PAGE_W = CELL_W * 2
PAGE_H = CELL_H * 2

margin = 10
gap = 8

# ------------------------------------------------------------
# CREATE VECTOR PDF
# ------------------------------------------------------------

out_pdf = outdir / "Figure1_CLOCK16_pan_cancer_combined.pdf"

combined = fitz.open()

page = combined.new_page(
    width=PAGE_W,
    height=PAGE_H
)

positions = {
    "A": (margin, margin),
    "B": (CELL_W + gap, margin),
    "C": (margin, CELL_H + gap),
    "D": (CELL_W + gap, CELL_H + gap),
}

for label, doc in docs:
    src = doc[0]
    x, y = positions[label]

    target = fitz.Rect(
        x,
        y,
        x + CELL_W - margin,
        y + CELL_H - margin
    )

    # Preserve vector content
    page.show_pdf_page(
        target,
        doc,
        0,
        keep_proportion=True
    )

    # Panel label
    page.insert_text(
        (x + 2, y + 14),
        label,
        fontsize=12,
        fontname="helv",
        color=(0, 0, 0),
    )

# ------------------------------------------------------------
# SAVE VECTOR PDF
# ------------------------------------------------------------

combined.save(
    out_pdf,
    garbage=4,
    deflate=True,
    clean=True
)

combined.close()

# ------------------------------------------------------------
# HIGH-RES PNG
# 600 DPI
# ------------------------------------------------------------

png_path = outdir / "Figure1_CLOCK16_pan_cancer_combined_600dpi.png"

doc = fitz.open(out_pdf)
page = doc[0]

pix = page.get_pixmap(
    matrix=fitz.Matrix(600 / 72, 600 / 72),
    alpha=False,
    annots=False
)

pix.save(png_path)

doc.close()

print("=" * 80)
print("FIGURE 1 — FINAL COMBINED")
print("=" * 80)
print()
print("Panels:")
for label, path in panels:
    print(f"  {label}: {path}")

print()
print("Vector PDF:")
print(out_pdf)

print()
print("600 DPI PNG:")
print(png_path)

print()
print("PDF size:", out_pdf.stat().st_size / 1024 / 1024, "MB")
print("PNG size:", png_path.stat().st_size / 1024 / 1024, "MB")

print()
print("DONE")
