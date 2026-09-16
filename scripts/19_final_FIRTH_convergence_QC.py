import os
import pandas as pd
import numpy as np

BASE = "results/clock16"

FILES = {
    "pooled": f"{BASE}/PAN_CANCER_pooled_FIRTH_LR.tsv",
    "burden_cancer": f"{BASE}/PAN_CANCER_burden_adjusted_FIRTH_LR.tsv",
    "gene": f"{BASE}/CLOCK16_gene_level_pooled_FIRTH_LR.tsv",
    "functional": f"{BASE}/CLOCK16_functional_CORRECTED_pooled_FIRTH_LR.tsv",
}

print("=" * 80)
print("FINAL CLOCK16 FIRTH ANALYSIS CONVERGENCE / RESULT QC")
print("=" * 80)

for name, path in FILES.items():

    print("\n" + "-" * 80)
    print(name.upper())
    print("-" * 80)
    print("File:", path)

    if not os.path.exists(path):
        print("STATUS: MISSING")
        continue

    df = pd.read_csv(path, sep="\t")

    print("Rows:", len(df))
    print("Columns:", list(df.columns))

    # Detect P-value column
    pcols = [c for c in df.columns if c in
             ["Firth_LR_P", "Firth_LR_p", "P", "p"]]

    if pcols:
        pcol = pcols[0]
        p = pd.to_numeric(df[pcol], errors="coerce")

        print("P-value column:", pcol)
        print("Valid P-values:", p.notna().sum())
        print("P = 0:", (p == 0).sum())
        print("P finite:", np.isfinite(p).sum())

    # FDR
    fdrcols = [c for c in df.columns if c.upper() == "FDR"]

    if fdrcols:
        fdr = pd.to_numeric(df[fdrcols[0]], errors="coerce")
        print("FDR < 0.05:", (fdr < 0.05).sum())

    # OR / CI
    if "OR_per_log1p_non_CLOCK16_mutation" in df.columns:

        OR = pd.to_numeric(
            df["OR_per_log1p_non_CLOCK16_mutation"],
            errors="coerce"
        )

        print("OR finite:", np.isfinite(OR).sum())
        print("OR <= 0:", (OR <= 0).sum())

    if "CI95_Low" in df.columns and "CI95_High" in df.columns:

        low = pd.to_numeric(df["CI95_Low"], errors="coerce")
        high = pd.to_numeric(df["CI95_High"], errors="coerce")

        print("Invalid CI:", ((low <= 0) | (high <= 0) | (low > high)).sum())

    # Direction
    if "Direction" in df.columns:
        print("Directions:")
        print(df["Direction"].value_counts(dropna=False).to_string())

print("\n" + "=" * 80)
print("FINAL QC COMPLETE")
print("=" * 80)
