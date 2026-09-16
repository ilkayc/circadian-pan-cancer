
library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

clock_file <- file.path(
  project,
  "results/clock16/REAL_primary_CLOCK16_mutations.tsv.gz"
)

outdir <- file.path(
  project,
  "results/clock16"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clock <- fread(clock_file)

cat("\n============================================================\n")
cat(" 12A — CLOCK16 FUNCTIONAL PREDICTION QC\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Restrict to nonsynonymous events
# ------------------------------------------------------------

nonsyn <- clock[
  effect_class == "nonsynonymous"
]

cat("Nonsynonymous events:", nrow(nonsyn), "\n")
cat(
  "Unique patients:",
  uniqueN(nonsyn$patient_id),
  "\n"
)
cat(
  "Unique patient-gene combinations:",
  uniqueN(
    nonsyn[, .(patient_id, gene)]
  ),
  "\n\n"
)

# ------------------------------------------------------------
# 2. Amino-acid annotation completeness
# ------------------------------------------------------------

aa_qc <- data.table(
  Annotation = c(
    "With_Amino_Acid_Change",
    "Without_Amino_Acid_Change"
  ),
  N = c(
    sum(
      !is.na(nonsyn$Amino_Acid_Change) &
      nonsyn$Amino_Acid_Change != ""
    ),
    sum(
      is.na(nonsyn$Amino_Acid_Change) |
      nonsyn$Amino_Acid_Change == ""
    )
  )
)

aa_qc[
  ,
  Percent := 100 * N / nrow(nonsyn)
]

cat("=== AMINO ACID ANNOTATION ===\n")
print(aa_qc)
cat("\n")

# ------------------------------------------------------------
# 3. SIFT availability
# ------------------------------------------------------------

sift_available <- !is.na(nonsyn$SIFT) &
  nonsyn$SIFT != ""

sift_qc <- data.table(
  Annotation = c(
    "SIFT_Available",
    "SIFT_Missing"
  ),
  N = c(
    sum(sift_available),
    sum(!sift_available)
  )
)

sift_qc[
  ,
  Percent := 100 * N / nrow(nonsyn)
]

cat("=== SIFT AVAILABILITY ===\n")
print(sift_qc)
cat("\n")

# ------------------------------------------------------------
# 4. PolyPhen availability
# ------------------------------------------------------------

poly_available <- !is.na(nonsyn$PolyPhen) &
  nonsyn$PolyPhen != ""

poly_qc <- data.table(
  Annotation = c(
    "PolyPhen_Available",
    "PolyPhen_Missing"
  ),
  N = c(
    sum(poly_available),
    sum(!poly_available)
  )
)

poly_qc[
  ,
  Percent := 100 * N / nrow(nonsyn)
]

cat("=== POLYPHEN AVAILABILITY ===\n")
print(poly_qc)
cat("\n")

# ------------------------------------------------------------
# 5. SIFT raw categories
# ------------------------------------------------------------

sift_raw <- nonsyn[
  sift_available,
  .(
    SIFT_Raw = SIFT,
    N = .N
  ),
  by = SIFT
]

setnames(
  sift_raw,
  "SIFT",
  "SIFT_Raw"
)

setorder(sift_raw, -N)

cat("=== SIFT RAW ANNOTATIONS ===\n")
print(sift_raw)
cat("\n")

# ------------------------------------------------------------
# 6. PolyPhen raw categories
# ------------------------------------------------------------

poly_raw <- nonsyn[
  poly_available,
  .(
    PolyPhen_Raw = PolyPhen,
    N = .N
  ),
  by = PolyPhen
]

setnames(
  poly_raw,
  "PolyPhen",
  "PolyPhen_Raw"
)

setorder(poly_raw, -N)

cat("=== POLYPHEN RAW ANNOTATIONS ===\n")
print(poly_raw)
cat("\n")

# ------------------------------------------------------------
# 7. Parse SIFT prediction category
# ------------------------------------------------------------

nonsyn[
  ,
  SIFT_Prediction := fifelse(
    grepl(
      "^deleterious",
      SIFT,
      ignore.case = TRUE
    ),
    "deleterious",
    fifelse(
      grepl(
        "^tolerated",
        SIFT,
        ignore.case = TRUE
      ),
      "tolerated",
      NA_character_
    )
  )
]

cat("=== SIFT PREDICTION CLASS ===\n")

print(
  nonsyn[
    ,
    .(
      N = .N
    ),
    by = SIFT_Prediction
  ][order(-N)]
)

cat("\n")

# ------------------------------------------------------------
# 8. Parse PolyPhen prediction category
# ------------------------------------------------------------

nonsyn[
  ,
  PolyPhen_Prediction := fifelse(
    grepl(
      "^probably_damaging",
      PolyPhen,
      ignore.case = TRUE
    ),
    "probably_damaging",
    fifelse(
      grepl(
        "^possibly_damaging",
        PolyPhen,
        ignore.case = TRUE
      ),
      "possibly_damaging",
      fifelse(
        grepl(
          "^benign",
          PolyPhen,
          ignore.case = TRUE
        ),
        "benign",
        NA_character_
      )
    )
  )
]

cat("=== POLYPHEN PREDICTION CLASS ===\n")

print(
  nonsyn[
    ,
    .(
      N = .N
    ),
    by = PolyPhen_Prediction
  ][order(-N)]
)

cat("\n")

# ------------------------------------------------------------
# 9. SIFT × PolyPhen cross-classification
# ------------------------------------------------------------

cat("=== SIFT × POLYPHEN ===\n")

prediction_cross <- nonsyn[
  ,
  .(N = .N),
  by = .(
    SIFT_Prediction,
    PolyPhen_Prediction
  )
][order(-N)]

print(prediction_cross)

cat("\n")

# ------------------------------------------------------------
# 10. Gene-level prediction landscape
# ------------------------------------------------------------

cat("=== GENE × FUNCTIONAL PREDICTION ===\n")

gene_prediction <- nonsyn[
  ,
  .(
    Nonsynonymous_Events = .N,
    SIFT_Available = sum(
      !is.na(SIFT_Prediction)
    ),
    SIFT_Deleterious = sum(
      SIFT_Prediction == "deleterious",
      na.rm = TRUE
    ),
    PolyPhen_Available = sum(
      !is.na(PolyPhen_Prediction)
    ),
    PolyPhen_Probably_Damaging = sum(
      PolyPhen_Prediction == "probably_damaging",
      na.rm = TRUE
    ),
    PolyPhen_Possibly_Damaging = sum(
      PolyPhen_Prediction == "possibly_damaging",
      na.rm = TRUE
    ),
    PolyPhen_Benign = sum(
      PolyPhen_Prediction == "benign",
      na.rm = TRUE
    )
  ),
  by = gene
]

setorder(
  gene_prediction,
  -Nonsynonymous_Events
)

print(gene_prediction)

cat("\n")

# ------------------------------------------------------------
# 11. Save outputs
# ------------------------------------------------------------

fwrite(
  aa_qc,
  file.path(
    outdir,
    "REAL_CLOCK16_amino_acid_QC.tsv"
  ),
  sep = "\t"
)

fwrite(
  sift_qc,
  file.path(
    outdir,
    "REAL_CLOCK16_SIFT_availability_QC.tsv"
  ),
  sep = "\t"
)

fwrite(
  poly_qc,
  file.path(
    outdir,
    "REAL_CLOCK16_PolyPhen_availability_QC.tsv"
  ),
  sep = "\t"
)

fwrite(
  sift_raw,
  file.path(
    outdir,
    "REAL_CLOCK16_SIFT_raw_annotations.tsv"
  ),
  sep = "\t"
)

fwrite(
  poly_raw,
  file.path(
    outdir,
    "REAL_CLOCK16_PolyPhen_raw_annotations.tsv"
  ),
  sep = "\t"
)

fwrite(
  prediction_cross,
  file.path(
    outdir,
    "REAL_CLOCK16_SIFT_PolyPhen_cross.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_prediction,
  file.path(
    outdir,
    "REAL_CLOCK16_gene_functional_prediction.tsv"
  ),
  sep = "\t"
)

fwrite(
  nonsyn,
  file.path(
    outdir,
    "REAL_CLOCK16_nonsynonymous_annotated.tsv.gz"
  ),
  sep = "\t"
)

cat("============================================================\n")
cat("12A FILES SAVED\n")
cat("============================================================\n")

