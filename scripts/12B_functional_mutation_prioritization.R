
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project <- path.expand("~/circadian_pan_cancer_Q2")

input_file <- file.path(
  project,
  "results/clock16/REAL_primary_CLOCK16_mutations.tsv.gz"
)

outdir <- file.path(
  project,
  "results/clock16"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("============================================================\n")
cat("12B — CLOCK16 FUNCTIONAL MUTATION PRIORITIZATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

x <- fread(input_file)

cat("Input rows:", nrow(x), "\n")

# ------------------------------------------------------------
# BASIC INPUT QC
# ------------------------------------------------------------

required_cols <- c(
  "sample",
  "gene",
  "effect",
  "Amino_Acid_Change",
  "SIFT",
  "PolyPhen",
  "patient_id",
  "effect_class"
)

missing_cols <- setdiff(required_cols, names(x))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# KEEP NONSYNONYMOUS EVENTS
# ------------------------------------------------------------

ns <- x[
  effect_class == "nonsynonymous"
]

cat("Nonsynonymous events:", nrow(ns), "\n")
cat("Unique patients:", uniqueN(ns$patient_id), "\n")
cat(
  "Unique patient-gene combinations:",
  uniqueN(ns[, .(patient_id, gene)]),
  "\n\n"
)

# ------------------------------------------------------------
# NORMALIZE SIFT
# ------------------------------------------------------------

ns[
  ,
  SIFT_Prediction := fifelse(
    is.na(SIFT) | SIFT == "",
    NA_character_,
    fifelse(
      grepl("^deleterious", SIFT, ignore.case = TRUE),
      "deleterious",
      fifelse(
        grepl("^tolerated", SIFT, ignore.case = TRUE),
        "tolerated",
        NA_character_
      )
    )
  )
]

# ------------------------------------------------------------
# NORMALIZE POLYPHEN
# ------------------------------------------------------------

ns[
  ,
  PolyPhen_Prediction := fifelse(
    is.na(PolyPhen) | PolyPhen == "",
    NA_character_,
    fifelse(
      grepl("^probably_damaging", PolyPhen, ignore.case = TRUE),
      "probably_damaging",
      fifelse(
        grepl("^possibly_damaging", PolyPhen, ignore.case = TRUE),
        "possibly_damaging",
        fifelse(
          grepl("^benign", PolyPhen, ignore.case = TRUE),
          "benign",
          NA_character_
        )
      )
    )
  )
]

# ------------------------------------------------------------
# FUNCTIONAL PRIORITY CLASS
# ------------------------------------------------------------

ns[
  ,
  Functional_Class := fcase(

    SIFT_Prediction == "deleterious" &
      PolyPhen_Prediction == "probably_damaging",
    "strong_damaging",

    SIFT_Prediction == "deleterious" &
      PolyPhen_Prediction == "possibly_damaging",
    "damaging",

    SIFT_Prediction == "deleterious" &
      PolyPhen_Prediction == "benign",
    "discordant_SIFT_deleterious",

    SIFT_Prediction == "tolerated" &
      PolyPhen_Prediction %in%
        c("probably_damaging", "possibly_damaging"),
    "discordant_PolyPhen_damaging",

    SIFT_Prediction == "tolerated" &
      PolyPhen_Prediction == "benign",
    "likely_benign",

    default = "prediction_missing"
  )
]

# ------------------------------------------------------------
# FUNCTIONAL CLASS QC
# ------------------------------------------------------------

cat("\n=== FUNCTIONAL CLASSIFICATION ===\n")

print(
  ns[
    ,
    .N,
    by = Functional_Class
  ][order(-N)]
)

# ------------------------------------------------------------
# STRONG DAMAGING EVENTS
# ------------------------------------------------------------

strong_damaging <- ns[
  Functional_Class == "strong_damaging"
]

damaging <- ns[
  Functional_Class %in%
    c("strong_damaging", "damaging")
]

cat("\n=== STRONG DAMAGING ===\n")
cat("Strong damaging events:", nrow(strong_damaging), "\n")
cat(
  "Strong damaging patients:",
  uniqueN(strong_damaging$patient_id),
  "\n"
)
cat(
  "Strong damaging patient-gene combinations:",
  uniqueN(
    strong_damaging[, .(patient_id, gene)]
  ),
  "\n"
)

cat("\n=== DAMAGING (SIFT DELETERIOUS) ===\n")
cat("Damaging events:", nrow(damaging), "\n")
cat(
  "Damaging patients:",
  uniqueN(damaging$patient_id),
  "\n"
)

# ------------------------------------------------------------
# FUNCTIONAL SUMMARY
# ------------------------------------------------------------

functional_summary <- ns[
  ,
  .(
    Nonsynonymous_Events = .N,
    Strong_Damaging = sum(
      Functional_Class == "strong_damaging"
    ),
    Damaging = sum(
      Functional_Class %in%
        c("strong_damaging", "damaging")
    ),
    Likely_Benign = sum(
      Functional_Class == "likely_benign"
    ),
    Prediction_Missing = sum(
      Functional_Class == "prediction_missing"
    )
  )
]

# ------------------------------------------------------------
# PRIORITY SUMMARY
# ------------------------------------------------------------

priority_summary <- ns[
  ,
  .(
    Events = .N,
    Patients = uniqueN(patient_id),
    Patient_Gene_Combinations =
      uniqueN(.SD[, .(patient_id, gene)])
  ),
  by = Functional_Class
][order(-Events)]

# ------------------------------------------------------------
# GENE-LEVEL FUNCTIONAL LANDSCAPE
# ------------------------------------------------------------

gene_functional <- ns[
  ,
  .(
    Nonsynonymous_Events = .N,

    Patients =
      uniqueN(patient_id),

    Strong_Damaging_Events =
      sum(
        Functional_Class ==
          "strong_damaging"
      ),

    Strong_Damaging_Patients =
      uniqueN(
        patient_id[
          Functional_Class ==
            "strong_damaging"
        ]
      ),

    Damaging_Events =
      sum(
        Functional_Class %in%
          c(
            "strong_damaging",
            "damaging"
          )
      ),

    Damaging_Patients =
      uniqueN(
        patient_id[
          Functional_Class %in%
            c(
              "strong_damaging",
              "damaging"
            )
        ]
      ),

    SIFT_Deleterious =
      sum(
        SIFT_Prediction ==
          "deleterious",
        na.rm = TRUE
      ),

    PolyPhen_Damaging =
      sum(
        PolyPhen_Prediction %in%
          c(
            "probably_damaging",
            "possibly_damaging"
          ),
        na.rm = TRUE
      )
  ),
  by = gene
]

setorder(
  gene_functional,
  -Strong_Damaging_Events,
  -Damaging_Events
)

# ------------------------------------------------------------
# PATIENT FUNCTIONAL BURDEN
# ------------------------------------------------------------

patient_functional <- ns[
  ,
  .(
    Nonsynonymous_Events = .N,

    Nonsynonymous_Genes =
      uniqueN(gene),

    Strong_Damaging_Events =
      sum(
        Functional_Class ==
          "strong_damaging"
      ),

    Damaging_Events =
      sum(
        Functional_Class %in%
          c(
            "strong_damaging",
            "damaging"
          )
      ),

    Strong_Damaging_Genes =
      uniqueN(
        gene[
          Functional_Class ==
            "strong_damaging"
        ]
      ),

    Damaging_Genes =
      uniqueN(
        gene[
          Functional_Class %in%
            c(
              "strong_damaging",
              "damaging"
            )
        ]
      )
  ),
  by = patient_id
]

# ------------------------------------------------------------
# SAVE FULL CLASSIFICATION
# ------------------------------------------------------------

fwrite(
  ns,
  file.path(
    outdir,
    "REAL_CLOCK16_nonsynonymous_functional_classification.tsv.gz"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# SAVE SUMMARY FILES
# ------------------------------------------------------------

fwrite(
  functional_summary,
  file.path(
    outdir,
    "REAL_CLOCK16_functional_class_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  priority_summary,
  file.path(
    outdir,
    "REAL_CLOCK16_functional_priority_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_functional,
  file.path(
    outdir,
    "REAL_CLOCK16_gene_functional_landscape.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# SAVE STRONG DAMAGING EVENTS
# ------------------------------------------------------------

fwrite(
  strong_damaging,
  file.path(
    outdir,
    "REAL_CLOCK16_strong_damaging_events.tsv.gz"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# SAVE GENE-LEVEL DAMAGING
# ------------------------------------------------------------

damaging_by_gene <- ns[
  Functional_Class %in%
    c(
      "strong_damaging",
      "damaging"
    ),
  .(
    Damaging_Events = .N,
    Damaging_Patients =
      uniqueN(patient_id)
  ),
  by = gene
][order(-Damaging_Events)]

fwrite(
  damaging_by_gene,
  file.path(
    outdir,
    "REAL_CLOCK16_damaging_by_gene.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# SAVE PATIENT FUNCTIONAL BURDEN
# ------------------------------------------------------------

fwrite(
  patient_functional,
  file.path(
    outdir,
    "REAL_CLOCK16_patient_functional_burden.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# FINAL QC
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("12B FINAL QC\n")
cat("============================================================\n")

cat(
  "Nonsynonymous events:",
  nrow(ns),
  "\n"
)

cat(
  "Strong damaging events:",
  nrow(strong_damaging),
  "\n"
)

cat(
  "Strong damaging patients:",
  uniqueN(strong_damaging$patient_id),
  "\n"
)

cat(
  "Damaging events:",
  nrow(damaging),
  "\n"
)

cat(
  "Damaging patients:",
  uniqueN(damaging$patient_id),
  "\n"
)

cat("\n=== FILES SAVED ===\n")

cat(
  "REAL_CLOCK16_nonsynonymous_functional_classification.tsv.gz\n",
  "REAL_CLOCK16_functional_class_summary.tsv\n",
  "REAL_CLOCK16_functional_priority_summary.tsv\n",
  "REAL_CLOCK16_gene_functional_landscape.tsv\n",
  "REAL_CLOCK16_strong_damaging_events.tsv.gz\n",
  "REAL_CLOCK16_damaging_by_gene.tsv\n",
  "REAL_CLOCK16_patient_functional_burden.tsv\n"
)

cat("\n12B COMPLETED SUCCESSFULLY.\n")

