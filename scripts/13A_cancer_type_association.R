
library(data.table)

cat("\n")
cat("============================================================\n")
cat("13A — CLOCK16 CANCER-TYPE ASSOCIATION ANALYSIS\n")
cat("============================================================\n\n")

base_dir <- "results/clock16"
denom_dir <- "results/denominator"

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

denom_file <- file.path(
  denom_dir,
  "primary_tumors_clock16_mutated.tsv"
)

any_mut_file <- file.path(
  base_dir,
  "REAL_primary_CLOCK16_mutations.tsv.gz"
)

nonsyn_file <- file.path(
  base_dir,
  "REAL_CLOCK16_nonsynonymous_functional_classification.tsv.gz"
)

if (!file.exists(denom_file))
  stop("Missing denominator file: ", denom_file)

if (!file.exists(any_mut_file))
  stop("Missing CLOCK16 mutation file: ", any_mut_file)

if (!file.exists(nonsyn_file))
  stop("Missing nonsynonymous functional file: ", nonsyn_file)


# ------------------------------------------------------------
# READ DATA
# ------------------------------------------------------------

denom <- fread(denom_file)
any_mut <- fread(any_mut_file)
nonsyn <- fread(nonsyn_file)

cat("Denominator rows:", nrow(denom), "\n")
cat("Any mutation rows:", nrow(any_mut), "\n")
cat("Nonsynonymous rows:", nrow(nonsyn), "\n\n")


# ------------------------------------------------------------
# HELPER: FIND COLUMN
# ------------------------------------------------------------

find_col <- function(dt, candidates, label) {

  hit <- candidates[candidates %in% names(dt)]

  if (length(hit) == 0) {
    stop(
      "Could not identify ", label,
      " column.\nAvailable columns:\n",
      paste(names(dt), collapse = ", ")
    )
  }

  hit[1]
}


# ------------------------------------------------------------
# IDENTIFY COLUMNS
# ------------------------------------------------------------

denom_patient_col <- find_col(
  denom,
  c("patient_id", "Patient_ID", "patient", "sample"),
  "patient/sample"
)

denom_cancer_col <- find_col(
  denom,
  c("cancer_type", "Cancer_Type", "cancer"),
  "cancer type"
)

any_patient_col <- find_col(
  any_mut,
  c("patient_id", "Patient_ID", "patient"),
  "patient ID in mutation file"
)

nonsyn_patient_col <- find_col(
  nonsyn,
  c("patient_id", "Patient_ID", "patient"),
  "patient ID in nonsynonymous file"
)

cat("Denominator patient column:", denom_patient_col, "\n")
cat("Denominator cancer column:", denom_cancer_col, "\n")
cat("Any mutation patient column:", any_patient_col, "\n")
cat("Nonsynonymous patient column:", nonsyn_patient_col, "\n\n")


# ------------------------------------------------------------
# PATIENT-LEVEL DENOMINATOR
# ------------------------------------------------------------

patient_cancer <- unique(
  denom[, .(
    patient_id = get(denom_patient_col),
    cancer_type = get(denom_cancer_col)
  )]
)

patient_cancer <- patient_cancer[
  !is.na(patient_id) &
  !is.na(cancer_type)
]

cat("Unique primary patients:", uniqueN(patient_cancer$patient_id), "\n")
cat("Cancer types:", uniqueN(patient_cancer$cancer_type), "\n\n")


# ------------------------------------------------------------
# ANY CLOCK16 MUTATION
# ------------------------------------------------------------

any_patients <- unique(
  any_mut[
    !is.na(get(any_patient_col)),
    get(any_patient_col)
  ]
)

any_patients <- any_patients[
  any_patients %in% patient_cancer$patient_id
]

cat("Any CLOCK16 mutated patients in denominator:",
    length(any_patients), "\n")


# ------------------------------------------------------------
# NONSYNONYMOUS CLOCK16 MUTATION
# ------------------------------------------------------------

nonsyn_patients <- unique(
  nonsyn[
    !is.na(get(nonsyn_patient_col)),
    get(nonsyn_patient_col)
  ]
)

nonsyn_patients <- nonsyn_patients[
  nonsyn_patients %in% patient_cancer$patient_id
]

cat("Nonsynonymous CLOCK16 mutated patients in denominator:",
    length(nonsyn_patients), "\n")


# ------------------------------------------------------------
# DAMAGING CLOCK16 MUTATION
# ------------------------------------------------------------

damaging_classes <- c(
  "strong_damaging",
  "damaging"
)

damaging_patients <- unique(
  nonsyn[
    Functional_Class %in% damaging_classes &
    !is.na(get(nonsyn_patient_col)),
    get(nonsyn_patient_col)
  ]
)

damaging_patients <- damaging_patients[
  damaging_patients %in% patient_cancer$patient_id
]

cat("Damaging CLOCK16 mutated patients in denominator:",
    length(damaging_patients), "\n\n")


# ------------------------------------------------------------
# ASSOCIATION FUNCTION
# ------------------------------------------------------------

run_association <- function(
  patient_cancer,
  mutated_patients,
  phenotype_name
) {

  cancer_types <- sort(
    unique(patient_cancer$cancer_type)
  )

  out <- vector(
    "list",
    length(cancer_types)
  )

  for (i in seq_along(cancer_types)) {

    ct <- cancer_types[i]

    ids <- patient_cancer[
      cancer_type == ct,
      patient_id
    ]

    total <- length(ids)

    mutated <- sum(
      ids %in% mutated_patients
    )

    nonmutated <- total - mutated

    outside_total <- nrow(patient_cancer) - total

    outside_mutated <- sum(
      patient_cancer$patient_id %in% mutated_patients &
      !(patient_cancer$cancer_type == ct)
    )

    outside_nonmutated <-
      outside_total - outside_mutated

    tab <- matrix(
      c(
        mutated,
        nonmutated,
        outside_mutated,
        outside_nonmutated
      ),
      nrow = 2,
      byrow = TRUE
    )

    fisher <- fisher.test(tab)

    prevalence <- 100 * mutated / total

    out[[i]] <- data.table(
      phenotype = phenotype_name,
      cancer_type = ct,
      total_primary_patients = total,
      mutated_patients = mutated,
      nonmutated_patients = nonmutated,
      prevalence_percent = prevalence,
      odds_ratio = unname(fisher$estimate),
      lower_95CI = fisher$conf.int[1],
      upper_95CI = fisher$conf.int[2],
      fisher_p = fisher$p.value
    )
  }

  result <- rbindlist(out)

  result[
    ,
    FDR := p.adjust(fisher_p, method = "BH")
  ]

  result[
    order(FDR, fisher_p)
  ]
}


# ------------------------------------------------------------
# RUN THREE ANALYSES
# ------------------------------------------------------------

cat("=== RUNNING ANY CLOCK16 ASSOCIATION ===\n")

res_any <- run_association(
  patient_cancer,
  any_patients,
  "Any_CLOCK16"
)

cat("Completed.\n\n")


cat("=== RUNNING NONSYNONYMOUS ASSOCIATION ===\n")

res_nonsyn <- run_association(
  patient_cancer,
  nonsyn_patients,
  "Nonsynonymous_CLOCK16"
)

cat("Completed.\n\n")


cat("=== RUNNING DAMAGING ASSOCIATION ===\n")

res_damaging <- run_association(
  patient_cancer,
  damaging_patients,
  "Damaging_CLOCK16"
)

cat("Completed.\n\n")


# ------------------------------------------------------------
# COMBINE
# ------------------------------------------------------------

all_results <- rbind(
  res_any,
  res_nonsyn,
  res_damaging
)


# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

fwrite(
  res_any,
  file.path(
    base_dir,
    "REAL_CLOCK16_cancer_association_any.tsv"
  ),
  sep = "\t"
)

fwrite(
  res_nonsyn,
  file.path(
    base_dir,
    "REAL_CLOCK16_cancer_association_nonsynonymous.tsv"
  ),
  sep = "\t"
)

fwrite(
  res_damaging,
  file.path(
    base_dir,
    "REAL_CLOCK16_cancer_association_damaging.tsv"
  ),
  sep = "\t"
)

fwrite(
  all_results,
  file.path(
    base_dir,
    "REAL_CLOCK16_cancer_association_ALL.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# SIGNIFICANT RESULTS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("13A — RESULTS SUMMARY\n")
cat("============================================================\n\n")


for (ph in unique(all_results$phenotype)) {

  cat("\n----------------------------------------\n")
  cat(ph, "\n")
  cat("----------------------------------------\n")

  tmp <- all_results[
    phenotype == ph &
    FDR < 0.05
  ]

  if (nrow(tmp) == 0) {

    cat("No cancer type significant after BH-FDR < 0.05.\n")

  } else {

    print(
      tmp[
        order(FDR, fisher_p),
        .(
          cancer_type,
          total_primary_patients,
          mutated_patients,
          prevalence_percent,
          odds_ratio,
          lower_95CI,
          upper_95CI,
          fisher_p,
          FDR
        )
      ]
    )
  }
}


# ------------------------------------------------------------
# TOP RESULTS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TOP 10 BY FDR — ALL PHENOTYPES\n")
cat("============================================================\n\n")

print(
  all_results[
    order(FDR, fisher_p)
  ][
    1:min(10, .N),
    .(
      phenotype,
      cancer_type,
      total_primary_patients,
      mutated_patients,
      prevalence_percent,
      odds_ratio,
      lower_95CI,
      upper_95CI,
      fisher_p,
      FDR
    )
  ]
)


# ------------------------------------------------------------
# FINAL QC
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("13A FINAL QC\n")
cat("============================================================\n")

cat(
  "Primary patients:",
  uniqueN(patient_cancer$patient_id),
  "\n"
)

cat(
  "Any CLOCK16 patients:",
  length(any_patients),
  "\n"
)

cat(
  "Nonsynonymous CLOCK16 patients:",
  length(nonsyn_patients),
  "\n"
)

cat(
  "Damaging CLOCK16 patients:",
  length(damaging_patients),
  "\n"
)

cat(
  "Cancer types tested:",
  uniqueN(patient_cancer$cancer_type),
  "\n"
)

cat("\n13A COMPLETED SUCCESSFULLY.\n")
cat("============================================================\n\n")

