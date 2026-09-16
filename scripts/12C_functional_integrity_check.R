
library(data.table)

cat("\n")
cat("============================================================\n")
cat("12C — CLOCK16 FUNCTIONAL INTEGRITY CHECK\n")
cat("============================================================\n\n")

base_dir <- "results/clock16"

input_file <- file.path(
  base_dir,
  "REAL_CLOCK16_nonsynonymous_functional_classification.tsv.gz"
)

if (!file.exists(input_file)) {
  stop("ERROR: 12B classification file not found.")
}

dt <- fread(input_file)

cat("Input file:\n")
cat(input_file, "\n\n")

cat("Rows:", nrow(dt), "\n")
cat("Columns:", ncol(dt), "\n\n")


# ------------------------------------------------------------
# 1. REQUIRED COLUMNS
# ------------------------------------------------------------

required_cols <- c(
  "patient_id",
  "gene",
  "Functional_Class"
)

missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0) {
  stop(
    "ERROR: Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("=== REQUIRED COLUMNS ===\n")
cat("All required columns present: YES\n\n")


# ------------------------------------------------------------
# 2. FUNCTIONAL CLASS COUNTS
# ------------------------------------------------------------

cat("=== FUNCTIONAL CLASS COUNTS ===\n")

class_counts <- dt[, .N, by = Functional_Class][order(-N)]

print(class_counts)

cat("\n")


# ------------------------------------------------------------
# 3. TOTAL NONSYNONYMOUS EVENTS
# ------------------------------------------------------------

cat("=== NONSYNONYMOUS EVENT CHECK ===\n")

n_total <- nrow(dt)

cat("Total events:", n_total, "\n")

if (n_total == 1536) {
  cat("PASS: 1536 nonsynonymous events retained.\n")
} else {
  cat("WARNING: Expected 1536 but found", n_total, "\n")
}

cat("\n")


# ------------------------------------------------------------
# 4. STRONG DAMAGING CHECK
# ------------------------------------------------------------

cat("=== STRONG DAMAGING CHECK ===\n")

strong <- dt[
  Functional_Class == "strong_damaging"
]

strong_events <- nrow(strong)
strong_patients <- uniqueN(strong$patient_id)

cat("Strong damaging events:", strong_events, "\n")
cat("Strong damaging patients:", strong_patients, "\n")

if (strong_events == 401) {
  cat("PASS: 401 strong damaging events.\n")
} else {
  cat("WARNING: Expected 401 strong damaging events.\n")
}

if (strong_patients == 294) {
  cat("PASS: 294 strong damaging patients.\n")
} else {
  cat("WARNING: Expected 294 strong damaging patients.\n")
}

cat("\n")


# ------------------------------------------------------------
# 5. DAMAGING CHECK
# ------------------------------------------------------------

cat("=== DAMAGING CHECK ===\n")

damaging_classes <- c(
  "strong_damaging",
  "damaging"
)

damaging <- dt[
  Functional_Class %in% damaging_classes
]

damaging_events <- nrow(damaging)
damaging_patients <- uniqueN(damaging$patient_id)

cat("Damaging events:", damaging_events, "\n")
cat("Damaging patients:", damaging_patients, "\n")

if (damaging_events == 598) {
  cat("PASS: 598 damaging events.\n")
} else {
  cat("WARNING: Expected 598 damaging events.\n")
}

if (damaging_patients == 395) {
  cat("PASS: 395 damaging patients.\n")
} else {
  cat("WARNING: Expected 395 damaging patients.\n")
}

cat("\n")


# ------------------------------------------------------------
# 6. CLASS SUM CHECK
# ------------------------------------------------------------

cat("=== CLASS SUM CHECK ===\n")

class_sum <- sum(class_counts$N)

cat("Sum of functional classes:", class_sum, "\n")
cat("Total rows:", n_total, "\n")

if (class_sum == n_total) {
  cat("PASS: Functional classes account for every event.\n")
} else {
  cat("FAIL: Functional classes do not account for every event.\n")
}

cat("\n")


# ------------------------------------------------------------
# 7. DUPLICATE EXACT ROW CHECK
# ------------------------------------------------------------

cat("=== EXACT DUPLICATE CHECK ===\n")

duplicate_rows <- sum(duplicated(dt))

cat("Exact duplicate rows:", duplicate_rows, "\n")

if (duplicate_rows == 0) {
  cat("PASS: No exact duplicate rows.\n")
} else {
  cat("WARNING: Exact duplicate rows detected.\n")
}

cat("\n")


# ------------------------------------------------------------
# 8. PATIENT-GENE COMBINATION CHECK
# ------------------------------------------------------------

cat("=== PATIENT-GENE CHECK ===\n")

if ("gene" %in% names(dt)) {

  patient_gene_n <- uniqueN(
    dt[, .(patient_id, gene)]
  )

  cat(
    "Unique patient-gene combinations:",
    patient_gene_n,
    "\n"
  )

} else {

  cat("gene column unavailable.\n")

}

cat("\n")


# ------------------------------------------------------------
# 9. MISSING FUNCTIONAL PREDICTION CHECK
# ------------------------------------------------------------

cat("=== MISSING PREDICTION CHECK ===\n")

missing_prediction <- dt[
  Functional_Class == "prediction_missing"
]

cat(
  "Prediction-missing events:",
  nrow(missing_prediction),
  "\n"
)

cat(
  "Prediction-missing patients:",
  uniqueN(missing_prediction$patient_id),
  "\n"
)

cat(
  "Important: these events were NOT assigned a functional class.\n"
)

cat("\n")


# ------------------------------------------------------------
# 10. IMPOSSIBLE CLASS CHECK
# ------------------------------------------------------------

cat("=== CLASS VALIDITY CHECK ===\n")

allowed_classes <- c(
  "likely_benign",
  "strong_damaging",
  "prediction_missing",
  "damaging",
  "discordant_SIFT_deleterious",
  "discordant_PolyPhen_damaging"
)

unexpected_classes <- setdiff(
  unique(dt$Functional_Class),
  allowed_classes
)

if (length(unexpected_classes) == 0) {

  cat("PASS: No unexpected functional classes.\n")

} else {

  cat(
    "WARNING: Unexpected classes:",
    paste(unexpected_classes, collapse = ", "),
    "\n"
  )

}

cat("\n")


# ------------------------------------------------------------
# 11. FINAL INTEGRITY SUMMARY
# ------------------------------------------------------------

cat("============================================================\n")
cat("12C FINAL INTEGRITY SUMMARY\n")
cat("============================================================\n")

all_ok <- (
  n_total == 1536 &&
  strong_events == 401 &&
  strong_patients == 294 &&
  damaging_events == 598 &&
  damaging_patients == 395 &&
  class_sum == n_total &&
  duplicate_rows == 0 &&
  length(unexpected_classes) == 0
)

if (all_ok) {

  cat("\n")
  cat("12C STATUS: PASS\n")
  cat("\n")
  cat("The 12B functional classification is internally consistent.\n")
  cat("No artificial filling of missing predictions was detected.\n")
  cat("No exact duplicate rows were detected.\n")
  cat("All 1536 nonsynonymous events are accounted for.\n")
  cat("\n")

} else {

  cat("\n")
  cat("12C STATUS: REVIEW REQUIRED\n")
  cat("\n")
  cat("One or more integrity checks did not match expectations.\n")
  cat("DO NOT proceed to downstream association analysis yet.\n")
  cat("\n")

}

cat("============================================================\n\n")

