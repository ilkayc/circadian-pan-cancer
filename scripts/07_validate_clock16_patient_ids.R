library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

primary_file <- file.path(
  project,
  "results/denominator/primary_samples_with_patient_id.tsv"
)

mutation_file <- file.path(
  project,
  "results/clock16/clock16_mutations.tsv.gz"
)

primary <- fread(primary_file)
clock <- fread(mutation_file)

cat("\n=== PATIENT ID VALIDATION ===\n")

# ------------------------------------------------------------
# PRIMARY DATA
# ------------------------------------------------------------

cat("\nPrimary columns:\n")
print(names(primary))

cat("\nPrimary patients:\n")
cat("Unique patient_id:", uniqueN(primary$patient_id), "\n")
cat("Unique samples:", uniqueN(primary$sample), "\n")


# ------------------------------------------------------------
# CLOCK16
# ------------------------------------------------------------

clock[
  ,
  patient_id_correct :=
    sub("^((TCGA-[^-]+-[^-]+)).*$", "\\1", sample)
]

clock[
  ,
  sample_type_code := substr(sample, 14, 15)
]

cat("\nCLOCK16:\n")
cat("Unique samples:", uniqueN(clock$sample), "\n")
cat(
  "Unique patients:",
  uniqueN(clock$patient_id_correct),
  "\n"
)

cat("\nSample type distribution:\n")

print(
  clock[
    ,
    .(
      Unique_Samples = uniqueN(sample),
      Unique_Patients = uniqueN(patient_id_correct)
    ),
    by = sample_type_code
  ][order(sample_type_code)]
)


# ------------------------------------------------------------
# PRIMARY CLOCK16 ONLY
# ------------------------------------------------------------

clock_primary <- clock[
  sample_type_code == "01"
]

cat("\n=== PRIMARY CLOCK16 ===\n")

cat(
  "Unique primary samples:",
  uniqueN(clock_primary$sample),
  "\n"
)

cat(
  "Unique primary patients:",
  uniqueN(clock_primary$patient_id_correct),
  "\n"
)


# ------------------------------------------------------------
# COMPARE WITH PRIMARY DENOMINATOR
# ------------------------------------------------------------

primary_ids <- unique(primary$patient_id)

clock_primary_ids <- unique(clock_primary$patient_id_correct)

cat(
  "\nCLOCK16 primary patients IN denominator:",
  length(intersect(clock_primary_ids, primary_ids)),
  "\n"
)

cat(
  "CLOCK16 primary patients NOT in denominator:",
  length(setdiff(clock_primary_ids, primary_ids)),
  "\n"
)

cat(
  "Denominator patients WITHOUT CLOCK16 mutation:",
  length(setdiff(primary_ids, clock_primary_ids)),
  "\n"
)


# ------------------------------------------------------------
# OVERALL PANEL PREVALENCE
# ------------------------------------------------------------

N <- length(primary_ids)

M <- length(intersect(clock_primary_ids, primary_ids))

cat("\n=== FINAL PRIMARY PANEL QC ===\n")

cat("Primary denominator:", N, "\n")
cat("CLOCK16-mutated primary patients:", M, "\n")
cat(
  "Prevalence:",
  round(100 * M / N, 3),
  "%\n"
)


# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

qc <- data.table(
  Primary_Denominator = N,
  CLOCK16_Primary_Patients = M,
  CLOCK16_Primary_Prevalence_Percent =
    round(100 * M / N, 3),
  CLOCK16_Primary_Samples =
    uniqueN(clock_primary$sample),
  CLOCK16_All_Samples =
    uniqueN(clock$sample),
  CLOCK16_All_Patients =
    uniqueN(clock$patient_id_correct)
)

fwrite(
  qc,
  file.path(
    project,
    "results/clock16/REAL_clock16_patient_QC.tsv"
  ),
  sep = "\t"
)

cat("\nQC file saved.\n")
