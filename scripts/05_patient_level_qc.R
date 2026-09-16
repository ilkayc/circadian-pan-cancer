
library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

primary_file <- file.path(
  project,
  "results/denominator/all_primary_tumor_samples.tsv"
)

clock_file <- file.path(
  project,
  "results/denominator/primary_tumors_clock16_mutated.tsv"
)

outdir <- file.path(
  project,
  "results/denominator"
)

primary <- fread(primary_file)
clock_primary <- fread(clock_file)

# TCGA patient identifier = first 3 components
# Example: TCGA-02-0003-01 -> TCGA-02-0003

primary[
  ,
  patient_id := sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample
  )
]

clock_primary[
  ,
  patient_id := sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample
  )
]

cat("=== PATIENT-LEVEL QC ===\n\n")

cat(
  "Primary tumor samples:",
  uniqueN(primary$sample),
  "\n"
)

cat(
  "Primary tumor patients:",
  uniqueN(primary$patient_id),
  "\n"
)

cat(
  "Clock-mutated primary samples:",
  uniqueN(clock_primary$sample),
  "\n"
)

cat(
  "Clock-mutated primary patients:",
  uniqueN(clock_primary$patient_id),
  "\n"
)

# Number of primary samples per patient
sample_per_patient <- primary[
  ,
  .(
    n_primary_samples = .N
  ),
  by = patient_id
]

cat("\nPrimary samples per patient:\n")
print(
  sample_per_patient[
    ,
    .N,
    by = n_primary_samples
  ][order(n_primary_samples)]
)

# Patients with multiple primary samples
multi_primary <- sample_per_patient[
  n_primary_samples > 1
]

cat(
  "\nPatients with >1 primary tumor sample:",
  nrow(multi_primary),
  "\n"
)

if (nrow(multi_primary) > 0) {
  print(
    head(multi_primary, 30)
  )
}

# Clock mutation status at patient level
clock_patient <- unique(
  clock_primary[
    ,
    .(
      patient_id
    )
  ]
)

cat(
  "\nPatient-level clock mutation prevalence:",
  round(
    100 *
    nrow(clock_patient) /
    uniqueN(primary$patient_id),
    3
  ),
    "%\n"
)

# Save
fwrite(
  primary,
  file.path(
    outdir,
    "primary_samples_with_patient_id.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_per_patient,
  file.path(
    outdir,
    "primary_sample_count_per_patient.tsv"
  ),
  sep = "\t"
)

fwrite(
  clock_patient,
  file.path(
    outdir,
    "clock_mutated_patients.tsv"
  ),
  sep = "\t"
)

cat("\n=== PATIENT-LEVEL QC COMPLETE ===\n")

