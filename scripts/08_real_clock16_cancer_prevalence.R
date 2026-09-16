library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

primary <- fread(
  file.path(
    project,
    "results/denominator/primary_samples_with_patient_id.tsv"
  )
)

clock <- fread(
  file.path(
    project,
    "results/clock16/clock16_mutations.tsv.gz"
  )
)

# Correct TCGA patient ID
clock[
  ,
  patient_id := sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample
  )
]

# Primary CLOCK16 patients only
clock_primary <- unique(
  clock[
    substr(sample, 14, 15) == "01",
    .(patient_id)
  ]
)

# Keep only patients in primary denominator
clock_primary <- clock_primary[
  patient_id %in% primary$patient_id
]

# Patient-level cancer table
patient_cancer <- unique(
  primary[
    ,
    .(patient_id, cancer_type)
  ]
)

# Denominator by cancer type
denom <- patient_cancer[
  ,
  .(Total_Primary_Patients = uniqueN(patient_id)),
  by = cancer_type
]

# Mutated patients by cancer type
mut <- patient_cancer[
  patient_id %in% clock_primary$patient_id,
  .(CLOCK16_Mutated_Patients = uniqueN(patient_id)),
  by = cancer_type
]

# Merge
res <- merge(
  denom,
  mut,
  by = "cancer_type",
  all.x = TRUE
)

res[
  is.na(CLOCK16_Mutated_Patients),
  CLOCK16_Mutated_Patients := 0L
]

# Prevalence and exact 95% CI
res[
  ,
  Prevalence_Percent :=
    100 * CLOCK16_Mutated_Patients / Total_Primary_Patients
]

res[
  ,
  Lower_95CI :=
    mapply(
      function(x, n) {
        100 * binom.test(x, n)$conf.int[1]
      },
      CLOCK16_Mutated_Patients,
      Total_Primary_Patients
    )
]

res[
  ,
  Upper_95CI :=
    mapply(
      function(x, n) {
        100 * binom.test(x, n)$conf.int[2]
      },
      CLOCK16_Mutated_Patients,
      Total_Primary_Patients
    )
]

# Round
res[
  ,
  `:=`(
    Prevalence_Percent = round(Prevalence_Percent, 3),
    Lower_95CI = round(Lower_95CI, 3),
    Upper_95CI = round(Upper_95CI, 3)
  )
]

# Order by prevalence
setorder(
  res,
  -Prevalence_Percent
)

cat("\n=== REAL CLOCK16 CANCER-TYPE PREVALENCE ===\n\n")
print(res, nrows = nrow(res))

# Save
fwrite(
  res,
  file.path(
    project,
    "results/clock16/REAL_clock16_cancer_prevalence.tsv"
  ),
  sep = "\t"
)

cat("\nSaved:\n")
cat(
  file.path(
    project,
    "results/clock16/REAL_clock16_cancer_prevalence.tsv"
  ),
  "\n"
)
