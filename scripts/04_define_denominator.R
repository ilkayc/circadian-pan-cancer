
library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

phenotype_file <- file.path(
  project,
  "data/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
)

maf_file <- file.path(
  project,
  "data/mc3.v0.2.8.PUBLIC.xena.gz"
)

clock_file <- file.path(
  project,
  "results/clock16/clock16_mutations.tsv.gz"
)

outdir <- file.path(
  project,
  "results",
  "denominator"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("=== DEFINING TCGA PRIMARY TUMOR DENOMINATOR ===\n\n")

# ---------------------------------------------------------
# 1. PHENOTYPE
# ---------------------------------------------------------

cat("Reading phenotype data...\n")

pheno <- fread(phenotype_file)

# Rename the disease column safely
setnames(
  pheno,
  "_primary_disease",
  "cancer_type"
)

cat(
  "Total phenotype samples:",
  nrow(pheno),
  "\n"
)

cat("\nSample types:\n")

sample_type_counts <- pheno[
  ,
  .N,
  by = .(
    sample_type_id,
    sample_type
  )
][order(-N)]

print(sample_type_counts)

# ---------------------------------------------------------
# 2. PRIMARY TUMORS
# ---------------------------------------------------------

primary <- unique(
  pheno[
    sample_type_id == 1,
    .(
      sample,
      sample_type_id,
      sample_type,
      cancer_type
    )
  ]
)

cat(
  "\nPrimary Tumor samples:",
  nrow(primary),
  "\n"
)

cat("\nCancer-type distribution:\n")

primary_cancer_counts <- primary[
  ,
  .N,
  by = cancer_type
][order(-N)]

print(
  head(primary_cancer_counts, 30)
)

# ---------------------------------------------------------
# 3. MC3 SAMPLES
# ---------------------------------------------------------

cat("\nReading MC3 sample IDs...\n")

mc3_samples <- unique(
  fread(
    maf_file,
    select = "sample"
  )
)

cat(
  "MC3 mutation-bearing samples:",
  nrow(mc3_samples),
  "\n"
)

# ---------------------------------------------------------
# 4. PRIMARY TUMOR ∩ MC3
# ---------------------------------------------------------

primary_mc3 <- merge(
  primary,
  mc3_samples,
  by = "sample"
)

cat(
  "\nPrimary tumors represented in MC3:",
  nrow(primary_mc3),
  "\n"
)

cat(
  "Primary tumors without MC3 mutation record:",
  nrow(primary) - nrow(primary_mc3),
  "\n"
)

# ---------------------------------------------------------
# 5. CLOCK16 MUTATIONS
# ---------------------------------------------------------

cat("\nReading clock16 mutations...\n")

clock16 <- fread(
  clock_file,
  select = c(
    "sample",
    "gene"
  )
)

clock16_samples <- unique(
  clock16[, .(sample)]
)

cat(
  "Tumors with >=1 clock-gene mutation:",
  nrow(clock16_samples),
  "\n"
)

# ---------------------------------------------------------
# 6. PRIMARY TUMORS WITH CLOCK MUTATION
# ---------------------------------------------------------

primary_clock16 <- merge(
  primary,
  clock16_samples,
  by = "sample"
)

cat(
  "\nPrimary tumors with >=1 clock-gene mutation:",
  nrow(primary_clock16),
  "\n"
)

# ---------------------------------------------------------
# 7. OVERALL PREVALENCE
# ---------------------------------------------------------

overall_prevalence <- data.table(
  Metric = c(
    "All Primary Tumor samples",
    "Primary tumors represented in MC3",
    "Primary tumors with >=1 clock-gene mutation"
  ),
  N = c(
    nrow(primary),
    nrow(primary_mc3),
    nrow(primary_clock16)
  )
)

overall_prevalence[
  ,
  Percent_of_all_primary :=
    round(
      100 * N / nrow(primary),
      3
    )
]

print(overall_prevalence)

# ---------------------------------------------------------
# 8. CANCER-SPECIFIC PREVALENCE
# ---------------------------------------------------------

cancer_denominator <- primary[
  ,
  .(
    All_Primary_Tumors = .N
  ),
  by = cancer_type
]

cancer_clock16 <- primary_clock16[
  ,
  .(
    Clock16_Mutated_Tumors = .N
  ),
  by = cancer_type
]

cancer_prevalence <- merge(
  cancer_denominator,
  cancer_clock16,
  by = "cancer_type",
  all.x = TRUE
)

cancer_prevalence[
  is.na(Clock16_Mutated_Tumors),
  Clock16_Mutated_Tumors := 0L
]

cancer_prevalence[
  ,
  Clock16_Prevalence_Percent :=
    round(
      100 *
      Clock16_Mutated_Tumors /
      All_Primary_Tumors,
      3
    )
]

setorder(
  cancer_prevalence,
  -Clock16_Mutated_Tumors
)

# ---------------------------------------------------------
# 9. SAVE
# ---------------------------------------------------------

fwrite(
  primary,
  file.path(
    outdir,
    "all_primary_tumor_samples.tsv"
  ),
  sep = "\t"
)

fwrite(
  primary_mc3,
  file.path(
    outdir,
    "primary_tumors_in_MC3.tsv"
  ),
  sep = "\t"
)

fwrite(
  primary_clock16,
  file.path(
    outdir,
    "primary_tumors_clock16_mutated.tsv"
  ),
  sep = "\t"
)

fwrite(
  overall_prevalence,
  file.path(
    outdir,
    "overall_denominator_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  cancer_prevalence,
  file.path(
    outdir,
    "cancer_type_denominator_prevalence.tsv"
  ),
  sep = "\t"
)

cat("\n=== DENOMINATOR ANALYSIS COMPLETE ===\n")

cat(
  "\nAll Primary Tumors:",
  nrow(primary),
  "\n"
)

cat(
  "Primary Tumors in MC3:",
  nrow(primary_mc3),
  "\n"
)

cat(
  "Primary Tumors with Clock Mutation:",
  nrow(primary_clock16),
  "\n"
)

