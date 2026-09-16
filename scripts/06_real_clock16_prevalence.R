library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

mutation_file <- file.path(
  project,
  "results/clock16/clock16_mutations.tsv.gz"
)

primary_file <- file.path(
  project,
  "results/denominator/primary_samples_with_patient_id.tsv"
)

outdir <- file.path(
  project,
  "results/clock16"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genes <- fread(
  file.path(project, "data/clock_genes_16.txt"),
  header = FALSE,
  col.names = "gene"
)$gene

primary <- fread(primary_file)

clock <- fread(mutation_file)

N <- uniqueN(primary$patient_id)

clock[
  ,
  patient_id := sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample
  )
]

clock <- clock[gene %in% genes]

# Unique patient-gene combinations
patient_gene <- unique(
  clock[, .(patient_id, gene)]
)

# Gene-level prevalence
prev <- patient_gene[
  ,
  .(Mutated_Patients = uniqueN(patient_id)),
  by = gene
]

# Ensure all 16 genes are present
prev <- merge(
  data.table(gene = genes),
  prev,
  by = "gene",
  all.x = TRUE
)

prev[
  is.na(Mutated_Patients),
  Mutated_Patients := 0L
]

# Exact binomial 95% CI
ci <- t(
  sapply(
    prev$Mutated_Patients,
    function(x) {
      bt <- binom.test(x, N)
      c(
        Lower_95CI = 100 * bt$conf.int[1],
        Upper_95CI = 100 * bt$conf.int[2]
      )
    }
  )
)

prev[
  ,
  Prevalence_Percent :=
    100 * Mutated_Patients / N
]

prev[
  ,
  Lower_95CI := ci[, 1]
]

prev[
  ,
  Upper_95CI := ci[, 2]
]

# Mutation event counts
events <- clock[
  ,
  .(Mutation_Events = .N),
  by = gene
]

prev <- merge(
  prev,
  events,
  by = "gene",
  all.x = TRUE
)

prev[
  ,
  Prevalence_Percent :=
    round(Prevalence_Percent, 3)
]

prev[
  ,
  Lower_95CI :=
    round(Lower_95CI, 3)
]

prev[
  ,
  Upper_95CI :=
    round(Upper_95CI, 3)
]

# Order by prevalence
setorder(prev, -Mutated_Patients)

cat("\n=== REAL CLOCK16 PREVALENCE ===\n")
cat("Denominator:", N, "unique primary patients\n\n")
print(prev)

# Overall panel
panel_patients <- unique(patient_gene$patient_id)

x <- length(panel_patients)

bt <- binom.test(x, N)

overall <- data.table(
  Total_Primary_Patients = N,
  CLOCK16_Mutated_Patients = x,
  Prevalence_Percent = round(100 * x / N, 3),
  Lower_95CI = round(100 * bt$conf.int[1], 3),
  Upper_95CI = round(100 * bt$conf.int[2], 3)
)

cat("\n=== REAL OVERALL PANEL PREVALENCE ===\n")
print(overall)

fwrite(
  prev,
  file.path(
    outdir,
    "REAL_clock16_gene_prevalence.tsv"
  ),
  sep = "\t"
)

fwrite(
  overall,
  file.path(
    outdir,
    "REAL_clock16_overall_prevalence.tsv"
  ),
  sep = "\t"
)

cat("\nFiles saved.\n")
