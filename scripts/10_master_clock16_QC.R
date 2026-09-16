
library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

clock_file <- file.path(
  project,
  "results/clock16/clock16_mutations.tsv.gz"
)

primary_file <- file.path(
  project,
  "results/denominator/primary_samples_with_patient_id.tsv"
)

genes_file <- file.path(
  project,
  "data/clock_genes_16.txt"
)

outdir <- file.path(
  project,
  "results/clock16"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD
# ============================================================

clock <- fread(clock_file)
primary <- fread(primary_file)
genes <- fread(
  genes_file,
  header = FALSE,
  col.names = "gene"
)$gene

cat("\n============================================================\n")
cat(" MASTER CLOCK16 QC\n")
cat("============================================================\n\n")

cat("CLOCK16 genes expected:", length(genes), "\n")
cat("CLOCK16 genes observed:", uniqueN(clock$gene), "\n")

missing_genes <- setdiff(genes, unique(clock$gene))
extra_genes   <- setdiff(unique(clock$gene), genes)

cat("Missing expected genes:", length(missing_genes), "\n")
cat("Unexpected genes:", length(extra_genes), "\n")

if (length(missing_genes) > 0) {
  print(missing_genes)
}

if (length(extra_genes) > 0) {
  print(extra_genes)
}

# ============================================================
# 2. PRIMARY DENOMINATOR QC
# ============================================================

N <- uniqueN(primary$patient_id)

cat("\n--- PRIMARY DENOMINATOR ---\n")
cat("Primary samples :", uniqueN(primary$sample), "\n")
cat("Primary patients:", N, "\n")

if (uniqueN(primary$sample) != N) {
  cat("WARNING: multiple primary samples per patient detected.\n")
} else {
  cat("PASS: one primary sample per patient.\n")
}

# ============================================================
# 3. ROBUST PATIENT ID CREATION
# ============================================================

clock[
  ,
  patient_id := sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample
  )
]

# ============================================================
# 4. SAMPLE TYPE
# ============================================================

clock[
  ,
  sample_type_code := substr(sample, 14, 15)
]

cat("\n--- CLOCK16 SAMPLE TYPES ---\n")

sample_qc <- clock[
  ,
  .(
    Events = .N,
    Samples = uniqueN(sample),
    Patients = uniqueN(patient_id)
  ),
  by = sample_type_code
][order(sample_type_code)]

print(sample_qc)

# ============================================================
# 5. PRIMARY CLOCK16
# ============================================================

clock_primary <- clock[
  sample_type_code == "01"
]

cat("\n--- PRIMARY CLOCK16 ---\n")

cat("Events   :", nrow(clock_primary), "\n")
cat("Samples  :", uniqueN(clock_primary$sample), "\n")
cat("Patients :", uniqueN(clock_primary$patient_id), "\n")

# ============================================================
# 6. PRIMARY IDs MUST EXIST IN DENOMINATOR
# ============================================================

primary_patients <- unique(primary$patient_id)

clock_primary_not_denominator <- setdiff(
  unique(clock_primary$patient_id),
  primary_patients
)

cat("\n--- DENOMINATOR MEMBERSHIP ---\n")

cat(
  "Primary CLOCK16 patients NOT in denominator:",
  length(clock_primary_not_denominator),
  "\n"
)

if (length(clock_primary_not_denominator) > 0) {
  print(head(clock_primary_not_denominator, 20))
}

# ============================================================
# 7. PRIMARY SAMPLE IDs MUST MATCH DENOMINATOR
# ============================================================

primary_samples <- unique(primary$sample)

clock_primary_not_in_primary_samples <- setdiff(
  unique(clock_primary$sample),
  primary_samples
)

cat(
  "Primary CLOCK16 samples NOT in denominator:",
  length(clock_primary_not_in_primary_samples),
  "\n"
)

# ============================================================
# 8. PATIENT-GENE UNIQUENESS
# ============================================================

patient_gene <- unique(
  clock_primary[
    ,
    .(patient_id, gene)
  ]
)

cat("\n--- PATIENT-GENE QC ---\n")

cat(
  "Raw primary CLOCK16 events:",
  nrow(clock_primary),
  "\n"
)

cat(
  "Unique patient-gene combinations:",
  nrow(patient_gene),
  "\n"
)

# ============================================================
# 9. GENES MUTATED PER PATIENT
# ============================================================

patient_gene_count <- patient_gene[
  ,
  .(
    CLOCK_genes_mutated = uniqueN(gene)
  ),
  by = patient_id
]

gene_count_distribution <- patient_gene_count[
  ,
  .N,
  by = CLOCK_genes_mutated
][order(CLOCK_genes_mutated)]

cat("\n--- NUMBER OF CLOCK GENES PER PATIENT ---\n")
print(gene_count_distribution)

cat(
  "\nPatients with >=2 CLOCK genes:",
  patient_gene_count[
    CLOCK_genes_mutated >= 2,
    .N
  ],
  "\n"
)

# ============================================================
# 10. EFFECT QC
# ============================================================

cat("\n--- EFFECT QC ---\n")

effect_qc <- clock_primary[
  ,
  .(
    Mutation_Events = .N,
    Unique_Patients = uniqueN(patient_id)
  ),
  by = effect
][order(-Mutation_Events)]

print(effect_qc)

# ============================================================
# 11. EFFECT CLASSIFICATION
# ============================================================

nonsynonymous_effects <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Nonstop_Mutation"
)

clock_primary[
  ,
  effect_class := fcase(
    effect %in% nonsynonymous_effects,
    "nonsynonymous",
    effect == "Silent",
    "silent",
    default = "other"
  )
]

effect_class_qc <- clock_primary[
  ,
  .(
    Mutation_Events = .N,
    Unique_Patients = uniqueN(patient_id)
  ),
  by = effect_class
][order(-Mutation_Events)]

cat("\n--- EFFECT CLASS ---\n")
print(effect_class_qc)

# ============================================================
# 12. PATIENT-LEVEL CLASS
# ============================================================

patient_effect <- clock_primary[
  ,
  .(
    Has_Nonsynonymous = any(effect_class == "nonsynonymous"),
    Has_Silent = any(effect_class == "silent"),
    Has_Other = any(effect_class == "other")
  ),
  by = patient_id
]

patient_effect[
  ,
  Patient_Class := fcase(
    Has_Nonsynonymous,
    "nonsynonymous",
    Has_Silent & !Has_Nonsynonymous & !Has_Other,
    "silent_only",
    default = "other_only"
  )
]

patient_class_qc <- patient_effect[
  ,
  .N,
  by = Patient_Class
][order(-N)]

cat("\n--- PATIENT-LEVEL CLASS ---\n")
print(patient_class_qc)

# ============================================================
# 13. OVERALL PRIMARY PANEL PREVALENCE
# ============================================================

M <- uniqueN(clock_primary$patient_id)

overall <- data.table(
  Primary_Patients = N,
  CLOCK16_Primary_Patients = M,
  CLOCK16_Prevalence_Percent = 100 * M / N
)

cat("\n--- OVERALL PRIMARY PREVALENCE ---\n")
print(overall)

# ============================================================
# 14. NONSYNONYMOUS PRIMARY PREVALENCE
# ============================================================

nonsyn_patients <- unique(
  clock_primary[
    effect_class == "nonsynonymous",
    patient_id
  ]
)

M_nonsyn <- length(nonsyn_patients)

bt_nonsyn <- binom.test(
  M_nonsyn,
  N
)

nonsyn_prev <- data.table(
  Primary_Patients = N,
  Nonsynonymous_CLOCK16_Patients = M_nonsyn,
  Prevalence_Percent = 100 * M_nonsyn / N,
  Lower_95CI = 100 * bt_nonsyn$conf.int[1],
  Upper_95CI = 100 * bt_nonsyn$conf.int[2]
)

cat("\n--- NONSYNONYMOUS PRIMARY PREVALENCE ---\n")
print(nonsyn_prev)

# ============================================================
# 15. GENE-LEVEL PRIMARY PREVALENCE
# ============================================================

gene_prev <- patient_gene[
  ,
  .(
    Mutated_Patients = uniqueN(patient_id)
  ),
  by = gene
]

gene_prev[
  ,
  Prevalence_Percent := 100 * Mutated_Patients / N
]

gene_prev[
  ,
  Mutation_Events := clock_primary[
    .SD,
    on = .(gene),
    .N,
    by = .EACHI
  ]$N
]

setorder(gene_prev, -Mutated_Patients)

cat("\n--- GENE-LEVEL PRIMARY PREVALENCE ---\n")
print(gene_prev)

# ============================================================
# 16. FINAL ASSERTIONS
# ============================================================

cat("\n============================================================\n")
cat(" FINAL ASSERTIONS\n")
cat("============================================================\n")

stopifnot(
  length(genes) == 16
)

stopifnot(
  uniqueN(primary$patient_id) == 10593
)

stopifnot(
  uniqueN(clock_primary$sample) ==
  uniqueN(clock_primary$patient_id)
)

stopifnot(
  length(clock_primary_not_denominator) == 0
)

stopifnot(
  length(clock_primary_not_in_primary_samples) == 0
)

stopifnot(
  M == 1165
)

cat("\nALL CORE ASSERTIONS PASSED.\n")

# ============================================================
# 17. SAVE MASTER QC
# ============================================================

fwrite(
  sample_qc,
  file.path(outdir, "REAL_MASTER_sample_type_QC.tsv"),
  sep = "\t"
)

fwrite(
  effect_qc,
  file.path(outdir, "REAL_MASTER_effect_QC.tsv"),
  sep = "\t"
)

fwrite(
  effect_class_qc,
  file.path(outdir, "REAL_MASTER_effect_class_QC.tsv"),
  sep = "\t"
)

fwrite(
  patient_class_qc,
  file.path(outdir, "REAL_MASTER_patient_class_QC.tsv"),
  sep = "\t"
)

fwrite(
  gene_count_distribution,
  file.path(outdir, "REAL_MASTER_genes_per_patient.tsv"),
  sep = "\t"
)

fwrite(
  overall,
  file.path(outdir, "REAL_MASTER_overall_prevalence.tsv"),
  sep = "\t"
)

fwrite(
  nonsyn_prev,
  file.path(outdir, "REAL_MASTER_nonsynonymous_prevalence.tsv"),
  sep = "\t"
)

fwrite(
  gene_prev,
  file.path(outdir, "REAL_MASTER_gene_prevalence.tsv"),
  sep = "\t"
)

cat("\nMASTER QC FILES SAVED.\n")

