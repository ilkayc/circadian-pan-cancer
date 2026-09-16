
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

outdir <- file.path(project, "results/clock16")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genes <- fread(
  file.path(project, "data/clock_genes_16.txt"),
  header = FALSE,
  col.names = "gene"
)$gene

clock <- fread(clock_file)
primary <- fread(primary_file)

cat("\n=== PRIMARY CLOCK16 MUTATION QC ===\n")

cat("CLOCK16 total events:", nrow(clock), "\n")
cat("CLOCK16 total samples:", uniqueN(clock$sample), "\n")

# Patient ID from TCGA barcode
clock[
  ,
  patient_id := sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample
  )
]

# Keep only primary tumor samples
clock_primary <- clock[
  substr(sample, 14, 15) == "01"
]

cat(
  "Primary CLOCK16 mutation events:",
  nrow(clock_primary),
  "\n"
)

cat(
  "Primary CLOCK16 samples:",
  uniqueN(clock_primary$sample),
  "\n"
)

cat(
  "Primary CLOCK16 patients:",
  uniqueN(clock_primary$patient_id),
  "\n"
)

# ------------------------------------------------------------
# Patient-gene combinations
# ------------------------------------------------------------

patient_gene <- unique(
  clock_primary[
    ,
    .(patient_id, gene)
  ]
)

# Number of CLOCK genes mutated per patient

patient_summary <- patient_gene[
  ,
  .(
    CLOCK_genes_mutated = uniqueN(gene)
  ),
  by = patient_id
]

cat("\n=== NUMBER OF CLOCK GENES MUTATED PER PATIENT ===\n")

print(
  patient_summary[
    ,
    .N,
    by = CLOCK_genes_mutated
  ][order(CLOCK_genes_mutated)]
)

# ------------------------------------------------------------
# Effect classification
# ------------------------------------------------------------

clock_primary[
  ,
  effect_class := fifelse(
    effect %in% c(
      "Missense_Mutation",
      "Nonsense_Mutation",
      "Frame_Shift_Del",
      "Frame_Shift_Ins",
      "Splice_Site",
      "In_Frame_Del",
      "In_Frame_Ins",
      "Nonstop_Mutation"
    ),
    "nonsynonymous",
    fifelse(
      effect == "Silent",
      "silent",
      "other"
    )
  )
]

cat("\n=== EFFECT CLASS ===\n")

print(
  clock_primary[
    ,
    .(
      Mutation_Events = .N,
      Unique_Patients = uniqueN(patient_id)
    ),
    by = effect_class
  ][order(-Mutation_Events)]
)

# ------------------------------------------------------------
# Patient-level effect classification
# ------------------------------------------------------------

patient_effect <- clock_primary[
  ,
  .(
    Has_Nonsynonymous =
      any(effect_class == "nonsynonymous"),
    Has_Silent =
      any(effect_class == "silent"),
    Has_Other =
      any(effect_class == "other")
  ),
  by = patient_id
]

patient_effect[
  ,
  Patient_Class := fcase(
    Has_Nonsynonymous, "nonsynonymous",
    Has_Silent & !Has_Nonsynonymous & !Has_Other,
      "silent_only",
    default = "other_only"
  )
]

cat("\n=== PATIENT-LEVEL MUTATION CLASS ===\n")

print(
  patient_effect[
    ,
    .N,
    by = Patient_Class
  ][order(-N)]
)

# ------------------------------------------------------------
# Multiple-gene patients
# ------------------------------------------------------------

multi_gene <- patient_summary[
  CLOCK_genes_mutated >= 2
]

cat(
  "\nPatients with >=2 CLOCK genes mutated:",
  nrow(multi_gene),
  "\n"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fwrite(
  clock_primary,
  file.path(
    outdir,
    "REAL_primary_CLOCK16_mutations.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  patient_summary,
  file.path(
    outdir,
    "REAL_primary_CLOCK16_patient_gene_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  patient_effect,
  file.path(
    outdir,
    "REAL_primary_CLOCK16_patient_effect_class.tsv"
  ),
  sep = "\t"
)

cat("\n=== FILES SAVED ===\n")

cat(
  file.path(
    outdir,
    "REAL_primary_CLOCK16_mutations.tsv.gz"
  ),
  "\n"
)

