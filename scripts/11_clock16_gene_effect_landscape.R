
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
# LOAD
# ============================================================

clock <- fread(clock_file)
primary <- fread(primary_file)

genes <- fread(
  genes_file,
  header = FALSE,
  col.names = "gene"
)$gene

N <- uniqueN(primary$patient_id)

# Patient ID
clock[
  ,
  patient_id := sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample
  )
]

# Primary only
clock_primary <- clock[
  substr(sample, 14, 15) == "01"
]

# Only the 16 predefined CLOCK genes
clock_primary <- clock_primary[
  gene %in% genes
]

# ============================================================
# EFFECT CLASSIFICATION
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

# ============================================================
# 1. GENE-LEVEL LANDSCAPE
# ============================================================

gene_landscape <- clock_primary[
  ,
  .(
    Mutation_Events = .N,
    Mutated_Patients = uniqueN(patient_id),
    Nonsynonymous_Events =
      sum(effect_class == "nonsynonymous"),
    Nonsynonymous_Patients =
      uniqueN(
        patient_id[effect_class == "nonsynonymous"]
      ),
    Silent_Events =
      sum(effect_class == "silent"),
    Silent_Patients =
      uniqueN(
        patient_id[effect_class == "silent"]
      ),
    Other_Events =
      sum(effect_class == "other"),
    Other_Patients =
      uniqueN(
        patient_id[effect_class == "other"]
      )
  ),
  by = gene
]

# Ensure all 16 genes appear
gene_landscape <- merge(
  data.table(gene = genes),
  gene_landscape,
  by = "gene",
  all.x = TRUE
)

# Replace missing counts with zero
count_cols <- setdiff(
  names(gene_landscape),
  "gene"
)

for (cc in count_cols) {
  set(
    gene_landscape,
    which(is.na(gene_landscape[[cc]])),
    cc,
    0L
  )
}

gene_landscape[
  ,
  Mutated_Prevalence_Percent :=
    100 * Mutated_Patients / N
]

gene_landscape[
  ,
  Nonsynonymous_Prevalence_Percent :=
    100 * Nonsynonymous_Patients / N
]

setorder(
  gene_landscape,
  -Nonsynonymous_Patients,
  -Mutated_Patients
)

cat("\n============================================================\n")
cat(" CLOCK16 GENE-LEVEL LANDSCAPE\n")
cat("============================================================\n\n")

print(gene_landscape)

# ============================================================
# 2. GENE × EFFECT
# ============================================================

gene_effect <- clock_primary[
  ,
  .(
    Mutation_Events = .N,
    Unique_Patients = uniqueN(patient_id)
  ),
  by = .(gene, effect)
]

setorder(
  gene_effect,
  gene,
  -Mutation_Events
)

cat("\n============================================================\n")
cat(" CLOCK16 GENE × EFFECT\n")
cat("============================================================\n\n")

print(gene_effect)

# ============================================================
# 3. GENE × EFFECT CLASS
# ============================================================

gene_effect_class <- clock_primary[
  ,
  .(
    Mutation_Events = .N,
    Unique_Patients = uniqueN(patient_id)
  ),
  by = .(gene, effect_class)
]

setorder(
  gene_effect_class,
  gene,
  -Mutation_Events
)

cat("\n============================================================\n")
cat(" CLOCK16 GENE × EFFECT CLASS\n")
cat("============================================================\n\n")

print(gene_effect_class)

# ============================================================
# 4. MULTIPLE MUTATIONS WITHIN SAME GENE/PATIENT
# ============================================================

patient_gene_events <- clock_primary[
  ,
  .(
    Mutation_Events = .N
  ),
  by = .(patient_id, gene)
]

multiple_same_gene <- patient_gene_events[
  Mutation_Events > 1
]

cat("\n============================================================\n")
cat(" MULTIPLE EVENTS WITHIN SAME PATIENT-GENE\n")
cat("============================================================\n\n")

cat(
  "Patient-gene combinations with >1 event:",
  nrow(multiple_same_gene),
  "\n"
)

cat(
  "Patients affected:",
  uniqueN(multiple_same_gene$patient_id),
  "\n"
)

print(
  multiple_same_gene[
    ,
    .N,
    by = gene
  ][order(-N)]
)

# ============================================================
# 5. AMINO ACID CHANGE QC
# ============================================================

nonsyn <- clock_primary[
  effect_class == "nonsynonymous"
]

cat("\n============================================================\n")
cat(" AMINO ACID CHANGE QC\n")
cat("============================================================\n\n")

cat(
  "Nonsynonymous events:",
  nrow(nonsyn),
  "\n"
)

cat(
  "Events with Amino_Acid_Change:",
  sum(
    !is.na(nonsyn$Amino_Acid_Change) &
    nonsyn$Amino_Acid_Change != ""
  ),
  "\n"
)

cat(
  "Events without Amino_Acid_Change:",
  sum(
    is.na(nonsyn$Amino_Acid_Change) |
    nonsyn$Amino_Acid_Change == ""
  ),
  "\n"
)

# ============================================================
# 6. SIFT QC
# ============================================================

cat("\n============================================================\n")
cat(" SIFT QC\n")
cat("============================================================\n\n")

if ("SIFT" %in% names(nonsyn)) {

  sift_nonempty <- nonsyn[
    !is.na(SIFT) &
    SIFT != ""
  ]

  cat(
    "Nonsynonymous events with SIFT:",
    nrow(sift_nonempty),
    "/",
    nrow(nonsyn),
    "\n"
  )

  if (nrow(sift_nonempty) > 0) {
    print(
      sift_nonempty[
        ,
        .N,
        by = SIFT
      ][order(-N)]
    )
  }

} else {

  cat("SIFT column not available.\n")

}

# ============================================================
# 7. POLYPHEN QC
# ============================================================

cat("\n============================================================\n")
cat(" PolyPhen QC\n")
cat("============================================================\n\n")

if ("PolyPhen" %in% names(nonsyn)) {

  pp_nonempty <- nonsyn[
    !is.na(PolyPhen) &
    PolyPhen != ""
  ]

  cat(
    "Nonsynonymous events with PolyPhen:",
    nrow(pp_nonempty),
    "/",
    nrow(nonsyn),
    "\n"
  )

  if (nrow(pp_nonempty) > 0) {
    print(
      pp_nonempty[
        ,
        .N,
        by = PolyPhen
      ][order(-N)]
    )
  }

} else {

  cat("PolyPhen column not available.\n")

}

# ============================================================
# 8. TOP AMINO ACID CHANGES
# ============================================================

if ("Amino_Acid_Change" %in% names(nonsyn)) {

  aa <- nonsyn[
    !is.na(Amino_Acid_Change) &
    Amino_Acid_Change != ""
  ]

  aa_summary <- aa[
    ,
    .(
      Mutation_Events = .N,
      Unique_Patients = uniqueN(patient_id)
    ),
    by = .(gene, Amino_Acid_Change)
  ][order(-Mutation_Events)]

} else {

  aa_summary <- data.table()

}

# ============================================================
# 9. SAVE
# ============================================================

fwrite(
  gene_landscape,
  file.path(
    outdir,
    "REAL_CLOCK16_gene_effect_landscape.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_effect,
  file.path(
    outdir,
    "REAL_CLOCK16_gene_effect.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_effect_class,
  file.path(
    outdir,
    "REAL_CLOCK16_gene_effect_class.tsv"
  ),
  sep = "\t"
)

fwrite(
  patient_gene_events,
  file.path(
    outdir,
    "REAL_CLOCK16_patient_gene_event_counts.tsv"
  ),
  sep = "\t"
)

fwrite(
  multiple_same_gene,
  file.path(
    outdir,
    "REAL_CLOCK16_multiple_same_gene_events.tsv"
  ),
  sep = "\t"
)

fwrite(
  aa_summary,
  file.path(
    outdir,
    "REAL_CLOCK16_amino_acid_changes.tsv"
  ),
  sep = "\t"
)

cat("\n============================================================\n")
cat(" FILES SAVED\n")
cat("============================================================\n")

cat(
  "REAL_CLOCK16_gene_effect_landscape.tsv\n",
  "REAL_CLOCK16_gene_effect.tsv\n",
  "REAL_CLOCK16_gene_effect_class.tsv\n",
  "REAL_CLOCK16_patient_gene_event_counts.tsv\n",
  "REAL_CLOCK16_multiple_same_gene_events.tsv\n",
  "REAL_CLOCK16_amino_acid_changes.tsv\n"
)

