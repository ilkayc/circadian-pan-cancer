
library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

maf_file <- file.path(
  project,
  "data/mc3.v0.2.8.PUBLIC.xena.gz"
)

gene_file <- file.path(
  project,
  "data/clock_genes_16.txt"
)

outdir <- file.path(project, "results", "clock16")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clock_genes <- fread(
  gene_file,
  header = FALSE,
  col.names = "gene"
)$gene

cat("16-gene panel:\n")
print(clock_genes)

cat("\nReading MC3 mutation data...\n")

maf <- fread(
  maf_file,
  select = c(
    "sample",
    "chr",
    "start",
    "end",
    "reference",
    "alt",
    "gene",
    "effect",
    "Amino_Acid_Change",
    "DNA_VAF",
    "SIFT",
    "PolyPhen"
  )
)

cat("\nTotal mutation events:", nrow(maf), "\n")

clock16 <- maf[
  gene %in% clock_genes
]

cat("Clock-gene mutation events:", nrow(clock16), "\n")

cat(
  "Unique mutated tumors:",
  uniqueN(clock16$sample),
  "\n"
)

cat(
  "Unique genes detected:",
  uniqueN(clock16$gene),
  "\n"
)

cat("\nGene-level mutation counts:\n")

gene_counts <- clock16[
  ,
  .(
    Mutation_Events = .N,
    Mutated_Tumors = uniqueN(sample),
    Missense = sum(effect == "Missense_Mutation", na.rm = TRUE),
    Nonsense = sum(effect == "Nonsense_Mutation", na.rm = TRUE),
    Frameshift_Del = sum(effect == "Frame_Shift_Del", na.rm = TRUE),
    Frameshift_Ins = sum(effect == "Frame_Shift_Ins", na.rm = TRUE),
    Splice_Site = sum(effect == "Splice_Site", na.rm = TRUE),
    InFrame_Del = sum(effect == "In_Frame_Del", na.rm = TRUE),
    InFrame_Ins = sum(effect == "In_Frame_Ins", na.rm = TRUE),
    Silent = sum(effect == "Silent", na.rm = TRUE)
  ),
  by = gene
]

setorder(gene_counts, -Mutated_Tumors)

print(gene_counts)

fwrite(
  gene_counts,
  file.path(outdir, "clock16_gene_counts.tsv"),
  sep = "\t"
)

cat("\nEffect distribution:\n")

effect_counts <- clock16[
  ,
  .(Mutation_Events = .N),
  by = effect
]

setorder(effect_counts, -Mutation_Events)

print(effect_counts)

fwrite(
  effect_counts,
  file.path(outdir, "clock16_effect_counts.tsv"),
  sep = "\t"
)

cat("\nSaving clock16 mutation table...\n")

fwrite(
  clock16,
  file.path(outdir, "clock16_mutations.tsv.gz"),
  sep = "\t"
)

cat("\nChecking missing genes:\n")

missing_genes <- setdiff(clock_genes, unique(clock16$gene))

print(missing_genes)

writeLines(
  missing_genes,
  file.path(outdir, "clock16_missing_genes.txt")
)

cat("\n=== CLOCK16 EXTRACTION COMPLETE ===\n")

