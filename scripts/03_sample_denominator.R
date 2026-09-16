
library(data.table)

project <- path.expand("~/circadian_pan_cancer_Q2")

maf_file <- file.path(
  project,
  "data/mc3.v0.2.8.PUBLIC.xena.gz"
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

cat("Reading sample IDs from MC3...\n")

maf_samples <- fread(
  maf_file,
  select = "sample"
)

maf_samples <- unique(maf_samples)

cat(
  "Unique samples with at least one mutation:",
  nrow(maf_samples),
  "\n"
)

cat("\nSample ID examples:\n")
print(head(maf_samples, 20))

cat("\nTCGA sample suffix distribution:\n")

maf_samples[
  ,
  sample_type_code := substr(sample, 14, 15)
]

print(
  maf_samples[
    ,
    .N,
    by = sample_type_code
  ][order(-N)]
)

fwrite(
  maf_samples,
  file.path(
    outdir,
    "mc3_unique_samples.tsv"
  ),
  sep = "\t"
)

cat(
  "\nSaved:",
  file.path(
    outdir,
    "mc3_unique_samples.tsv"
  ),
    "\n"
)

