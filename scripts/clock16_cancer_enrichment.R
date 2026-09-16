
primary_file <- "results/denominator/primary_samples_with_patient_id.tsv"
mut_file <- "results/clock16/clock16_primary_mutations_with_cancer.tsv"

# ============================================================
# PRIMARY DENOMINATOR
# ============================================================

primary <- read.delim(
    primary_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# CORRECT COLUMNS:
# patient_id = column 5
# cancer_type = column 4

primary_cancer <- unique(
    primary[, c(5,4)]
)

colnames(primary_cancer) <- c(
    "patient_id",
    "Cancer_Type"
)

# ============================================================
# CLOCK16 MUTATIONS
# ============================================================

mut <- read.delim(
    mut_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

colnames(mut) <- c(
    "patient_id",
    "sample",
    "gene",
    "effect",
    "amino_acid_change",
    "VAF",
    "Cancer_Type"
)

# Unique mutated patients
mut_patients <- unique(
    mut[, c("patient_id", "Cancer_Type")]
)

mutated_ids <- unique(mut_patients$patient_id)

# ============================================================
# CHECK
# ============================================================

cat("Primary patients:", nrow(primary_cancer), "\n")
cat("CLOCK16 mutated patients:", length(mutated_ids), "\n")
cat(
    "Matched mutated patients:",
    sum(primary_cancer$patient_id %in% mutated_ids),
    "\n"
)

# ============================================================
# CANCER-TYPE ENRICHMENT
# ============================================================

results <- data.frame()

for (ct in sort(unique(primary_cancer$Cancer_Type))) {

    sub <- primary_cancer[
        primary_cancer$Cancer_Type == ct,
    ]

    N <- nrow(sub)

    M <- sum(sub$patient_id %in% mutated_ids)

    # 2x2:
    #                    mutated   non-mutated
    # this cancer           M        N-M
    # other cancers        C-M       rest

    C <- length(mutated_ids)

    a <- M
    b <- N - M
    c <- C - M
    d <- nrow(primary_cancer) - N - c

    ft <- fisher.test(
        matrix(
            c(a,b,c,d),
            nrow = 2,
            byrow = TRUE
        )
    )

    results <- rbind(
        results,
        data.frame(
            Cancer_Type = ct,
            Primary_Patients = N,
            CLOCK16_Mutated_Patients = M,
            Prevalence_Percent = 100*M/N,
            Odds_Ratio = unname(ft$estimate),
            P_value = ft$p.value
        )
    )
}

# ============================================================
# FDR
# ============================================================

results$FDR <- p.adjust(
    results$P_value,
    method = "BH"
)

overall_prev <- length(mutated_ids) / nrow(primary_cancer)

results$Direction <- ifelse(
    results$Prevalence_Percent / 100 > overall_prev,
    "Enriched",
    "Depleted"
)

results <- results[
    order(results$FDR, results$P_value),
]

# ============================================================
# SAVE
# ============================================================

write.table(
    results,
    "results/clock16/clock16_cancer_enrichment.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

cat("\n=== CLOCK16 CANCER ENRICHMENT ===\n\n")
print(results, row.names = FALSE)

cat("\n=== FDR < 0.05 ===\n\n")

sig <- results[results$FDR < 0.05, ]

if (nrow(sig) == 0) {
    cat("No cancer type survives FDR < 0.05.\n")
} else {
    print(sig, row.names = FALSE)
}

