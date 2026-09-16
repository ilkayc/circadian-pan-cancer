
library(data.table)

cat("\n")
cat("============================================================\n")
cat("13A — CLOCK16 CANCER-TYPE ASSOCIATION v2\n")
cat("============================================================\n\n")


# ============================================================
# 1. INPUT FILES
# ============================================================

primary_file <- "results/denominator/primary_samples_with_patient_id.tsv"

any_file <- "results/clock16/REAL_primary_CLOCK16_patient_gene_summary.tsv"

nonsyn_file <- "results/clock16/REAL_CLOCK16_nonsynonymous_functional_classification.tsv.gz"

if (!file.exists(primary_file))
    stop("Missing primary patient denominator: ", primary_file)

if (!file.exists(any_file))
    stop("Missing any-mutation patient file: ", any_file)

if (!file.exists(nonsyn_file))
    stop("Missing nonsynonymous file: ", nonsyn_file)


# ============================================================
# 2. READ DATA
# ============================================================

primary <- fread(primary_file)

any_mut <- fread(any_file)

nonsyn <- fread(nonsyn_file)


cat("Primary rows:", nrow(primary), "\n")
cat("Any-mutation rows:", nrow(any_mut), "\n")
cat("Nonsynonymous rows:", nrow(nonsyn), "\n\n")


# ============================================================
# 3. CHECK PRIMARY DENOMINATOR
# ============================================================

required_primary <- c(
    "patient_id",
    "cancer_type"
)

missing_primary <- setdiff(
    required_primary,
    names(primary)
)

if (length(missing_primary) > 0) {

    stop(
        "Missing primary columns: ",
        paste(missing_primary, collapse = ", ")
    )

}


# Keep only primary patient + cancer type

patient_cancer <- unique(
    primary[
        !is.na(patient_id) &
        !is.na(cancer_type),
        .(
            patient_id,
            cancer_type
        )
    ]
)


cat("Unique primary patients:", uniqueN(patient_cancer$patient_id), "\n")
cat("Cancer types:", uniqueN(patient_cancer$cancer_type), "\n")


# ============================================================
# 4. CRITICAL PATIENT-LEVEL QC
# ============================================================

patient_counts <- patient_cancer[
    ,
    .N,
    by = patient_id
]

if (any(patient_counts$N > 1)) {

    problematic <- patient_counts[N > 1]

    stop(
        "ERROR: Some patients have multiple cancer types. ",
        nrow(problematic),
        " problematic patients."
    )

}


primary_n <- uniqueN(patient_cancer$patient_id)

if (primary_n != 10593) {

    stop(
        "ERROR: Expected 10,593 primary patients but found ",
        primary_n
    )

}


cat("PASS: Primary denominator = 10,593 patients.\n\n")


# ============================================================
# 5. ANY CLOCK16 MUTATION
# ============================================================

if (!"patient_id" %in% names(any_mut)) {

    stop("patient_id missing from any-mutation file.")

}

any_patients <- unique(
    any_mut[
        !is.na(patient_id),
        patient_id
    ]
)


# Restrict to denominator

any_patients <- intersect(
    any_patients,
    patient_cancer$patient_id
)


cat("Any CLOCK16 patients:", length(any_patients), "\n")


if (length(any_patients) != 1165) {

    stop(
        "ERROR: Expected 1,165 any CLOCK16 patients but found ",
        length(any_patients)
    )

}


cat("PASS: Any CLOCK16 = 1,165 patients.\n\n")


# ============================================================
# 6. NONSYNONYMOUS CLOCK16
# ============================================================

if (!"patient_id" %in% names(nonsyn)) {

    stop("patient_id missing from nonsynonymous file.")

}


nonsyn_patients <- unique(
    nonsyn[
        !is.na(patient_id),
        patient_id
    ]
)


nonsyn_patients <- intersect(
    nonsyn_patients,
    patient_cancer$patient_id
)


cat(
    "Nonsynonymous CLOCK16 patients:",
    length(nonsyn_patients),
    "\n"
)


if (length(nonsyn_patients) != 854) {

    stop(
        "ERROR: Expected 854 nonsynonymous patients but found ",
        length(nonsyn_patients)
    )

}


cat("PASS: Nonsynonymous CLOCK16 = 854 patients.\n\n")


# ============================================================
# 7. DAMAGING CLOCK16
# ============================================================

if (!"Functional_Class" %in% names(nonsyn)) {

    stop("Functional_Class missing from nonsynonymous file.")

}


damaging_classes <- c(
    "strong_damaging",
    "damaging"
)


damaging_patients <- unique(
    nonsyn[
        Functional_Class %in% damaging_classes &
        !is.na(patient_id),
        patient_id
    ]
)


damaging_patients <- intersect(
    damaging_patients,
    patient_cancer$patient_id
)


cat(
    "Damaging CLOCK16 patients:",
    length(damaging_patients),
    "\n"
)


if (length(damaging_patients) != 395) {

    stop(
        "ERROR: Expected 395 damaging patients but found ",
        length(damaging_patients)
    )

}


cat("PASS: Damaging CLOCK16 = 395 patients.\n\n")


# ============================================================
# 8. MUTATION STATUS TABLE
# ============================================================

patient_cancer[
    ,
    any_CLOCK16 := as.integer(
        patient_id %in% any_patients
    )
]

patient_cancer[
    ,
    nonsynonymous_CLOCK16 := as.integer(
        patient_id %in% nonsyn_patients
    )
]

patient_cancer[
    ,
    damaging_CLOCK16 := as.integer(
        patient_id %in% damaging_patients
    )
]


# ============================================================
# 9. GLOBAL STATUS QC
# ============================================================

cat("============================================================\n")
cat("GLOBAL MUTATION STATUS QC\n")
cat("============================================================\n\n")

cat(
    "Any CLOCK16:",
    sum(patient_cancer$any_CLOCK16),
    "\n"
)

cat(
    "Nonsynonymous CLOCK16:",
    sum(patient_cancer$nonsynonymous_CLOCK16),
    "\n"
)

cat(
    "Damaging CLOCK16:",
    sum(patient_cancer$damaging_CLOCK16),
    "\n\n"
)


if (sum(patient_cancer$any_CLOCK16) != 1165)
    stop("Global any CLOCK16 count mismatch.")

if (sum(patient_cancer$nonsynonymous_CLOCK16) != 854)
    stop("Global nonsynonymous count mismatch.")

if (sum(patient_cancer$damaging_CLOCK16) != 395)
    stop("Global damaging count mismatch.")


cat("PASS: Global patient-level counts match expected values.\n\n")


# ============================================================
# 10. ASSOCIATION FUNCTION
# ============================================================

run_association <- function(
    dat,
    status_col,
    phenotype_name
) {

    cancer_types <- sort(
        unique(dat$cancer_type)
    )

    results <- vector(
        "list",
        length(cancer_types)
    )


    for (i in seq_along(cancer_types)) {

        ct <- cancer_types[i]

        in_cancer <- dat$cancer_type == ct

        mutated <- dat[[status_col]] == 1

        a <- sum(in_cancer & mutated)
        b <- sum(in_cancer & !mutated)
        c <- sum(!in_cancer & mutated)
        d <- sum(!in_cancer & !mutated)

        total <- a + b

        prevalence <- 100 * a / total


        tab <- matrix(
            c(a, b, c, d),
            nrow = 2,
            byrow = TRUE
        )


        fisher <- fisher.test(tab)


        results[[i]] <- data.table(

            phenotype = phenotype_name,

            cancer_type = ct,

            total_primary_patients = total,

            mutated_patients = a,

            nonmutated_patients = b,

            prevalence_percent = prevalence,

            odds_ratio = unname(fisher$estimate),

            lower_95CI = fisher$conf.int[1],

            upper_95CI = fisher$conf.int[2],

            fisher_p = fisher$p.value

        )

    }


    result <- rbindlist(results)


    result[
        ,
        FDR := p.adjust(
            fisher_p,
            method = "BH"
        )
    ]


    result[
        order(FDR, fisher_p)
    ]

}


# ============================================================
# 11. RUN ASSOCIATIONS
# ============================================================

cat("============================================================\n")
cat("RUNNING ASSOCIATION TESTS\n")
cat("============================================================\n\n")


res_any <- run_association(
    patient_cancer,
    "any_CLOCK16",
    "Any_CLOCK16"
)

cat("Any CLOCK16: completed.\n")


res_nonsyn <- run_association(
    patient_cancer,
    "nonsynonymous_CLOCK16",
    "Nonsynonymous_CLOCK16"
)

cat("Nonsynonymous CLOCK16: completed.\n")


res_damaging <- run_association(
    patient_cancer,
    "damaging_CLOCK16",
    "Damaging_CLOCK16"
)

cat("Damaging CLOCK16: completed.\n\n")


# ============================================================
# 12. COMBINE
# ============================================================

all_results <- rbind(
    res_any,
    res_nonsyn,
    res_damaging
)


# ============================================================
# 13. SAVE
# ============================================================

fwrite(
    res_any,
    "results/clock16/REAL_CLOCK16_cancer_association_any_v2.tsv",
    sep = "\t"
)

fwrite(
    res_nonsyn,
    "results/clock16/REAL_CLOCK16_cancer_association_nonsynonymous_v2.tsv",
    sep = "\t"
)

fwrite(
    res_damaging,
    "results/clock16/REAL_CLOCK16_cancer_association_damaging_v2.tsv",
    sep = "\t"
)

fwrite(
    all_results,
    "results/clock16/REAL_CLOCK16_cancer_association_ALL_v2.tsv",
    sep = "\t"
)


# ============================================================
# 14. SIGNIFICANT RESULTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("13A RESULTS — FDR < 0.05\n")
cat("============================================================\n")


for (ph in unique(all_results$phenotype)) {

    cat("\n")
    cat("----------------------------------------\n")
    cat(ph, "\n")
    cat("----------------------------------------\n")

    sig <- all_results[
        phenotype == ph &
        FDR < 0.05
    ]

    if (nrow(sig) == 0) {

        cat("No cancer type significant after BH-FDR < 0.05.\n")

    } else {

        print(
            sig[
                order(FDR, fisher_p),
                .(
                    cancer_type,
                    total_primary_patients,
                    mutated_patients,
                    prevalence_percent,
                    odds_ratio,
                    lower_95CI,
                    upper_95CI,
                    fisher_p,
                    FDR
                )
            ]
        )

    }

}


# ============================================================
# 15. TOP 15 RESULTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TOP 15 RESULTS BY FDR\n")
cat("============================================================\n\n")


print(
    all_results[
        order(FDR, fisher_p)
    ][
        1:min(15, .N),
        .(
            phenotype,
            cancer_type,
            total_primary_patients,
            mutated_patients,
            prevalence_percent,
            odds_ratio,
            lower_95CI,
            upper_95CI,
            fisher_p,
            FDR
        )
    ]
)


# ============================================================
# 16. FINAL QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("13A FINAL QC\n")
cat("============================================================\n\n")


cat(
    "Primary patients:",
    uniqueN(patient_cancer$patient_id),
    "\n"
)

cat(
    "Cancer types:",
    uniqueN(patient_cancer$cancer_type),
    "\n"
)

cat(
    "Any CLOCK16:",
    sum(patient_cancer$any_CLOCK16),
    "\n"
)

cat(
    "Nonsynonymous CLOCK16:",
    sum(patient_cancer$nonsynonymous_CLOCK16),
    "\n"
)

cat(
    "Damaging CLOCK16:",
    sum(patient_cancer$damaging_CLOCK16),
    "\n"
)


if (
    uniqueN(patient_cancer$patient_id) == 10593 &&
    sum(patient_cancer$any_CLOCK16) == 1165 &&
    sum(patient_cancer$nonsynonymous_CLOCK16) == 854 &&
    sum(patient_cancer$damaging_CLOCK16) == 395
) {

    cat("\n")
    cat("13A STATUS: PASS\n")
    cat("\n")
    cat("Correct 10,593-patient denominator used.\n")
    cat("Mutation status assigned at patient level.\n")
    cat("All expected mutation counts reproduced.\n")

} else {

    stop(
        "\n13A STATUS: FAIL — QC mismatch."
    )

}


cat("\n")
cat("============================================================\n")
cat("13A COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

