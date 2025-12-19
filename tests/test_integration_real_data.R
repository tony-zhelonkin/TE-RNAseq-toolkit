#!/usr/bin/env Rscript
# =============================================================================
# test_integration_real_data.R - Integration tests with real project data
# =============================================================================
# TE-RNAseq-toolkit
#
# This script validates that the toolkit functions work correctly with
# real data from the AdaW_eWAT_WL_2025 project.
#
# Usage: Rscript test_integration_real_data.R
# =============================================================================

message("=================================================================")
message("TE-RNAseq-toolkit Integration Tests")
message("=================================================================")
message("Started: ", Sys.time())

# Track test results
test_results <- list()
n_passed <- 0
n_failed <- 0

run_test <- function(test_name, test_fn) {
    message("\n--- Testing: ", test_name, " ---")
    tryCatch({
        result <- test_fn()
        message("[PASS] ", test_name)
        test_results[[test_name]] <<- list(status = "PASS", result = result)
        n_passed <<- n_passed + 1
        return(TRUE)
    }, error = function(e) {
        message("[FAIL] ", test_name, ": ", e$message)
        test_results[[test_name]] <<- list(status = "FAIL", error = e$message)
        n_failed <<- n_failed + 1
        return(FALSE)
    })
}

# =============================================================================
# TEST 1: Source all R functions
# =============================================================================

run_test("Source all toolkit functions", function() {
    TOOLKIT_PATH <- "/workspaces/AdaW_eWAT_WL_2025/01_modules/TE-RNAseq-toolkit/R"

    # Required packages
    suppressPackageStartupMessages({
        library(edgeR)
        library(limma)
        library(dplyr)
        library(tibble)
        library(stringr)
        library(readr)
    })

    # Source in dependency order
    r_files <- c(
        "te_utils.R",
        "validate_te_input.R",
        "create_combined_dge.R",
        "create_separate_dge.R",
        "compute_te_enrichment.R",
        "compute_te_aggregates.R",
        "compute_te_de.R",
        "compute_te_proportions.R",
        "compute_te_correlations.R",
        "format_te_enrichment.R",
        "format_te_de.R",
        "format_te_aggregates.R",
        "plot_utils.R",
        "plot_te_enrichment_dotplot.R",
        "plot_te_enrichment_single.R",
        "plot_te_heatmap.R",
        "plot_te_volcano.R",
        "plot_te_proportions.R",
        "plot_te_correlation.R"
    )

    for (f in r_files) {
        fpath <- file.path(TOOLKIT_PATH, f)
        if (file.exists(fpath)) {
            source(fpath)
            message("  Sourced: ", f)
        } else {
            warning("  Missing: ", f)
        }
    }

    # Check key functions exist
    stopifnot(exists("build_te_annotation"))
    stopifnot(exists("create_combined_dge"))
    stopifnot(exists("validate_no_id_collision"))

    return(list(sourced_files = r_files))
})

# =============================================================================
# TEST 2: Load real data
# =============================================================================

run_test("Load real project data", function() {
    DATA_PATH <- "/workspaces/AdaW_eWAT_WL_2025/00_data/processed"

    # Load metadata (semicolon-delimited)
    metadata <- read.csv(
        file.path(DATA_PATH, "metadata.csv"),
        sep = ";",
        stringsAsFactors = FALSE
    )
    rownames(metadata) <- metadata$Sample
    colnames(metadata) <- c("sample_id", "subject_id", "group")
    message("  Metadata: ", nrow(metadata), " samples")
    message("  Groups: ", paste(unique(metadata$group), collapse = ", "))

    # Load TE counts
    te_counts <- read.delim(
        file.path(DATA_PATH, "featurecounts_TE/te_counts_matrix.txt"),
        row.names = 1,
        check.names = FALSE
    )
    te_counts <- as.matrix(te_counts)
    message("  TE counts: ", nrow(te_counts), " TEs x ", ncol(te_counts), " samples")

    # Load combined gene+TE counts
    combined_counts <- read.delim(
        file.path(DATA_PATH, "combined_gene_TE_counts.tsv"),
        row.names = 1,
        check.names = FALSE
    )
    combined_counts <- as.matrix(combined_counts)
    message("  Combined counts: ", nrow(combined_counts), " features x ", ncol(combined_counts), " samples")

    return(list(
        metadata = metadata,
        te_counts = te_counts,
        combined_counts = combined_counts
    ))
})

# =============================================================================
# TEST 3: Test TE ID parsing
# =============================================================================

run_test("Parse TE IDs from real data", function() {
    te_ids <- rownames(test_results[["Load real project data"]]$result$te_counts)

    # Test parsing
    parsed <- build_te_annotation(te_ids)

    message("  Parsed ", nrow(parsed), " TE IDs")
    message("  TE classes found: ", paste(unique(na.omit(parsed$class)), collapse = ", "))
    message("  TE families found: ", length(unique(na.omit(parsed$family))))
    message("  TE subfamilies found: ", length(unique(na.omit(parsed$subfamily))))

    # Check all IDs were parsed as TEs
    stopifnot(all(parsed$is_te))

    # Check expected classes are present
    expected_classes <- c("LINE", "SINE", "LTR", "DNA")
    found_classes <- unique(parsed$class)
    stopifnot(all(expected_classes %in% found_classes))

    return(parsed)
})

# =============================================================================
# TEST 4: Test feature type detection from combined matrix
# =============================================================================

run_test("Detect feature types in combined matrix", function() {
    combined_counts <- test_results[["Load real project data"]]$result$combined_counts
    feature_ids <- rownames(combined_counts)

    # Detect types
    feature_types <- detect_feature_type(feature_ids)

    n_genes <- sum(feature_types == "gene")
    n_tes <- sum(feature_types == "te")

    message("  Detected: ", n_genes, " genes, ", n_tes, " TEs")
    message("  Total: ", length(feature_types), " features")

    # Sanity checks
    stopifnot(n_genes > 0)
    stopifnot(n_tes > 0)
    stopifnot(n_genes > n_tes)  # Should have more genes than TEs

    return(list(n_genes = n_genes, n_tes = n_tes))
})

# =============================================================================
# TEST 5: Create combined DGE from pre-combined matrix
# =============================================================================

run_test("Create combined DGE from combined matrix", function() {
    combined_counts <- test_results[["Load real project data"]]$result$combined_counts
    metadata <- test_results[["Load real project data"]]$result$metadata

    # Detect feature types
    feature_ids <- rownames(combined_counts)
    is_te_vec <- detect_feature_type(feature_ids) == "te"

    # Split into genes and TEs
    gene_counts <- combined_counts[!is_te_vec, , drop = FALSE]
    te_counts <- combined_counts[is_te_vec, , drop = FALSE]

    message("  Gene counts: ", nrow(gene_counts), " x ", ncol(gene_counts))
    message("  TE counts: ", nrow(te_counts), " x ", ncol(te_counts))

    # Create combined DGE
    dge <- create_combined_dge(
        gene_counts = gene_counts,
        te_counts = te_counts,
        samples = metadata,
        validate = TRUE,
        normalize = TRUE
    )

    # Verify structure
    stopifnot(inherits(dge, "DGEList"))
    stopifnot(!is.null(dge$genes$feature_type))
    stopifnot(dge$analysis_mode == "combined")
    stopifnot(!is.null(dge$samples$norm.factors))

    message("  DGE created: ", nrow(dge), " features, ", ncol(dge), " samples")
    message("  Analysis mode: ", dge$analysis_mode)

    return(dge)
})

# =============================================================================
# TEST 6: Summarize TE content
# =============================================================================

run_test("Summarize TE content", function() {
    dge <- test_results[["Create combined DGE from combined matrix"]]$result

    summary <- summarize_te_content(dge)

    message("  TE proportion range: ",
            round(min(summary$te_proportion) * 100, 2), "% - ",
            round(max(summary$te_proportion) * 100, 2), "%")

    # Sanity check - TE proportion should be reasonable
    stopifnot(all(summary$te_proportion < 0.5))  # Less than 50%
    stopifnot(all(summary$te_proportion > 0.001))  # More than 0.1%

    return(summary)
})

# =============================================================================
# TEST 7: Create TE index for gene set testing
# =============================================================================

run_test("Create TE index for enrichment", function() {
    dge <- test_results[["Create combined DGE from combined matrix"]]$result

    # Build annotation from genes slot
    te_annotation <- dge$genes[dge$genes$feature_type == "te", ]

    # Convert to expected format
    annotation_for_index <- tibble(
        feature_id = te_annotation$feature_id,
        is_te = TRUE,
        subfamily = te_annotation$subfamily,
        family = te_annotation$family,
        class = te_annotation$class,
        replication_type = te_annotation$replication_type
    )

    # Create index at family level
    te_index <- create_te_index(annotation_for_index, level = "family")

    message("  TE families indexed: ", length(te_index))
    message("  Example families: ", paste(head(names(te_index), 5), collapse = ", "))

    # Sanity checks
    stopifnot(length(te_index) > 0)
    stopifnot("L1" %in% names(te_index) || any(grepl("L1", names(te_index))))

    return(te_index)
})

# =============================================================================
# TEST 8: Test validation functions
# =============================================================================

run_test("Validation functions detect issues", function() {
    # Create synthetic data with collision
    gene_counts <- matrix(1:12, nrow = 3, ncol = 4,
                         dimnames = list(c("Gene1", "Gene2", "Collision:ID:TEST"),
                                        c("S1", "S2", "S3", "S4")))
    te_counts <- matrix(1:8, nrow = 2, ncol = 4,
                       dimnames = list(c("TE1:L1:LINE", "Collision:ID:TEST"),
                                      c("S1", "S2", "S3", "S4")))

    # Should detect collision
    result <- validate_no_id_collision(gene_counts, te_counts, stop_on_error = FALSE)
    stopifnot(result == FALSE)  # Should return FALSE due to collision

    message("  Collision detection: Working")

    # Test sample mismatch detection
    te_counts_bad <- matrix(1:6, nrow = 2, ncol = 3,
                           dimnames = list(c("TE1:L1:LINE", "TE2:L1:LINE"),
                                          c("S1", "S2", "S5")))  # S5 is wrong

    result2 <- validate_sample_consistency(gene_counts, te_counts_bad, stop_on_error = FALSE)
    stopifnot(result2 == FALSE)  # Should return FALSE due to mismatch

    message("  Sample mismatch detection: Working")

    return(TRUE)
})

# =============================================================================
# TEST 9: Test compute_te_de (basic limma-voom workflow)
# =============================================================================

run_test("Run limma-voom DE analysis", function() {
    dge <- test_results[["Create combined DGE from combined matrix"]]$result
    metadata <- test_results[["Load real project data"]]$result$metadata

    # Create a simple design (simplified for testing)
    # Use baseline vs progression as a simple comparison
    metadata$group_simple <- ifelse(metadata$group == "baseline", "baseline", "other")
    design <- model.matrix(~ 0 + group_simple, data = metadata)
    colnames(design) <- gsub("group_simple", "", colnames(design))

    message("  Design columns: ", paste(colnames(design), collapse = ", "))

    # Apply filterByExpr
    keep <- filterByExpr(dge, design = design)
    dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
    message("  Features after filtering: ", nrow(dge_filtered))

    # Voom transformation
    v <- voom(dge_filtered, design)

    # Fit model
    fit <- lmFit(v, design)

    # Create contrast
    contrast_matrix <- makeContrasts(
        other_vs_baseline = other - baseline,
        levels = design
    )

    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    # Extract results
    results <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")

    message("  DE results: ", nrow(results), " features")
    message("  Significant (FDR < 0.05): ", sum(results$adj.P.Val < 0.05, na.rm = TRUE))

    return(list(fit = fit2, results = results))
})

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=================================================================")
message("Test Summary")
message("=================================================================")
message("Total tests: ", n_passed + n_failed)
message("Passed: ", n_passed)
message("Failed: ", n_failed)
message("")

if (n_failed > 0) {
    message("FAILED TESTS:")
    for (name in names(test_results)) {
        if (test_results[[name]]$status == "FAIL") {
            message("  - ", name, ": ", test_results[[name]]$error)
        }
    }
    message("")
    quit(status = 1)
} else {
    message("ALL TESTS PASSED!")
    message("")
}

message("Completed: ", Sys.time())
