# validate_te_input.R - Input validation functions for TE analysis
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Validation functions to ensure data integrity before analysis.
#              These functions fail fast with informative error messages.
#
# Dependencies: None (base R only)

# =============================================================================
# ID VALIDATION
# =============================================================================

#' Validate No ID Collision Between Gene and TE Matrices
#'
#' Checks that no feature IDs appear in both the gene and TE count matrices.
#' ID collisions indicate annotation overlap, which would cause double-counting.
#'
#' @param gene_counts Gene count matrix (features x samples)
#' @param te_counts TE count matrix (features x samples)
#' @param stop_on_error If TRUE (default), throws error on collision.
#'                      If FALSE, returns FALSE instead.
#' @return TRUE if valid (invisible), or FALSE if stop_on_error=FALSE and collision found
#'
#' @details
#' ID collision can occur when:
#' - Exonic TEs were not removed from the TE annotation
#' - Gene IDs use a format similar to TE labels
#' - The same feature was quantified in both matrices
#'
#' @examples
#' \dontrun{
#' validate_no_id_collision(gene_counts, te_counts)
#' }
#'
#' @export
validate_no_id_collision <- function(gene_counts, te_counts, stop_on_error = TRUE) {
    gene_ids <- rownames(gene_counts)
    te_ids <- rownames(te_counts)

    if (is.null(gene_ids) || is.null(te_ids)) {
        if (stop_on_error) {
            stop("Count matrices must have rownames (feature IDs)")
        }
        return(FALSE)
    }

    collision <- intersect(gene_ids, te_ids)

    if (length(collision) > 0) {
        msg <- paste0(
            "ID collision detected between gene and TE matrices!\n",
            "Found ", length(collision), " duplicate IDs:\n  ",
            paste(head(collision, 10), collapse = "\n  "),
            if (length(collision) > 10) paste0("\n  ... and ", length(collision) - 10, " more") else "",
            "\n\nThis usually means TE annotation overlaps with gene annotation.\n",
            "Solutions:\n",
            "  1. Use separate analysis mode: create_separate_te_dge()\n",
            "  2. Filter overlapping TEs from your SAF annotation\n",
            "  3. Review preprocessing pipeline"
        )

        if (stop_on_error) {
            stop(msg, call. = FALSE)
        } else {
            warning(msg, call. = FALSE)
            return(FALSE)
        }
    }

    message("[OK] No ID collision between gene and TE matrices (",
            length(gene_ids), " genes, ", length(te_ids), " TEs)")
    invisible(TRUE)
}


# =============================================================================
# SAMPLE VALIDATION
# =============================================================================

#' Validate Sample Consistency Between Matrices
#'
#' Checks that gene and TE count matrices have identical samples in the same order.
#'
#' @param gene_counts Gene count matrix (features x samples)
#' @param te_counts TE count matrix (features x samples)
#' @param stop_on_error If TRUE (default), throws error on mismatch.
#' @return TRUE if valid (invisible), or FALSE if stop_on_error=FALSE and mismatch found
#'
#' @examples
#' \dontrun{
#' validate_sample_consistency(gene_counts, te_counts)
#' }
#'
#' @export
validate_sample_consistency <- function(gene_counts, te_counts, stop_on_error = TRUE) {
    gene_samples <- colnames(gene_counts)
    te_samples <- colnames(te_counts)

    if (is.null(gene_samples) || is.null(te_samples)) {
        if (stop_on_error) {
            stop("Count matrices must have colnames (sample IDs)", call. = FALSE)
        }
        return(FALSE)
    }

    # Check for exact match (including order)
    if (!identical(gene_samples, te_samples)) {
        missing_in_te <- setdiff(gene_samples, te_samples)
        missing_in_gene <- setdiff(te_samples, gene_samples)
        order_mismatch <- length(missing_in_te) == 0 && length(missing_in_gene) == 0

        msg_parts <- c("Sample mismatch between gene and TE matrices!")

        if (length(missing_in_te) > 0) {
            msg_parts <- c(msg_parts,
                          paste0("Missing in TE matrix: ", paste(missing_in_te, collapse = ", ")))
        }

        if (length(missing_in_gene) > 0) {
            msg_parts <- c(msg_parts,
                          paste0("Missing in Gene matrix: ", paste(missing_in_gene, collapse = ", ")))
        }

        if (order_mismatch) {
            msg_parts <- c(msg_parts,
                          "Samples present in both but in different order.",
                          "Reorder columns to match before combining.")
        }

        msg <- paste(msg_parts, collapse = "\n")

        if (stop_on_error) {
            stop(msg, call. = FALSE)
        } else {
            warning(msg, call. = FALSE)
            return(FALSE)
        }
    }

    message("[OK] Sample consistency validated (", length(gene_samples), " samples)")
    invisible(TRUE)
}


#' Validate Sample Metadata Matches Count Matrix
#'
#' @param counts Count matrix (features x samples)
#' @param samples Sample metadata data.frame with rownames = sample IDs
#' @param stop_on_error If TRUE (default), throws error on mismatch
#' @return TRUE if valid (invisible)
#'
#' @export
validate_samples_metadata <- function(counts, samples, stop_on_error = TRUE) {
    count_samples <- colnames(counts)
    meta_samples <- rownames(samples)

    if (is.null(meta_samples)) {
        # Try to use a common ID column if rownames not set
        if ("sample_id" %in% colnames(samples)) {
            meta_samples <- samples$sample_id
        } else if ("ID" %in% colnames(samples)) {
            meta_samples <- samples$ID
        } else {
            if (stop_on_error) {
                stop("Sample metadata must have rownames or a 'sample_id'/'ID' column",
                     call. = FALSE)
            }
            return(FALSE)
        }
    }

    missing_in_meta <- setdiff(count_samples, meta_samples)
    missing_in_counts <- setdiff(meta_samples, count_samples)

    if (length(missing_in_meta) > 0 || length(missing_in_counts) > 0) {
        msg_parts <- c("Sample metadata does not match count matrix!")

        if (length(missing_in_meta) > 0) {
            msg_parts <- c(msg_parts,
                          paste0("Missing in metadata: ",
                                paste(head(missing_in_meta, 5), collapse = ", "),
                                if (length(missing_in_meta) > 5)
                                    paste0(" ... (", length(missing_in_meta), " total)") else ""))
        }

        if (length(missing_in_counts) > 0) {
            msg_parts <- c(msg_parts,
                          paste0("Missing in counts: ",
                                paste(head(missing_in_counts, 5), collapse = ", "),
                                if (length(missing_in_counts) > 5)
                                    paste0(" ... (", length(missing_in_counts), " total)") else ""))
        }

        msg <- paste(msg_parts, collapse = "\n")

        if (stop_on_error) {
            stop(msg, call. = FALSE)
        } else {
            warning(msg, call. = FALSE)
            return(FALSE)
        }
    }

    message("[OK] Sample metadata matches count matrix (", length(count_samples), " samples)")
    invisible(TRUE)
}


# =============================================================================
# LIBRARY SIZE VALIDATION
# =============================================================================

#' Validate Library Size Ratio (Warning Only)
#'
#' Checks if TE proportion is within expected range. High TE proportions
#' may indicate preprocessing issues or extreme TE derepression.
#'
#' @param gene_counts Gene count matrix (features x samples)
#' @param te_counts TE count matrix (features x samples)
#' @param warn_threshold Threshold for warning (default 0.3 = 30%)
#' @param error_threshold Threshold for error (default 0.5 = 50%). Set to 1 to disable.
#' @return Invisible vector of TE proportions per sample
#'
#' @details
#' Normal TE proportions are typically 2-10% of total mapped reads.
#' Higher proportions may indicate:
#' - Exonic TEs not removed (double-counting)
#' - Extreme TE derepression (biologically interesting)
#' - Preprocessing/annotation issues
#'
#' @examples
#' \dontrun{
#' te_props <- validate_library_size_ratio(gene_counts, te_counts)
#' }
#'
#' @export
validate_library_size_ratio <- function(gene_counts, te_counts,
                                        warn_threshold = 0.3,
                                        error_threshold = 0.5) {
    gene_lib_sizes <- colSums(gene_counts)
    te_lib_sizes <- colSums(te_counts)
    total_lib_sizes <- gene_lib_sizes + te_lib_sizes
    te_proportion <- te_lib_sizes / total_lib_sizes

    # Check for errors (likely preprocessing issues)
    if (any(te_proportion > error_threshold)) {
        problem_samples <- names(te_proportion)[te_proportion > error_threshold]
        stop(
            "Extremely high TE proportion detected!\n",
            "Samples above ", error_threshold * 100, "% threshold:\n  ",
            paste(problem_samples, collapse = ", "), "\n",
            "Max TE proportion: ", round(max(te_proportion) * 100, 1), "%\n\n",
            "This likely indicates a preprocessing issue:\n",
            "  - Check if exonic TEs were removed from annotation\n",
            "  - Verify gene and TE counts are not double-counted\n",
            "  - Review featureCounts parameters",
            call. = FALSE
        )
    }

    # Check for warnings (unusual but potentially valid)
    if (any(te_proportion > warn_threshold)) {
        problem_samples <- names(te_proportion)[te_proportion > warn_threshold]
        warning(
            "[CAUTION] High TE proportion detected in some samples:\n",
            "Samples above ", warn_threshold * 100, "% threshold: ",
            paste(problem_samples, collapse = ", "), "\n",
            "Max TE proportion: ", round(max(te_proportion) * 100, 1), "%\n",
            "This may indicate preprocessing issues or extreme TE derepression.\n",
            "Review your annotation and counting pipeline if unexpected.",
            call. = FALSE
        )
    } else {
        message("[OK] TE proportions within normal range ",
                "(range: ", round(min(te_proportion) * 100, 1), "% - ",
                round(max(te_proportion) * 100, 1), "%)")
    }

    invisible(te_proportion)
}


# =============================================================================
# COUNT MATRIX VALIDATION
# =============================================================================

#' Validate Count Matrix Structure
#'
#' Checks that a count matrix has valid structure for analysis.
#'
#' @param counts Count matrix
#' @param name Name of the matrix (for error messages)
#' @param min_features Minimum number of features required
#' @param min_samples Minimum number of samples required
#' @param stop_on_error If TRUE (default), throws error on invalid structure
#' @return TRUE if valid (invisible)
#'
#' @export
validate_count_matrix <- function(counts, name = "Count matrix",
                                  min_features = 1, min_samples = 2,
                                  stop_on_error = TRUE) {
    issues <- character()

    # Check it's a matrix or data.frame
    if (!is.matrix(counts) && !is.data.frame(counts)) {
        issues <- c(issues, paste0(name, " must be a matrix or data.frame"))
    }

    # Check dimensions
    if (nrow(counts) < min_features) {
        issues <- c(issues, paste0(name, " has ", nrow(counts),
                                  " features (minimum: ", min_features, ")"))
    }

    if (ncol(counts) < min_samples) {
        issues <- c(issues, paste0(name, " has ", ncol(counts),
                                  " samples (minimum: ", min_samples, ")"))
    }

    # Check for names
    if (is.null(rownames(counts))) {
        issues <- c(issues, paste0(name, " must have rownames (feature IDs)"))
    }

    if (is.null(colnames(counts))) {
        issues <- c(issues, paste0(name, " must have colnames (sample IDs)"))
    }

    # Check for numeric data
    if (!all(sapply(counts, is.numeric))) {
        issues <- c(issues, paste0(name, " must contain numeric data"))
    }

    # Check for negative values
    if (any(counts < 0, na.rm = TRUE)) {
        issues <- c(issues, paste0(name, " contains negative values (not valid for counts)"))
    }

    # Check for NAs
    if (any(is.na(counts))) {
        n_na <- sum(is.na(counts))
        issues <- c(issues, paste0(name, " contains ", n_na, " NA values"))
    }

    if (length(issues) > 0) {
        msg <- paste(c(paste0(name, " validation failed:"), issues), collapse = "\n  - ")
        if (stop_on_error) {
            stop(msg, call. = FALSE)
        } else {
            warning(msg, call. = FALSE)
            return(FALSE)
        }
    }

    message("[OK] ", name, " structure validated (",
            nrow(counts), " features x ", ncol(counts), " samples)")
    invisible(TRUE)
}


# =============================================================================
# COMPREHENSIVE VALIDATION
# =============================================================================

#' Run All Validations for Combined Analysis
#'
#' Convenience function that runs all relevant validations for combined
#' gene+TE analysis.
#'
#' @param gene_counts Gene count matrix
#' @param te_counts TE count matrix
#' @param samples Sample metadata data.frame (optional)
#' @param stop_on_error If TRUE (default), stops on first error
#' @return List with validation results
#'
#' @examples
#' \dontrun{
#' validation <- validate_combined_input(gene_counts, te_counts, samples)
#' }
#'
#' @export
validate_combined_input <- function(gene_counts, te_counts, samples = NULL,
                                    stop_on_error = TRUE) {
    message("=== Validating input for combined analysis ===")

    results <- list()

    # Validate individual matrices
    results$gene_matrix <- validate_count_matrix(gene_counts, "Gene counts",
                                                 stop_on_error = stop_on_error)
    results$te_matrix <- validate_count_matrix(te_counts, "TE counts",
                                               stop_on_error = stop_on_error)

    # Validate compatibility
    results$no_collision <- validate_no_id_collision(gene_counts, te_counts,
                                                     stop_on_error = stop_on_error)
    results$sample_consistency <- validate_sample_consistency(gene_counts, te_counts,
                                                              stop_on_error = stop_on_error)

    # Validate library sizes
    results$te_proportion <- validate_library_size_ratio(gene_counts, te_counts)

    # Validate metadata if provided
    if (!is.null(samples)) {
        results$metadata <- validate_samples_metadata(gene_counts, samples,
                                                      stop_on_error = stop_on_error)
    }

    message("=== All validations passed ===")
    invisible(results)
}


#' Run All Validations for Separate Analysis
#'
#' @param te_counts TE count matrix
#' @param samples Sample metadata data.frame (optional)
#' @param stop_on_error If TRUE (default), stops on first error
#' @return List with validation results
#'
#' @export
validate_separate_input <- function(te_counts, samples = NULL,
                                    stop_on_error = TRUE) {
    message("=== Validating input for separate TE analysis ===")

    results <- list()

    # Validate TE matrix
    results$te_matrix <- validate_count_matrix(te_counts, "TE counts",
                                               stop_on_error = stop_on_error)

    # Validate metadata if provided
    if (!is.null(samples)) {
        results$metadata <- validate_samples_metadata(te_counts, samples,
                                                      stop_on_error = stop_on_error)
    }

    message("=== All validations passed ===")
    invisible(results)
}
