# compute_te_enrichment.R - Compute TE family/subfamily enrichment statistics
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Pure computation functions for TE enrichment analysis.
#              Uses limma's geneSetTest for competitive testing.
#              NO visualization - see plot_te_enrichment_dotplot.R for plotting.
#
# Dependencies: limma, dplyr, tibble, te_utils.R

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

#' Compute TE Family/Subfamily Enrichment Statistics
#'
#' Tests whether TEs grouped by subfamily, family, or class are enriched
#' among differentially expressed features. Uses limma's geneSetTest
#' for competitive gene set testing.
#'
#' @param fit A limma MArrayLM fit object (from eBayes)
#' @param dge A DGEList object with TE annotation in $genes slot.
#'            Must have columns: feature_type, subfamily, family, class
#' @param te_level Level of TE grouping: "subfamily", "family", or "class"
#' @param contrasts Character vector of contrast names to test. NULL = all.
#' @param p_adjust_method Method for p-value adjustment (default: "BH")
#' @param adjust_across How to adjust p-values: "all" (across all tests) or
#'                      "by_contrast" (within each contrast)
#' @param exclude_groups Character vector of TE groups to exclude
#'                       (e.g., c("Unspecified", "Unknown"))
#' @param alternative Test alternative: "either" (default), "up", "down", "mixed"
#' @param min_size Minimum number of TEs required in a group (default: 3)
#'
#' @return A list with components:
#'   \item{results}{Tibble with enrichment statistics (no visualization)}
#'   \item{metadata}{List with analysis parameters}
#'
#' @details
#' The function performs competitive gene set testing using t-statistics
#' from the limma fit. For each TE group and contrast:
#' - Extracts t-statistics for all TEs in the group
#' - Tests if these t-statistics are more extreme than expected
#' - Computes average logFC and direction
#'
#' This is a PURE COMPUTATION function. It returns data, not plots.
#' Use \code{\link{plot_te_enrichment_dotplot}} for visualization.
#'
#' @section P-value Adjustment:
#' When adjust_across = "all", p-values are adjusted across all
#' (TE group x contrast) combinations. This is more stringent but
#' controls the overall FDR.
#'
#' When adjust_across = "by_contrast", p-values are adjusted within
#' each contrast. This is less stringent but may be appropriate when
#' contrasts are independent hypotheses.
#'
#' @examples
#' \dontrun{
#' # Load data
#' dge <- readRDS("checkpoints/1.5_te_dge_combined.rds")
#' fit <- readRDS("checkpoints/1.5_te_fit.rds")
#'
#' # Compute enrichment at subfamily level
#' enrichment <- compute_te_enrichment(
#'     fit = fit,
#'     dge = dge,
#'     te_level = "subfamily",
#'     exclude_groups = c("Unspecified", "Unknown")
#' )
#'
#' # Access results
#' enrichment$results
#' enrichment$metadata
#'
#' # Test specific contrasts only
#' enrichment <- compute_te_enrichment(
#'     fit = fit,
#'     dge = dge,
#'     te_level = "family",
#'     contrasts = c("Treatment_vs_Control", "TimeB_vs_TimeA")
#' )
#' }
#'
#' @seealso \code{\link{format_te_enrichment}} for master table generation
#' @seealso \code{\link{plot_te_enrichment_dotplot}} for visualization
#'
#' @export
compute_te_enrichment <- function(fit,
                                  dge,
                                  te_level = c("subfamily", "family", "class"),
                                  contrasts = NULL,
                                  p_adjust_method = "BH",
                                  adjust_across = c("all", "by_contrast"),
                                  exclude_groups = NULL,
                                  alternative = c("either", "up", "down", "mixed"),
                                  min_size = 3) {

    # Load limma
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("limma package is required. Install with BiocManager::install('limma')")
    }

    # Match arguments
    te_level <- match.arg(te_level)
    adjust_across <- match.arg(adjust_across)
    alternative <- match.arg(alternative)

    message("=== Computing TE Enrichment (", te_level, " level) ===")

    # -------------------------------------------------------------------------
    # Validate inputs
    # -------------------------------------------------------------------------

    # Check fit object
    if (!inherits(fit, "MArrayLM")) {
        stop("fit must be a limma MArrayLM object (from eBayes)", call. = FALSE)
    }

    if (is.null(fit$t)) {
        stop("fit object must have t-statistics. Did you run eBayes()?", call. = FALSE)
    }

    # Check DGE has required annotation
    required_cols <- c("feature_type", te_level)
    missing_cols <- setdiff(required_cols, colnames(dge$genes))
    if (length(missing_cols) > 0) {
        stop("DGE$genes missing required columns: ",
             paste(missing_cols, collapse = ", "),
             "\nRun create_combined_dge() or create_separate_te_dge() first.",
             call. = FALSE)
    }

    # Check feature alignment
    if (!identical(rownames(fit), rownames(dge))) {
        stop("Rownames of fit and dge do not match. ",
             "Ensure fit was created from the same DGE.", call. = FALSE)
    }

    # -------------------------------------------------------------------------
    # Determine contrasts to test
    # -------------------------------------------------------------------------
    available_contrasts <- colnames(fit$coefficients)

    if (is.null(contrasts)) {
        contrasts <- available_contrasts
        message("[INFO] Testing all ", length(contrasts), " contrasts")
    } else {
        missing_contrasts <- setdiff(contrasts, available_contrasts)
        if (length(missing_contrasts) > 0) {
            stop("Contrasts not found in fit: ",
                 paste(missing_contrasts, collapse = ", "), call. = FALSE)
        }
        message("[INFO] Testing ", length(contrasts), " specified contrasts")
    }

    # -------------------------------------------------------------------------
    # Build TE group index
    # -------------------------------------------------------------------------
    message("[INDEX] Building TE group index at ", te_level, " level...")

    # Get annotation for TEs only
    te_mask <- dge$genes$feature_type == "te"
    te_annotation <- dge$genes[te_mask, ]

    # Get unique groups
    all_groups <- unique(te_annotation[[te_level]])
    all_groups <- all_groups[!is.na(all_groups)]

    # Exclude specified groups
    if (!is.null(exclude_groups)) {
        all_groups <- setdiff(all_groups, exclude_groups)
        message("[FILTER] Excluded ", length(exclude_groups), " groups: ",
                paste(exclude_groups, collapse = ", "))
    }

    message("[INDEX] Found ", length(all_groups), " TE groups to test")

    # Create index: for each group, which row indices in the fit?
    group_indices <- lapply(all_groups, function(grp) {
        which(dge$genes[[te_level]] == grp)
    })
    names(group_indices) <- all_groups

    # Filter by minimum size
    group_sizes <- sapply(group_indices, length)
    keep_groups <- names(group_sizes)[group_sizes >= min_size]

    if (length(keep_groups) < length(all_groups)) {
        n_filtered <- length(all_groups) - length(keep_groups)
        message("[FILTER] Removed ", n_filtered, " groups with < ", min_size, " TEs")
    }

    group_indices <- group_indices[keep_groups]
    message("[INDEX] Testing ", length(group_indices), " groups with >= ",
            min_size, " TEs each")

    # -------------------------------------------------------------------------
    # Compute enrichment for each group x contrast
    # -------------------------------------------------------------------------
    message("[COMPUTE] Running geneSetTest for ",
            length(group_indices), " groups x ", length(contrasts), " contrasts...")

    results_list <- list()
    total_tests <- length(group_indices) * length(contrasts)
    test_count <- 0

    for (contrast in contrasts) {
        t_stats <- fit$t[, contrast]
        logfc <- fit$coefficients[, contrast]

        for (grp in names(group_indices)) {
            test_count <- test_count + 1

            idx <- group_indices[[grp]]

            # Gene set test
            pval <- limma::geneSetTest(
                index = idx,
                statistics = t_stats,
                alternative = alternative,
                type = "t"
            )

            # Compute summary statistics for the group
            grp_t <- t_stats[idx]
            grp_logfc <- logfc[idx]

            results_list[[test_count]] <- tibble(
                te_group = grp,
                te_level = te_level,
                contrast = contrast,
                pval_raw = pval,
                t_mean = mean(grp_t, na.rm = TRUE),
                t_median = median(grp_t, na.rm = TRUE),
                avg_logfc = mean(grp_logfc, na.rm = TRUE),
                median_logfc = median(grp_logfc, na.rm = TRUE),
                n_elements = length(idx),
                n_up = sum(grp_logfc > 0 & fit$p.value[idx, contrast] < 0.05, na.rm = TRUE),
                n_down = sum(grp_logfc < 0 & fit$p.value[idx, contrast] < 0.05, na.rm = TRUE)
            )
        }
    }

    # Combine results
    results <- bind_rows(results_list)

    message("[COMPUTE] Completed ", nrow(results), " tests")

    # -------------------------------------------------------------------------
    # P-value adjustment
    # -------------------------------------------------------------------------
    message("[ADJUST] Adjusting p-values (", p_adjust_method, ", across=", adjust_across, ")...")

    if (adjust_across == "all") {
        results <- results %>%
            mutate(padj = p.adjust(pval_raw, method = p_adjust_method))
    } else {
        # Adjust within each contrast
        results <- results %>%
            group_by(contrast) %>%
            mutate(padj = p.adjust(pval_raw, method = p_adjust_method)) %>%
            ungroup()
    }

    # -------------------------------------------------------------------------
    # Add derived columns
    # -------------------------------------------------------------------------
    results <- results %>%
        mutate(
            direction = case_when(
                padj >= 0.05 ~ "ns",
                avg_logfc > 0 ~ "up",
                avg_logfc < 0 ~ "down",
                TRUE ~ "ns"
            ),
            prop_significant = (n_up + n_down) / n_elements
        )

    # Add TE class and replication type mapping
    class_map <- map_te_to_class(dge$genes, from_level = te_level)
    replication_map <- map_te_to_replication(dge$genes, from_level = te_level)

    results <- results %>%
        mutate(
            te_class = class_map[te_group],
            replication_type = replication_map[te_group]
        ) %>%
        relocate(te_class, replication_type, .after = te_level)

    # -------------------------------------------------------------------------
    # Build metadata
    # -------------------------------------------------------------------------
    metadata <- list(
        te_level = te_level,
        contrasts = contrasts,
        n_contrasts = length(contrasts),
        n_groups = length(group_indices),
        n_tests = nrow(results),
        p_adjust_method = p_adjust_method,
        adjust_across = adjust_across,
        alternative = alternative,
        min_size = min_size,
        excluded_groups = exclude_groups,
        analysis_mode = ifelse(is.null(dge$analysis_mode), "unknown", dge$analysis_mode),
        computed_at = Sys.time(),
        toolkit_version = "2.0.0"
    )

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    n_sig <- sum(results$padj < 0.05)
    message("=== TE Enrichment Complete ===")
    message("  Tests: ", nrow(results))
    message("  Significant (padj < 0.05): ", n_sig, " (",
            round(n_sig / nrow(results) * 100, 1), "%)")
    message("  Up-regulated groups: ", sum(results$direction == "up"))
    message("  Down-regulated groups: ", sum(results$direction == "down"))

    return(list(
        results = results,
        metadata = metadata
    ))
}


#' Compute TE Enrichment for Multiple Levels
#'
#' Convenience function to compute enrichment at multiple TE levels
#' (subfamily, family, class) in one call.
#'
#' @param fit A limma MArrayLM fit object
#' @param dge A DGEList object with TE annotation
#' @param te_levels Character vector of levels to test (default: all three)
#' @param ... Additional arguments passed to compute_te_enrichment
#'
#' @return Named list of enrichment results, one per level
#'
#' @examples
#' \dontrun{
#' all_enrichment <- compute_te_enrichment_multi(
#'     fit = fit,
#'     dge = dge,
#'     te_levels = c("subfamily", "family"),
#'     exclude_groups = c("Unspecified")
#' )
#'
#' # Access results
#' all_enrichment$subfamily$results
#' all_enrichment$family$results
#' }
#'
#' @export
compute_te_enrichment_multi <- function(fit,
                                        dge,
                                        te_levels = c("subfamily", "family", "class"),
                                        ...) {

    results_list <- list()

    for (level in te_levels) {
        message("\n")
        results_list[[level]] <- compute_te_enrichment(
            fit = fit,
            dge = dge,
            te_level = level,
            ...
        )
    }

    return(results_list)
}


#' Quick summary of enrichment results
#'
#' @param enrichment Output from compute_te_enrichment
#' @param alpha Significance threshold (default: 0.05)
#' @return Invisible NULL (prints summary)
#'
#' @export
summarize_te_enrichment <- function(enrichment, alpha = 0.05) {

    results <- enrichment$results
    metadata <- enrichment$metadata

    cat("\n=== TE Enrichment Summary ===\n")
    cat("Level:", metadata$te_level, "\n")
    cat("Contrasts:", metadata$n_contrasts, "\n")
    cat("Groups tested:", metadata$n_groups, "\n")
    cat("Total tests:", metadata$n_tests, "\n")
    cat("P-adjustment:", metadata$p_adjust_method, "(", metadata$adjust_across, ")\n")
    cat("Analysis mode:", metadata$analysis_mode, "\n\n")

    # Summary by contrast
    cat("Results by contrast (padj <", alpha, "):\n")
    summary_by_contrast <- results %>%
        group_by(contrast) %>%
        summarize(
            n_sig = sum(padj < alpha),
            n_up = sum(direction == "up"),
            n_down = sum(direction == "down"),
            .groups = "drop"
        )

    print(as.data.frame(summary_by_contrast), row.names = FALSE)

    # Top significant groups
    cat("\nTop 10 most significant groups:\n")
    top_sig <- results %>%
        arrange(padj) %>%
        head(10) %>%
        select(te_group, contrast, padj, avg_logfc, direction, n_elements)

    print(as.data.frame(top_sig), row.names = FALSE)

    invisible(NULL)
}
