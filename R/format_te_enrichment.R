# format_te_enrichment.R - Format TE enrichment results to master table schema
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Formats output from compute_te_enrichment() into a standardized
#              master table structure for export and cross-language compatibility.
#
# Dependencies: dplyr, tibble

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

#' Format TE Enrichment Results to Master Table Schema
#'
#' Transforms the raw output of compute_te_enrichment() into a standardized
#' tibble conforming to the project's master table schema.
#'
#' @param enrichment_results Output object from compute_te_enrichment()
#' @param alpha Significance threshold for final filtering (optional)
#'
#' @return A tibble with columns:
#'   - te_group: TE group name
#'   - te_level: Level of grouping (subfamily, family, class)
#'   - te_class: TE class (LINE, SINE, etc.)
#'   - replication_type: cytoplasmic or nuclear
#'   - contrast: Contrast name
#'   - pvalue: Raw p-value
#'   - padj: Adjusted p-value
#'   - enrichment_score: Mean t-statistic (for visualization size)
#'   - avg_logfc: Average logFC of group members
#'   - direction: up, down, or ns
#'   - n_elements: Number of TEs in the group
#'   - prop_significant: Proportion of group members individually significant
#'   - analysis_mode: combined or separate
#'
#' @export
format_te_enrichment <- function(enrichment_results, alpha = 0.05) {

    # Extract components
    results <- enrichment_results$results
    meta <- enrichment_results$metadata

    # Validate input
    required_cols <- c("te_group", "te_level", "contrast", "pval_raw",
                       "t_mean", "avg_logfc", "direction")
    if (!all(required_cols %in% colnames(results))) {
        stop("Input results missing required columns. Run compute_te_enrichment() first.")
    }

    # Select and rename columns to match schema
    formatted <- results %>%
        dplyr::select(
            te_group,
            te_level,
            te_class,
            replication_type,
            contrast,
            pvalue = pval_raw,
            padj,
            enrichment_score = t_mean,
            avg_logfc,
            direction,
            n_elements,
            prop_significant
        ) %>%
        dplyr::mutate(
            analysis_mode = meta$analysis_mode,
            # Ensure direction is consistent with alpha if provided
            direction = case_when(
                padj >= alpha ~ "ns",
                avg_logfc > 0 ~ "up",
                avg_logfc < 0 ~ "down",
                TRUE ~ "ns"
            )
        )

    return(formatted)
}
