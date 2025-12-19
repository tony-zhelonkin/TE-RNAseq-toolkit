# format_te_aggregates.R - Format TE aggregates to master table schema
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Formats output from compute_te_aggregates() into a standardized
#              long-format master table.
#
# Dependencies: dplyr, tidyr, tibble

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
})

#' Format TE Aggregates to Master Table Schema
#'
#' Transforms the matrix output of compute_te_aggregates() into a long-format
#' tibble suitable for plotting and export.
#'
#' @param aggregates_results Output object from compute_te_aggregates()
#'
#' @return A tibble with columns:
#'   - te_group: TE group name
#'   - te_level: Level of grouping
#'   - sample_id: Sample identifier
#'   - expression: Expression value (counts, cpm, or logcpm)
#'   - expression_type: Type of expression value
#'   - te_class: TE class (if available)
#'   - analysis_mode: combined or separate
#'
#' @export
format_te_aggregates <- function(aggregates_results) {

    # Extract components
    mat <- aggregates_results$matrix
    annot <- aggregates_results$annotation
    meta <- aggregates_results$metadata

    # Convert matrix to long format
    long_df <- mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column("te_group") %>%
        tidyr::pivot_longer(
            cols = -te_group,
            names_to = "sample_id",
            values_to = "expression"
        )

    # Join with annotation
    formatted <- long_df %>%
        dplyr::left_join(annot, by = "te_group") %>%
        dplyr::mutate(
            te_level = meta$te_level,
            expression_type = meta$value_type,
            analysis_mode = meta$analysis_mode
        ) %>%
        dplyr::select(
            te_group,
            te_level,
            sample_id,
            expression,
            expression_type,
            te_class = class,
            replication_type,
            analysis_mode
        )

    return(formatted)
}
