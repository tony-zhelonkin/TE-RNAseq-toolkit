# plot_te_heatmap.R - Heatmap visualization for TE expression
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates heatmaps from TE expression aggregates.
#              Wraps pheatmap for standardized TE visualization.
#
# Dependencies: pheatmap, RColorBrewer, dplyr, plot_utils.R

suppressPackageStartupMessages({
    library(pheatmap)
    library(dplyr)
    library(RColorBrewer)
})

#' Create TE Expression Heatmap
#'
#' Visualizes TE expression (aggregated or individual) as a heatmap.
#'
#' @param aggregates Output from compute_te_aggregates() or a numeric matrix.
#' @param samples Sample metadata data.frame (optional).
#' @param annotation_col Variables in samples to annotate columns by.
#' @param annotation_row Variables in aggregates$annotation to annotate rows by
#'                       (e.g., "class", "replication_type").
#' @param scale Scale values? "row", "column", or "none" (default: "row").
#' @param cluster_rows Cluster rows? (default: TRUE)
#' @param cluster_cols Cluster columns? (default: TRUE)
#' @param show_rownames Show row names? (default: TRUE)
#' @param show_colnames Show column names? (default: TRUE)
#' @param color_palette Color palette name ("viridis", "magma", "rdbu").
#' @param main_title Plot title.
#'
#' @return A pheatmap object
#'
#' @export
plot_te_heatmap <- function(aggregates,
                            samples = NULL,
                            annotation_col = NULL,
                            annotation_row = c("class"),
                            scale = "row",
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            show_rownames = TRUE,
                            show_colnames = TRUE,
                            color_palette = "viridis",
                            main_title = "TE Expression Heatmap") {

    # Extract matrix
    if (is.list(aggregates) && "matrix" %in% names(aggregates)) {
        mat <- aggregates$matrix
        row_meta <- aggregates$annotation
    } else {
        mat <- aggregates
        row_meta <- NULL
    }

    # Prepare column annotation
    col_annot_df <- NA
    if (!is.null(samples) && !is.null(annotation_col)) {
        # Check intersection
        valid_cols <- intersect(annotation_col, colnames(samples))
        if (length(valid_cols) > 0) {
            col_annot_df <- samples[colnames(mat), valid_cols, drop = FALSE]
        }
    }

    # Prepare row annotation
    row_annot_df <- NA
    if (!is.null(row_meta) && !is.null(annotation_row)) {
        valid_rows <- intersect(annotation_row, colnames(row_meta))
        if (length(valid_rows) > 0) {
            # Need rownames to match matrix
            temp_df <- row_meta %>%
                select(te_group, all_of(valid_rows)) %>%
                as.data.frame()
            rownames(temp_df) <- NULL # Strip existing rownames to satisfy column_to_rownames
            
            row_annot_df <- temp_df %>%
                tibble::column_to_rownames("te_group")
            
            # Ensure order matches matrix
            row_annot_df <- row_annot_df[rownames(mat), , drop = FALSE]
        }
    }

    # Colors
    if (color_palette == "rdbu") {
        colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
    } else if (color_palette == "magma") {
        colors <- viridis::magma(100)
    } else {
        colors <- viridis::viridis(100)
    }

    # Plot
    pheatmap::pheatmap(
        mat,
        scale = scale,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = show_rownames,
        show_colnames = show_colnames,
        annotation_col = col_annot_df,
        annotation_row = row_annot_df,
        color = colors,
        main = main_title,
        border_color = NA
    )
}
