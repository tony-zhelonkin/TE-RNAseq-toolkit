# plot_te_correlation.R - Visualization for TE-pathway correlations
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates heatmaps for TE-pathway correlation results.
#
# Dependencies: pheatmap, RColorBrewer, dplyr, plot_utils.R

suppressPackageStartupMessages({
    library(pheatmap)
    library(dplyr)
    library(RColorBrewer)
})

#' Plot TE-Pathway Correlation Heatmap
#'
#' Visualizes the correlation matrix between TEs and pathways.
#'
#' @param cor_results Output from compute_te_correlations().
#' @param pval_threshold Adjusted p-value threshold for highlighting significant correlations.
#' @param show_sig_only If TRUE, only shows rows/cols with at least one significant correlation.
#' @param ... Additional arguments passed to pheatmap.
#'
#' @return A pheatmap object.
#'
#' @export
plot_te_correlation <- function(cor_results, 
                                pval_threshold = 0.05, 
                                show_sig_only = FALSE,
                                ...) {
    
    mat <- cor_results$cor_matrix
    padj <- cor_results$padj_matrix
    
    if (show_sig_only) {
        # Keep rows/cols with any sig correlation
        sig_rows <- rowSums(padj < pval_threshold, na.rm=TRUE) > 0
        sig_cols <- colSums(padj < pval_threshold, na.rm=TRUE) > 0
        
        mat <- mat[sig_rows, sig_cols, drop=FALSE]
        padj <- padj[sig_rows, sig_cols, drop=FALSE]
        
        if (nrow(mat) == 0 || ncol(mat) == 0) {
            stop("No significant correlations found to plot with show_sig_only=TRUE")
        }
    }
    
    # Create significance markers (* for p < 0.05)
    display_numbers <- matrix("", nrow = nrow(padj), ncol = ncol(padj))
    display_numbers[padj < 0.05] <- "*"
    display_numbers[padj < 0.01] <- "**"
    display_numbers[padj < 0.001] <- "***"
    
    # Colors
    colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
    
    # Plot
    pheatmap::pheatmap(
        mat,
        color = colors,
        display_numbers = display_numbers,
        number_color = "black",
        fontsize_number = 10,
        main = paste("TE-Pathway Correlations (", cor_results$metadata$method, ")"),
        ...
    )
}
