# plot_te_volcano.R - Volcano plot for TE differential expression
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates volcano plots highlighting TEs against gene background
#              or TE-only volcano plots.
#
# Dependencies: ggplot2, dplyr, ggrepel, plot_utils.R

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
})

#' Create TE Volcano Plot
#'
#' Visualizes differential expression results. Can highlight TEs against a
#' background of genes (if combined analysis) or show TEs only.
#'
#' @param de_results Master table from format_te_de() or combined gene+TE table.
#'                   Must contain: logfc, padj, feature_type (optional).
#' @param highlight_te If TRUE, highlights TEs in color and genes in grey.
#' @param te_labels Character vector of TE IDs to label textually.
#' @param pval_threshold Significance threshold (default: 0.05).
#' @param logfc_threshold LogFC threshold (default: 1.0).
#' @param title Plot title.
#'
#' @return A ggplot object
#'
#' @export
plot_te_volcano <- function(de_results,
                            highlight_te = TRUE,
                            te_labels = NULL,
                            pval_threshold = 0.05,
                            logfc_threshold = 1.0,
                            title = "Volcano Plot") {

    # Check required columns
    required <- c("logfc", "padj")
    if (!all(required %in% colnames(de_results))) {
        stop("de_results missing required columns (logfc, padj)")
    }

    # Prepare data
    plot_data <- de_results %>%
        mutate(
            is_sig = padj < pval_threshold & abs(logfc) > logfc_threshold,
            nlog10p = -log10(padj)
        )

    # Determine colors
    if (highlight_te && "feature_type" %in% colnames(plot_data)) {
        plot_data <- plot_data %>%
            mutate(
                color_group = case_when(
                    !is_sig ~ "NS",
                    feature_type == "te" ~ "TE (Sig)",
                    feature_type == "gene" ~ "Gene (Sig)",
                    TRUE ~ "Sig"
                )
            )
        
        colors <- c(
            "NS" = "grey80",
            "TE (Sig)" = "red",
            "Gene (Sig)" = "grey40",
            "Sig" = "black"
        )
    } else {
        plot_data <- plot_data %>%
            mutate(
                color_group = case_when(
                    !is_sig ~ "NS",
                    logfc > 0 ~ "Up",
                    logfc < 0 ~ "Down"
                )
            )
        
        colors <- c(
            "NS" = "grey80",
            "Up" = "firebrick",
            "Down" = "steelblue"
        )
    }

    # Base plot
    p <- ggplot(plot_data, aes(x = logfc, y = nlog10p)) +
        geom_point(aes(color = color_group), alpha = 0.7, size = 1.5) +
        scale_color_manual(values = colors) +
        geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
                   linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(pval_threshold),
                   linetype = "dashed", color = "grey50") +
        theme_te() +
        labs(
            title = title,
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value",
            color = "Status"
        )

    # Labels
    if (!is.null(te_labels)) {
        label_data <- plot_data %>%
            filter(te_id %in% te_labels | (is_sig & feature_type == "te"))
        
        p <- p + geom_text_repel(
            data = label_data,
            aes(label = te_id),
            size = 3,
            max.overlaps = 20
        )
    }

    return(p)
}
