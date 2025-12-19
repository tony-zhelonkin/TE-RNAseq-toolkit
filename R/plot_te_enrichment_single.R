# plot_te_enrichment_single.R - Single contrast enrichment visualization
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates a detailed visualization for a single contrast,
#              showing the relationship between genes and TEs.
#              (Replaces the old 'te_enrichment_single_contrast.R')
#
# Dependencies: ggplot2, dplyr, plot_utils.R

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

#' Create TE Enrichment Plot for Single Contrast
#'
#' Visualizes TE enrichment for a specific contrast, usually comparing
#' two groups (e.g. Up vs Down regulated) or a barcode plot style.
#' 
#' CURRENT IMPLEMENTATION: Bar plot of Normalized Enrichment Scores (NES)
#' or t-statistics for a single contrast, colored by significance.
#'
#' @param enrichment_table Master table from format_te_enrichment().
#' @param contrast_name The specific contrast to plot.
#' @param pval_threshold Significance threshold.
#' @param te_levels_to_show Filter for specific levels (e.g. "family").
#' 
#' @return A ggplot object.
#'
#' @export
plot_te_enrichment_single <- function(enrichment_table,
                                      contrast_name,
                                      pval_threshold = 0.05,
                                      te_levels_to_show = NULL) {
    
    # Filter data
    df <- enrichment_table %>%
        filter(contrast == contrast_name)
    
    if (!is.null(te_levels_to_show)) {
        df <- df %>% filter(te_level %in% te_levels_to_show)
    }
    
    if (nrow(df) == 0) {
        stop("No data found for contrast: ", contrast_name)
    }
    
    # Create significance category
    df <- df %>%
        mutate(
            Significance = case_when(
                padj < 0.001 ~ "***",
                padj < 0.01 ~ "**",
                padj < 0.05 ~ "*",
                TRUE ~ "ns"
            ),
            is_sig = padj < pval_threshold
        ) %>%
        arrange(desc(enrichment_score)) # Sort by score
    
    # Set factor level for ordering in plot
    df$te_group <- factor(df$te_group, levels = df$te_group)
    
    p <- ggplot(df, aes(x = te_group, y = enrichment_score, fill = is_sig)) +
        geom_col() +
        geom_text(aes(label = Significance), vjust = ifelse(df$enrichment_score > 0, -0.2, 1.2)) +
        scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70"), 
                          labels = c("TRUE" = "Sig", "FALSE" = "NS")) +
        theme_te() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            title = paste("TE Enrichment:", contrast_name),
            y = "Enrichment Score (t-stat)",
            x = "TE Group",
            fill = paste("FDR <", pval_threshold)
        )
    
    return(p)
}
