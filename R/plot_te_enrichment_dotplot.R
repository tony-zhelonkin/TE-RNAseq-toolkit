# plot_te_enrichment_dotplot.R - Visualization for TE enrichment
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates dotplots for TE enrichment results.
#              Visualizes enrichment score, significance, and direction
#              across multiple contrasts.
#
# Dependencies: ggplot2, dplyr, plot_utils.R

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

#' Create TE Enrichment Dotplot
#'
#' Visualizes TE enrichment results (from format_te_enrichment) as a dotplot.
#' Point size represents enrichment score (magnitude), color represents
#' significance (or logFC), and shape can represent direction.
#'
#' @param enrichment_table Master table from format_te_enrichment() or tibble/df.
#' @param contrast_groups Optional named list for faceting contrasts
#'                        (e.g., list(Time = c("c1", "c2"))).
#' @param contrast_labels Optional named vector for renaming contrasts in plot.
#' @param te_levels_to_show Character vector of TE levels to include (default: all).
#' @param pval_threshold Threshold for significance (default: 0.05).
#' @param color_by What to map to color: "padj" (significance) or "logfc" (direction).
#' @param size_by What to map to size: "enrichment_score" (t-stat) or "n_elements".
#' @param facet_by Variables to facet by (default: c("te_class", "group")).
#'
#' @return A ggplot object
#'
#' @export
plot_te_enrichment_dotplot <- function(enrichment_table,
                                       contrast_groups = NULL,
                                       contrast_labels = NULL,
                                       te_levels_to_show = NULL,
                                       pval_threshold = 0.05,
                                       color_by = c("padj", "logfc"),
                                       size_by = c("enrichment_score", "n_elements"),
                                       facet_by = c("te_class", "group")) {

    color_by <- match.arg(color_by)
    size_by <- match.arg(size_by)

    # Filter by level if requested
    if (!is.null(te_levels_to_show)) {
        enrichment_table <- enrichment_table %>%
            filter(te_level %in% te_levels_to_show)
    }

    # Add contrast grouping
    if (!is.null(contrast_groups)) {
        enrichment_table <- enrichment_table %>%
            mutate(group = map_contrast_to_group(contrast, contrast_groups))
    } else {
        enrichment_table$group <- "Comparisons"
    }

    # Rename contrasts if labels provided
    if (!is.null(contrast_labels)) {
        # Check coverage
        missing <- setdiff(unique(enrichment_table$contrast), names(contrast_labels))
        if (length(missing) > 0) {
            warning("Some contrasts missing labels: ", paste(missing, collapse = ", "))
        }
        enrichment_table <- enrichment_table %>%
            mutate(contrast_display = coalesce(contrast_labels[contrast], contrast))
    } else {
        enrichment_table$contrast_display <- enrichment_table$contrast
    }

    # Prepare data for plotting
    plot_data <- enrichment_table %>%
        mutate(
            is_sig = padj < pval_threshold,
            # Create a significance score for coloring (-log10 padj * sign)
            sig_score = -log10(padj) * sign(avg_logfc),
            # Cap extreme values for better visualization
            sig_score = pmax(pmin(sig_score, 10), -10)
        )

    # Base plot
    p <- ggplot(plot_data, aes(x = contrast_display, y = te_group))

    # Add points
    if (color_by == "padj") {
        p <- p + geom_point(aes(size = !!sym(size_by), color = sig_score)) +
            scale_color_gradient2(
                low = "blue", mid = "white", high = "red",
                midpoint = 0, name = "Signed -log10(FDR)"
            )
    } else {
        p <- p + geom_point(aes(size = !!sym(size_by), color = avg_logfc)) +
            scale_color_gradient2(
                low = "blue", mid = "white", high = "red",
                midpoint = 0, name = "Avg LogFC"
            )
    }

    # Add significance outline (if using logfc color) or shape
    if (color_by == "logfc") {
        p <- p + geom_point(data = filter(plot_data, is_sig),
                            aes(size = !!sym(size_by)),
                            shape = 21, color = "black", fill = NA, stroke = 1)
    }

    # Faceting
    if (!is.null(facet_by)) {
        # Only facet by variables present in data
        facet_vars <- intersect(facet_by, colnames(plot_data))
        if (length(facet_vars) > 0) {
            facet_formula <- as.formula(paste("~", paste(facet_vars, collapse = "+")))
            p <- p + facet_grid(facet_formula, scales = "free", space = "free")
        }
    }

    # Styling
    p <- p +
        theme_te() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank()
        ) +
        labs(
            title = "TE Enrichment Analysis",
            x = "Contrast",
            y = "TE Group",
            size = ifelse(size_by == "enrichment_score", "Enrichment Score", "Count")
        )

    return(p)
}
