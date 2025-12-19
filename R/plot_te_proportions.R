# plot_te_proportions.R - Visualization for TE proportions
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Visualizes TE proportions (Total vs TE) and TE composition.
#
# Dependencies: ggplot2, dplyr, plot_utils.R

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

#' Plot Total TE Proportion
#'
#' Creates a boxplot or barplot of the proportion of mapped reads assigned to TEs.
#'
#' @param props_output Output from compute_te_proportions().
#' @param group_by Column name to group samples by (e.g., "group", "condition").
#' @param type Plot type: "boxplot" or "barplot" (with error bars).
#' @param color_palette Palette name (passed to get_te_colors).
#'
#' @return A ggplot object.
#'
#' @export
plot_te_total_prop <- function(props_output, group_by = "group", type = "boxplot", color_palette = "qualitative") {
    
    df <- props_output$samples
    
    if (!group_by %in% colnames(df)) {
        stop("group_by variable '", group_by, "' not found in data.")
    }
    
    p <- ggplot(df, aes(x = !!sym(group_by), y = te_prop, fill = !!sym(group_by)))
    
    if (type == "boxplot") {
        p <- p + geom_boxplot(alpha = 0.8, outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.5)
    } else if (type == "barplot") {
        # Summarize for barplot
        summ <- df %>%
            group_by(!!sym(group_by)) %>%
            summarise(
                mean = mean(te_prop),
                sd = sd(te_prop),
                .groups = "drop"
            )
        p <- ggplot(summ, aes(x = !!sym(group_by), y = mean, fill = !!sym(group_by))) +
            geom_col(alpha = 0.8) +
            geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = mean + sd), width = 0.2)
    }
    
    p <- p +
        scale_fill_manual(values = get_te_colors("qualitative", n = length(unique(df[[group_by]])))) +
        scale_y_continuous(labels = scales::percent) +
        theme_te() +
        labs(
            title = "Proportion of Library Mapped to TEs",
            y = "TE Proportion (% of Total Counts)",
            x = group_by
        )
    
    return(p)
}

#' Plot TE Composition (Stacked Bar)
#'
#' Visualizes the composition of the TE fraction (e.g., by Class).
#'
#' @param props_output Output from compute_te_proportions().
#' @param level Level to plot: "class" (currently only class supported by compute function).
#' @param group_by Column name to group/facet by.
#' @param average_by_group If TRUE, averages composition within groups. If FALSE, shows individual samples.
#'
#' @return A ggplot object.
#'
#' @export
plot_te_composition <- function(props_output, level = "class", group_by = NULL, average_by_group = FALSE) {
    
    if (is.null(props_output$te_composition[[level]])) {
        stop("Composition data for level '", level, "' not found.")
    }
    
    df <- props_output$te_composition[[level]]
    
    # If averaging
    if (average_by_group && !is.null(group_by)) {
        df <- df %>%
            group_by(!!sym(group_by), te_class) %>%
            summarise(prop_of_te = mean(prop_of_te), .groups = "drop")
        
        x_var <- group_by
    } else {
        x_var <- "sample_id"
    }
    
    p <- ggplot(df, aes(x = !!sym(x_var), y = prop_of_te, fill = te_class)) +
        geom_col(position = "fill", width = 0.8) +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = get_te_colors("qualitative", n = length(unique(df$te_class)))) +
        theme_te() +
        labs(
            title = paste("TE Composition by", tools::toTitleCase(level)),
            y = "Proportion of TE Fraction",
            x = ifelse(average_by_group, group_by, "Sample"),
            fill = tools::toTitleCase(level)
        )
    
    if (!average_by_group) {
        p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        # Facet if group provided
        if (!is.null(group_by)) {
            p <- p + facet_grid(as.formula(paste("~", group_by)), scales = "free_x", space = "free_x")
        }
    }
    
    return(p)
}
