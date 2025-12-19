# plot_utils.R - Shared plotting utilities for TE-RNAseq-toolkit
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Common themes, color palettes, and helper functions for
#              TE visualizations. Ensures consistent aesthetics.
#
# Dependencies: ggplot2, RColorBrewer, viridis

suppressPackageStartupMessages({
    library(ggplot2)
    library(RColorBrewer)
    library(viridis)
})

#' TE Toolkit Theme
#'
#' A standardized ggplot2 theme for publication-quality figures.
#'
#' @param base_size Base font size (default: 12)
#' @param base_family Base font family (default: "sans")
#'
#' @return A ggplot theme object
#'
#' @export
theme_te <- function(base_size = 12, base_family = "sans") {
    theme_minimal(base_size = base_size, base_family = base_family) +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.text = element_text(color = "black"),
            axis.title = element_text(face = "bold"),
            strip.background = element_rect(fill = "grey90", color = NA),
            strip.text = element_text(face = "bold"),
            legend.position = "right",
            plot.title = element_text(face = "bold", hjust = 0),
            plot.subtitle = element_text(hjust = 0)
        )
}

#' Get TE Color Palette
#'
#' Returns standardized color palettes for TE visualizations.
#'
#' @param type Type of palette: "enrichment" (diverging), "heatmap" (sequential),
#'             "qualitative" (categorical)
#' @param name Name of specific palette (optional)
#' @param n Number of colors needed (for qualitative palettes)
#'
#' @return Vector of colors
#'
#' @export
get_te_colors <- function(type = c("enrichment", "heatmap", "qualitative"),
                          name = NULL, n = NULL) {
    type <- match.arg(type)

    if (type == "enrichment") {
        # Diverging palette for logFC or enrichment score
        # Default: Blue-White-Red or similar
        if (is.null(name) || name == "rdbu") {
            return(rev(brewer.pal(11, "RdBu")))
        } else if (name == "piyg") {
            return(rev(brewer.pal(11, "PiYG")))
        } else if (name == "plasma") {
            return(viridis::plasma(11))
        }
    } else if (type == "heatmap") {
        # Sequential palette for expression
        if (is.null(name) || name == "viridis") {
            return(viridis::viridis(100))
        } else if (name == "magma") {
            return(viridis::magma(100))
        } else if (name == "blues") {
            return(colorRampPalette(brewer.pal(9, "Blues"))(100))
        }
    } else if (type == "qualitative") {
        if (is.null(n)) n <- 8
        if (is.null(name) || name == "set2") {
            return(colorRampPalette(brewer.pal(8, "Set2"))(n))
        } else if (name == "paired") {
            return(colorRampPalette(brewer.pal(12, "Paired"))(n))
        }
    }

    # Fallback
    return(viridis::viridis(10))
}

#' Map contrasts to groups
#'
#' Helper to map contrast names to user-defined groups for faceting.
#'
#' @param contrasts Character vector of contrast names
#' @param contrast_groups Named list mapping groups to contrasts
#'                        e.g., list(Time = c("c1", "c2"), Treat = c("c3"))
#'
#' @return Character vector of group names corresponding to contrasts
#'
#' @export
map_contrast_to_group <- function(contrasts, contrast_groups) {
    if (is.null(contrast_groups)) {
        return(rep("All Contrasts", length(contrasts)))
    }

    groups <- rep("Other", length(contrasts))
    names(groups) <- contrasts

    for (grp in names(contrast_groups)) {
        matches <- intersect(contrasts, contrast_groups[[grp]])
        if (length(matches) > 0) {
            groups[matches] <- grp
        }
    }

    unname(groups)
}
