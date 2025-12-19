# Function to visualise baseline TE expression and enrichment differences 
# (for single contrast only, preferably the contrast designated for baseline differences)
create_te_pair_plot <- function(
    DGE_te,
    fit.t,
    contrast_name = "Th1_vs_Th17_37",
    te_level = "subfamily",
    pval_threshold = 0.05,
    p_adjust_method = "BH",
    exclude_groups = c("Unspecified", "Unknown"),
    replication_location = c("all", "cytoplasmic", "nuclear"),
    color_scheme = "plasma", # "default", "viridis", "plasma", "inferno", "blue_yellow", "gray"
    color_by = c("logFC", "pvalue")  # NEW PARAMETER: control what's used for coloring
) {
    replication_location <- match.arg(replication_location)
    color_scheme <- match.arg(color_scheme)
    color_by <- match.arg(color_by)  # Match the color_by parameter
    
    # Define TE families by replication location
    cytoplasmic_families <- c("LINE", "SINE", "LTR", "Retroposon")
    nuclear_families <- c("DNA", "RC")
    
    # Extract TE information
    te_rownames <- rownames(DGE_te$counts)
    parsed <- lapply(te_rownames, parse_te_info)
    te_groups <- sapply(parsed, "[[", te_level)
    te_families <- sapply(parsed, "[[", "family")
    
    # Get results for the single contrast
    t_stats <- fit.t$t[, contrast_name]
    log_fc <- fit.t$coefficients[, contrast_name]
    p_values <- fit.t$p.value[, contrast_name]
    
    # Create data frame with all results
    results_df <- data.frame(
        te_id = te_rownames,
        te_group = te_groups,
        te_family = te_families,
        t_stat = t_stats,
        log_fc = log_fc,
        p_value = p_values,
        stringsAsFactors = FALSE
    )
    
    # Filter by replication location if specified
    if(replication_location == "cytoplasmic") {
        results_df <- results_df[results_df$te_family %in% cytoplasmic_families, ]
    } else if(replication_location == "nuclear") {
        results_df <- results_df[results_df$te_family %in% nuclear_families, ]
    }
    
    # Exclude specified groups
    if(!is.null(exclude_groups)) {
        results_df <- results_df[!results_df$te_group %in% exclude_groups, ]
    }
    
    # Group by subfamily/family and calculate summary statistics
    group_stats <- results_df %>%
        group_by(te_group) %>%
        summarize(
            mean_logFC = mean(log_fc, na.rm = TRUE),
            mean_t = mean(t_stat, na.rm = TRUE),
            # Calculate p-value using geneSetTest
            enrichment_p = geneSetTest(
                index = which(te_groups %in% te_group),
                statistics = t_stats,
                alternative = "either"
            ),
            n_elements = n()
        ) %>%
        mutate(
            adj_pval = p.adjust(enrichment_p, method = p_adjust_method),
            significant = adj_pval < pval_threshold,
            cell_direction = ifelse(mean_logFC > 0, "Th1", "Th17"),
            neg_log_p = -log10(adj_pval)  # Add negative log p-value for coloring
        )
    
    # Choose color scheme based on parameter
    if(color_by == "logFC") {
        # Original behavior - color by logFC
        color_scale <- switch(color_scheme,
            "default" = scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0),
            "viridis" = scale_color_viridis_c(option = "viridis", direction = -1),
            "plasma" = scale_color_viridis_c(option = "plasma", direction = -1),
            "inferno" = scale_color_viridis_c(option = "inferno", direction = -1),
            "blue_yellow" = scale_color_gradient2(low = "#000066", mid = "white", high = "#EEEE00", midpoint = 0),
            "gray" = scale_color_gradient2(low = "black", mid = "white", high = "darkgray", midpoint = 0)
        )
        
        # Define color mapping
        color_mapping <- aes(size = neg_log_p, color = mean_logFC)
        color_name <- "Mean LogFC"
        
        # For logFC coloring, we want both size and color legends
        size_scale <- scale_size_continuous(range = c(3, 8))
        size_name <- bquote("-log10("~italic(P[adj])~")")
    } else {
        # New behavior - color by p-value
        color_scale <- switch(color_scheme,
            "default" = scale_color_gradient(low = "lightblue", high = "darkblue"),
            "viridis" = scale_color_viridis_c(option = "viridis"),
            "plasma" = scale_color_viridis_c(option = "plasma"),
            "inferno" = scale_color_viridis_c(option = "inferno"),
            "blue_yellow" = scale_color_gradient(low = "yellow", high = "blue"),
            "gray" = scale_color_gradient(low = "lightgray", high = "black")
        )
        
        # Define color mapping
        color_mapping <- aes(size = neg_log_p, color = neg_log_p)
        color_name <- bquote("-log10("~italic(P[adj])~")")
        
        # For p-value coloring, hide the size legend as it's redundant
        size_scale <- scale_size_continuous(range = c(3, 8), guide = "none")
        size_name <- NULL
    }

    # Create plot
    p <- ggplot(group_stats, aes(x = cell_direction, y = te_group)) +
        geom_point(color_mapping) +
        geom_point(
            data = subset(group_stats, significant),
            aes(size = neg_log_p),
            shape = 21, fill = NA, color = "black", stroke = 1
        ) +
        color_scale +
        size_scale +  # Apply the appropriate size scale (with or without guide)
        theme_minimal() +
        labs(
            title = paste("TE", te_level, "Expression: Th1 vs Th17 at Baseline"),
            subtitle = paste("Outline for adj. p <", pval_threshold),
            x = "Higher Expression in",
            y = paste("TE", stringr::str_to_title(te_level)),
            color = color_name,
            size = size_name  # Will be NULL when color_by = "pvalue"
        )
    
    list(plot = p, data = group_stats)
}