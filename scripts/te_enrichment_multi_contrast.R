create_te_enrichment_plot <- function(
    DGE_te, 
    contrast_matrix, 
    fit.t,
    te_level = "subfamily", 
    pval_threshold = 0.05,
    selected_contrasts = NULL,
    contrast_names = NULL,
    custom_title = NULL,
    custom_subtitle = NULL,
    exclude_features = c("gene"),
    exclude_groups = NULL,
    p_adjust_method = "BH",
    adjust_across = c("by_contrast", "all"),
    color_scheme = "plasma", # "default", "viridis", "plasma", "inferno", "blue_yellow", "gray"
    replication_location = c("all", "cytoplasmic", "nuclear")  # NEW PARAMETER

) {  
    # Match arguments
    adjust_across <- match.arg(adjust_across)
    color_scheme <- match.arg(color_scheme)
    replication_location <- match.arg(replication_location)

    # Define TE families by replication location
    cytoplasmic_families <- c("LINE", "SINE", "LTR", "Retroposon")
    nuclear_families <- c("DNA", "RC")
    
    # Default title and subtitle
    loc_text <- if(replication_location != "all") paste0(" (", replication_location, " replication)") else ""
    plot_title <- if(!is.null(custom_title)) custom_title else paste0("TE ", te_level, " Enrichment", loc_text)
    plot_subtitle <- if(!is.null(custom_subtitle)) custom_subtitle else 
        paste("Outline for adj. p <", pval_threshold, "(method:", p_adjust_method, ")")
    
    # Extract TE information
    te_rownames <- rownames(DGE_te$counts)
    parsed <- lapply(te_rownames, parse_te_info)

    # Extract subfamily/family information
    te_groups <- sapply(parsed, "[[", te_level)  # e.g. subfamily
    unique_groups <- unique(te_groups)
    
    # Create index list for TE groups
    set_list <- lapply(unique_groups, function(group) which(te_groups == group))
    names(set_list) <- unique_groups
    
    # Coeff names
    coef_names <- colnames(contrast_matrix)
    
    # Helper to calculate enrichment stats for a (contrast, group)
    calculate_enrichment <- function(cn, group) {
        idx <- set_list[[group]]
        tstats <- fit.t$t[, cn]
        raw_pval <- geneSetTest(index = idx, statistics = tstats, alternative = "either")
        avg_logFC <- mean(fit.t$coefficients[idx, cn])
        direction <- if(avg_logFC > 0) "up" else "down"
        
        list(pval = raw_pval, avg_logFC = avg_logFC, direction = direction)
    }
    
    # 1) Expand to all subfamily × contrast combos
    all_gst_results <- expand.grid(
        te_group  = unique_groups,
        contrast  = coef_names,
        stringsAsFactors = FALSE
    ) %>%
    rowwise() %>%
    mutate(
        results   = list(calculate_enrichment(contrast, te_group)),
        pval      = results$pval,
        avg_logFC = results$avg_logFC,
        direction = results$direction
    ) %>%
    ungroup() %>%
    select(te_group, contrast, pval, avg_logFC, direction)
    
    # 2) First, filter by contrasts and exclude features
    plot_data <- all_gst_results %>%
        { if(!is.null(selected_contrasts)) filter(., contrast %in% selected_contrasts) else . } %>%
        filter(!te_group %in% exclude_features) %>%
        { if(!is.null(exclude_groups)) filter(., !te_group %in% exclude_groups) else . }

    # 3) Apply multiple testing correction (still on all TEs)
    if(adjust_across == "by_contrast") {
        plot_data <- plot_data %>%
          group_by(contrast) %>%
          mutate(pval_adj = p.adjust(pval, method = p_adjust_method)) %>%
          ungroup()
    } else {
        # adjust_across == "all"
        plot_data <- plot_data %>%
          mutate(pval_adj = p.adjust(pval, method = p_adjust_method))
    }
    
    # 4) Add columns for plotting
    plot_data <- plot_data %>%
        mutate(
            celltype = case_when(
                str_starts(contrast, "Th1_")  ~ "Th1",
                str_starts(contrast, "Th17_") ~ "Th17",
                TRUE                          ~ "Other"
            ),
            negLogP = -log10(pval_adj),   # use the adjusted p-value for significance
            sign    = ifelse(avg_logFC > 0, 1, -1)
        )
    
    # 5) If user supplies new names for the contrasts
    if (!is.null(contrast_names) && is.list(contrast_names)) {
        plot_data <- plot_data %>%
            mutate(
                contrast = str_replace_all(contrast, unlist(contrast_names))
            )
    }
    
    # 6) NOW filter by replication location for PLOTTING ONLY
    # First, we need to map each te_group to its family
    if(replication_location != "all") {
        # Create a lookup from te_group to family
        group_to_family <- list()
        for(i in seq_along(parsed)) {
            group <- parsed[[i]][[te_level]]
            family <- parsed[[i]][["family"]]
            group_to_family[[group]] <- family
        }
        
        # Determine which te_groups to keep based on their family
        if(replication_location == "cytoplasmic") {
            keep_groups <- names(group_to_family)[sapply(group_to_family, function(x) x %in% cytoplasmic_families)]
            plot_data <- plot_data %>% filter(te_group %in% keep_groups)
        } else if(replication_location == "nuclear") {
            keep_groups <- names(group_to_family)[sapply(group_to_family, function(x) x %in% nuclear_families)]
            plot_data <- plot_data %>% filter(te_group %in% keep_groups)
        }
    }

    # Choose color scheme based on parameter
    color_scale <- switch(color_scheme,
        "default" = scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0),
        "viridis" = scale_color_viridis_c(option = "viridis", direction = -1),  # Higher values darker
        "plasma" = scale_color_viridis_c(option = "plasma", direction = -1),
        "inferno" = scale_color_viridis_c(option = "inferno", direction = -1),
        "blue_yellow" = scale_color_gradient2(low = "#000066", mid = "white", high = "#EEEE00", midpoint = 0),
        "gray" = scale_color_gradient2(low = "black", mid = "white", high = "darkgray", midpoint = 0)
    )
    
   # 7) Build the ggplot
    p <- ggplot(plot_data, aes(x = contrast, y = te_group)) +
        # base bubble layer
        geom_point(aes(size = negLogP, color = avg_logFC)) +
        
        # black outline for points that pass threshold on adjusted p-value
        geom_point(
            data = subset(plot_data, pval_adj <= pval_threshold),
            aes(size = negLogP),
            shape = 21, fill = NA, color = "black", stroke = 1
        ) +
        
        color_scale +  # Apply the chosen color scale
        scale_size_continuous(range = c(1, 8)) +
        
        # separate facets for Th1 vs. Th17, removing unused x levels
        facet_wrap(~ celltype, scales = "free_x") +
        theme_minimal() +
        labs(
            title    = plot_title,
            subtitle = plot_subtitle,
            x        = "Contrast",
            y        = paste("TE", stringr::str_to_title(te_level)),
            color    = "Mean LogFC",
            size     = bquote("-log10("~italic(P[adj])~")")
        )
    
    # Return both the plot and the final data
    list(plot = p, enrichment_data = plot_data)
}