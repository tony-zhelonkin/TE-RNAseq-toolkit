# compute_te_proportions.R - Compute TE proportions in library
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Computes the proportion of reads assigned to TEs vs genes,
#              and the composition of the TE fraction itself.
#
# Dependencies: edgeR, dplyr, tibble, te_utils.R

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

#' Compute TE Proportions
#'
#' Calculates the proportion of library size attributed to TEs versus genes,
#' and optionally breaks down the TE fraction by class, family, or subfamily.
#'
#' @param dge A DGEList object.
#' @param group_by Metadata column in dge$samples to group summaries by (optional).
#'
#' @return A list containing:
#'   - \code{samples}: Tibble with per-sample proportions (total_counts, gene_counts, te_counts, te_prop).
#'   - \code{groups}: Tibble with mean proportions per group (if group_by provided).
#'   - \code{te_composition}: List of tibbles showing TE class/family composition.
#'
#' @export
compute_te_proportions <- function(dge, group_by = NULL) {

    # Validate input
    if (!inherits(dge, "DGEList")) {
        stop("dge must be a DGEList object")
    }
    
    # Ensure feature type annotation exists
    if (is.null(dge$genes$feature_type)) {
        # Try to detect if not present
        dge$genes$feature_type <- detect_feature_type(rownames(dge))
    }
    
    # 1. Calculate Per-Sample TE vs Gene Proportions
    # -----------------------------------------------
    counts <- dge$counts
    is_te <- dge$genes$feature_type == "te"
    
    gene_counts <- colSums(counts[!is_te, , drop=FALSE])
    te_counts   <- colSums(counts[is_te, , drop=FALSE])
    total_counts <- gene_counts + te_counts
    
    sample_props <- tibble::tibble(
        sample_id = colnames(counts),
        total_counts = total_counts,
        gene_counts = gene_counts,
        te_counts = te_counts,
        te_prop = te_counts / total_counts,
        gene_prop = gene_counts / total_counts
    )
    
    # Add sample metadata
    if (!is.null(dge$samples)) {
        # Avoid duplicating columns
        meta_cols <- setdiff(colnames(dge$samples), colnames(sample_props))
        if ("sample_id" %in% colnames(dge$samples)) meta_cols <- setdiff(meta_cols, "sample_id")
        
        # Merge safely by position since dge$samples rows match counts columns
        meta_df <- dge$samples %>% 
            as.data.frame() %>%
            rownames_to_column("sample_id_match") %>%
            select(all_of(c("sample_id_match", meta_cols)))
        
        sample_props <- sample_props %>%
            bind_cols(meta_df %>% select(-sample_id_match))
    }
    
    # 2. Group Summary (if requested)
    # -------------------------------
    group_summary <- NULL
    if (!is.null(group_by)) {
        if (group_by %in% colnames(sample_props)) {
            group_summary <- sample_props %>%
                group_by(across(all_of(group_by))) %>%
                summarise(
                    mean_te_prop = mean(te_prop, na.rm=TRUE),
                    sd_te_prop = sd(te_prop, na.rm=TRUE),
                    n_samples = n(),
                    .groups = "drop"
                )
        } else {
            warning("group_by column '", group_by, "' not found in samples.")
        }
    }
    
    # 3. TE Composition (Class Level)
    # -------------------------------
    # Breakdown the TE fraction itself into classes (LINE, SINE, LTR, etc.)
    # We use CPM for this to normalize for library size, or raw counts?
    # Raw counts are better for "proportion of reads".
    
    te_comp_list <- list()
    
    if ("class" %in% colnames(dge$genes)) {
        te_counts_mat <- counts[is_te, , drop=FALSE]
        te_classes <- dge$genes$class[is_te]
        
        # Aggregate by class
        class_counts <- rowsum(te_counts_mat, group = te_classes, reorder = FALSE)
        
        # Normalize by total TE counts per sample (to get % of TE fraction)
        # Handle division by zero if no TEs
        te_totals <- colSums(te_counts_mat)
        te_totals[te_totals == 0] <- 1 
        
        class_props <- t(t(class_counts) / te_totals)
        
        # Convert to tidy format
        te_comp_list$class <- class_props %>%
            as.data.frame() %>%
            rownames_to_column("te_class") %>%
            tidyr::pivot_longer(-te_class, names_to = "sample_id", values_to = "prop_of_te") %>%
            left_join(sample_props %>% select(sample_id, all_of(group_by)), by="sample_id")
    }
    
    list(
        samples = sample_props,
        groups = group_summary,
        te_composition = te_comp_list
    )
}
