# compute_te_aggregates.R - Aggregate TE expression by family/class
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Pure computation functions for aggregating TE expression
#              at subfamily, family, or class level. Used for expression
#              summaries and heatmap preparation.
#              NO visualization - see plot_te_heatmap.R for plotting.
#
# Dependencies: edgeR, dplyr, tibble, te_utils.R

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

#' Aggregate TE Expression by Group
#'
#' Aggregates TE expression (counts or logCPM) at the subfamily, family,
#' or class level. This is useful for:
#' - Creating expression summaries
#' - Reducing dimensionality for heatmaps
#' - Computing group-level statistics
#'
#' @param dge A DGEList object with TE annotation in $genes slot.
#'            Must have columns: feature_type, subfamily, family, class
#' @param te_level Level of aggregation: "subfamily", "family", or "class"
#' @param method Aggregation method: "sum" (default) or "mean"
#' @param value_type Type of values to aggregate: "counts" (raw), "cpm",
#'                   or "logcpm" (default)
#' @param prior_count Prior count for logCPM calculation (default: 2)
#' @param exclude_groups Character vector of TE groups to exclude
#'                       (e.g., c("Unspecified", "Unknown"))
#' @param min_features Minimum features required in group (default: 1)
#' @param te_only If TRUE (default), only aggregate TE features.
#'                If FALSE, include all features (genes will be NA groups).
#'
#' @return A list with components:
#'   \item{matrix}{Aggregated expression matrix (groups x samples)}
#'   \item{annotation}{Tibble with group metadata}
#'   \item{metadata}{List with analysis parameters}
#'
#' @details
#' The function performs these steps:
#' 1. Subsets to TE features (if te_only=TRUE)
#' 2. Transforms counts to requested value_type
#' 3. Groups features by te_level
#' 4. Aggregates using method (sum or mean)
#'
#' For heatmaps, "logcpm" with "mean" is recommended as it:
#' - Normalizes for library size
#' - Handles zeros gracefully
#' - Gives interpretable per-element expression
#'
#' For differential expression at group level, use "counts" with "sum"
#' and create a new DGEList for analysis.
#'
#' @section Value Types:
#' - counts: Raw counts (integer). Sum = total reads in group.
#' - cpm: Counts per million. Better for visualization than raw counts.
#' - logcpm: log2(CPM + prior). Recommended for heatmaps.
#'
#' @examples
#' \dontrun{
#' # Load data
#' dge <- readRDS("checkpoints/1.5_te_dge_combined.rds")
#'
#' # Aggregate at family level for heatmap
#' agg <- compute_te_aggregates(
#'     dge = dge,
#'     te_level = "family",
#'     method = "mean",
#'     value_type = "logcpm",
#'     exclude_groups = c("Unspecified", "Unknown")
#' )
#'
#' # Access results
#' agg$matrix        # Expression matrix
#' agg$annotation    # Group metadata
#'
#' # Aggregate counts for group-level DE
#' agg_counts <- compute_te_aggregates(
#'     dge = dge,
#'     te_level = "family",
#'     method = "sum",
#'     value_type = "counts"
#' )
#' }
#'
#' @seealso \code{\link{create_aggregated_dge}} for creating DGEList from aggregates
#' @seealso \code{\link{plot_te_heatmap}} for visualization
#'
#' @export
compute_te_aggregates <- function(dge,
                                  te_level = c("subfamily", "family", "class"),
                                  method = c("sum", "mean"),
                                  value_type = c("logcpm", "cpm", "counts"),
                                  prior_count = 2,
                                  exclude_groups = NULL,
                                  min_features = 1,
                                  te_only = TRUE) {

    # Load edgeR for cpm/logCPM
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR package is required. Install with BiocManager::install('edgeR')")
    }

    # Match arguments
    te_level <- match.arg(te_level)
    method <- match.arg(method)
    value_type <- match.arg(value_type)

    message("=== Computing TE Aggregates (", te_level, " level) ===")

    # -------------------------------------------------------------------------
    # Validate inputs
    # -------------------------------------------------------------------------
    if (!inherits(dge, "DGEList")) {
        stop("dge must be a DGEList object", call. = FALSE)
    }

    required_cols <- c("feature_type", te_level)
    missing_cols <- setdiff(required_cols, colnames(dge$genes))
    if (length(missing_cols) > 0) {
        stop("DGE$genes missing required columns: ",
             paste(missing_cols, collapse = ", "),
             "\nRun create_combined_dge() or create_separate_te_dge() first.",
             call. = FALSE)
    }

    # -------------------------------------------------------------------------
    # Subset to TEs if requested
    # -------------------------------------------------------------------------
    if (te_only) {
        te_mask <- dge$genes$feature_type == "te"
        if (sum(te_mask) == 0) {
            stop("No TE features found in DGEList", call. = FALSE)
        }
        dge_subset <- dge[te_mask, ]
        message("[SUBSET] Using ", nrow(dge_subset), " TE features")
    } else {
        dge_subset <- dge
        message("[SUBSET] Using all ", nrow(dge_subset), " features")
    }

    # -------------------------------------------------------------------------
    # Transform to requested value type
    # -------------------------------------------------------------------------
    message("[TRANSFORM] Computing ", value_type, " values...")

    if (value_type == "counts") {
        expr_matrix <- dge_subset$counts
    } else if (value_type == "cpm") {
        expr_matrix <- edgeR::cpm(dge_subset, log = FALSE)
    } else if (value_type == "logcpm") {
        expr_matrix <- edgeR::cpm(dge_subset, log = TRUE, prior.count = prior_count)
    }

    # -------------------------------------------------------------------------
    # Build group mapping
    # -------------------------------------------------------------------------
    message("[GROUP] Building ", te_level, " groups...")

    annotation <- dge_subset$genes
    groups <- annotation[[te_level]]

    # Handle NAs
    na_count <- sum(is.na(groups))
    if (na_count > 0) {
        message("[FILTER] Removing ", na_count, " features with NA group")
        keep <- !is.na(groups)
        expr_matrix <- expr_matrix[keep, , drop = FALSE]
        annotation <- annotation[keep, , drop = FALSE]
        groups <- groups[keep]
    }

    # Get unique groups
    unique_groups <- unique(groups)

    # Exclude specified groups
    if (!is.null(exclude_groups)) {
        exclude_found <- intersect(unique_groups, exclude_groups)
        if (length(exclude_found) > 0) {
            message("[FILTER] Excluding ", length(exclude_found), " groups: ",
                    paste(exclude_found, collapse = ", "))
            keep <- !(groups %in% exclude_groups)
            expr_matrix <- expr_matrix[keep, , drop = FALSE]
            annotation <- annotation[keep, , drop = FALSE]
            groups <- groups[keep]
            unique_groups <- setdiff(unique_groups, exclude_groups)
        }
    }

    # Filter by minimum features
    group_sizes <- table(groups)
    small_groups <- names(group_sizes)[group_sizes < min_features]
    if (length(small_groups) > 0) {
        message("[FILTER] Removing ", length(small_groups),
                " groups with < ", min_features, " features")
        keep <- !(groups %in% small_groups)
        expr_matrix <- expr_matrix[keep, , drop = FALSE]
        annotation <- annotation[keep, , drop = FALSE]
        groups <- groups[keep]
        unique_groups <- setdiff(unique_groups, small_groups)
    }

    message("[GROUP] Aggregating ", length(unique_groups), " groups from ",
            nrow(expr_matrix), " features")

    # -------------------------------------------------------------------------
    # Aggregate expression
    # -------------------------------------------------------------------------
    message("[AGGREGATE] Using ", method, " aggregation...")

    # Pre-allocate result matrix
    agg_matrix <- matrix(
        NA_real_,
        nrow = length(unique_groups),
        ncol = ncol(expr_matrix),
        dimnames = list(unique_groups, colnames(expr_matrix))
    )

    # Aggregate each group
    for (grp in unique_groups) {
        idx <- which(groups == grp)
        if (length(idx) == 1) {
            agg_matrix[grp, ] <- expr_matrix[idx, ]
        } else {
            if (method == "sum") {
                agg_matrix[grp, ] <- colSums(expr_matrix[idx, , drop = FALSE])
            } else {
                agg_matrix[grp, ] <- colMeans(expr_matrix[idx, , drop = FALSE])
            }
        }
    }

    # -------------------------------------------------------------------------
    # Build group annotation
    # -------------------------------------------------------------------------
    message("[ANNOTATE] Building group annotation...")

    # Get class and replication type for each group
    group_annotation <- annotation %>%
        select(all_of(c(te_level, "class", "replication_type"))) %>%
        distinct() %>%
        rename(te_group = !!sym(te_level))

    # If aggregating at class level, class = te_group
    if (te_level == "class") {
        group_annotation <- group_annotation %>%
            mutate(class = te_group)
    }

    # Count features per group
    feature_counts <- as.data.frame(table(groups))
    colnames(feature_counts) <- c("te_group", "n_features")
    feature_counts$te_group <- as.character(feature_counts$te_group)

    # Merge
    group_annotation <- group_annotation %>%
        left_join(feature_counts, by = "te_group") %>%
        arrange(te_group)

    # Ensure order matches matrix
    group_annotation <- group_annotation[match(rownames(agg_matrix), group_annotation$te_group), ]

    # -------------------------------------------------------------------------
    # Build metadata
    # -------------------------------------------------------------------------
    metadata <- list(
        te_level = te_level,
        method = method,
        value_type = value_type,
        prior_count = if (value_type == "logcpm") prior_count else NA,
        n_groups = nrow(agg_matrix),
        n_samples = ncol(agg_matrix),
        n_features_input = nrow(dge_subset),
        te_only = te_only,
        excluded_groups = exclude_groups,
        min_features = min_features,
        analysis_mode = ifelse(is.null(dge$analysis_mode), "unknown", dge$analysis_mode),
        computed_at = Sys.time(),
        toolkit_version = "2.0.0"
    )

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    message("=== TE Aggregation Complete ===")
    message("  Groups: ", nrow(agg_matrix))
    message("  Samples: ", ncol(agg_matrix))
    message("  Value type: ", value_type)
    message("  Aggregation: ", method)

    # Summary by class
    if ("class" %in% colnames(group_annotation)) {
        class_summary <- table(group_annotation$class)
        message("  Groups by class:")
        for (cls in names(class_summary)) {
            message("    - ", cls, ": ", class_summary[[cls]])
        }
    }

    return(list(
        matrix = agg_matrix,
        annotation = group_annotation,
        metadata = metadata
    ))
}


#' Create DGEList from Aggregated TE Counts
#'
#' Creates a new DGEList from aggregated TE counts, suitable for
#' group-level differential expression analysis.
#'
#' @param dge Original DGEList object (for sample metadata and lib sizes)
#' @param te_level Level of aggregation: "subfamily", "family", or "class"
#' @param exclude_groups Character vector of TE groups to exclude
#' @param min_features Minimum features required in group (default: 3)
#' @param normalize If TRUE (default), recalculate TMM normalization factors
#'
#' @return A DGEList object with aggregated counts
#'
#' @details
#' This function:
#' 1. Aggregates raw counts by sum at the specified level
#' 2. Creates a new DGEList with the aggregated counts
#' 3. Preserves sample metadata from original DGEList
#' 4. Optionally recalculates normalization factors
#'
#' Use this when you want to perform differential expression analysis
#' at the TE family or class level rather than element level.
#'
#' @examples
#' \dontrun{
#' # Create family-level DGEList for DE analysis
#' dge_family <- create_aggregated_dge(
#'     dge = dge,
#'     te_level = "family",
#'     exclude_groups = c("Unspecified")
#' )
#'
#' # Proceed with standard limma-voom workflow
#' design <- model.matrix(~ 0 + group, data = dge_family$samples)
#' v <- voom(dge_family, design)
#' fit <- lmFit(v, design)
#' }
#'
#' @export
create_aggregated_dge <- function(dge,
                                  te_level = c("subfamily", "family", "class"),
                                  exclude_groups = NULL,
                                  min_features = 3,
                                  normalize = TRUE) {

    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR package is required. Install with BiocManager::install('edgeR')")
    }

    te_level <- match.arg(te_level)

    message("=== Creating Aggregated DGEList (", te_level, " level) ===")

    # Aggregate counts
    agg <- compute_te_aggregates(
        dge = dge,
        te_level = te_level,
        method = "sum",
        value_type = "counts",
        exclude_groups = exclude_groups,
        min_features = min_features,
        te_only = TRUE
    )

    # Create annotation data.frame
    genes_df <- data.frame(
        feature_id = rownames(agg$matrix),
        te_level = te_level,
        te_group = agg$annotation$te_group,
        class = agg$annotation$class,
        replication_type = agg$annotation$replication_type,
        n_elements = agg$annotation$n_features,
        row.names = rownames(agg$matrix),
        stringsAsFactors = FALSE
    )

    # Create new DGEList
    dge_agg <- edgeR::DGEList(
        counts = agg$matrix,
        samples = dge$samples,
        genes = genes_df
    )

    # Add metadata
    dge_agg$analysis_mode <- "aggregated"
    dge_agg$aggregation_level <- te_level
    dge_agg$te_toolkit_version <- "2.0.0"

    # Normalize
    if (normalize) {
        message("[NORMALIZE] Calculating TMM normalization factors...")
        dge_agg <- edgeR::calcNormFactors(dge_agg, method = "TMM")
        message("[NORMALIZE] Normalization factors range: ",
                round(min(dge_agg$samples$norm.factors), 3), " - ",
                round(max(dge_agg$samples$norm.factors), 3))
    }

    message("=== Aggregated DGEList created ===")
    message("  Groups: ", nrow(dge_agg))
    message("  Samples: ", ncol(dge_agg))
    message("  Level: ", te_level)

    return(dge_agg)
}


#' Compute Aggregated Expression Statistics
#'
#' Computes summary statistics for aggregated TE expression across samples
#' or conditions. Useful for exploratory analysis and reporting.
#'
#' @param aggregates Output from compute_te_aggregates()
#' @param samples Optional sample metadata for grouping
#' @param group_var Column name in samples for grouping (if NULL, all samples pooled)
#'
#' @return Tibble with summary statistics per TE group
#'
#' @export
summarize_te_aggregates <- function(aggregates, samples = NULL, group_var = NULL) {

    mat <- aggregates$matrix
    annot <- aggregates$annotation

    # Compute overall statistics
    overall_stats <- tibble(
        te_group = rownames(mat),
        mean_expr = rowMeans(mat),
        sd_expr = apply(mat, 1, sd),
        min_expr = apply(mat, 1, min),
        max_expr = apply(mat, 1, max),
        median_expr = apply(mat, 1, median)
    )

    # Merge with annotation
    result <- annot %>%
        left_join(overall_stats, by = "te_group")

    # If group variable provided, add per-group means
    if (!is.null(samples) && !is.null(group_var)) {
        if (!group_var %in% colnames(samples)) {
            warning("group_var '", group_var, "' not found in samples")
        } else {
            groups <- unique(samples[[group_var]])
            for (grp in groups) {
                grp_samples <- rownames(samples)[samples[[group_var]] == grp]
                grp_samples <- intersect(grp_samples, colnames(mat))
                if (length(grp_samples) > 0) {
                    col_name <- paste0("mean_", grp)
                    result[[col_name]] <- rowMeans(mat[, grp_samples, drop = FALSE])
                }
            }
        }
    }

    result
}


#' Quick summary of aggregation results
#'
#' @param aggregates Output from compute_te_aggregates()
#' @return Invisible NULL (prints summary)
#'
#' @export
print_aggregates_summary <- function(aggregates) {

    mat <- aggregates$matrix
    annot <- aggregates$annotation
    meta <- aggregates$metadata

    cat("\n=== TE Aggregates Summary ===\n")
    cat("Level:", meta$te_level, "\n")
    cat("Groups:", meta$n_groups, "\n")
    cat("Samples:", meta$n_samples, "\n")
    cat("Value type:", meta$value_type, "\n")
    cat("Aggregation:", meta$method, "\n")
    cat("Input features:", meta$n_features_input, "\n")
    cat("Analysis mode:", meta$analysis_mode, "\n\n")

    # Top expressed groups
    cat("Top 10 groups by mean expression:\n")
    mean_expr <- rowMeans(mat)
    top_groups <- head(sort(mean_expr, decreasing = TRUE), 10)
    for (grp in names(top_groups)) {
        n_feat <- annot$n_features[annot$te_group == grp]
        cat(sprintf("  %-30s: %.2f (n=%d)\n", grp, top_groups[grp], n_feat))
    }

    # Summary by class
    if ("class" %in% colnames(annot)) {
        cat("\nGroups by TE class:\n")
        class_counts <- table(annot$class)
        for (cls in names(class_counts)) {
            cat("  -", cls, ":", class_counts[[cls]], "groups\n")
        }
    }

    invisible(NULL)
}
