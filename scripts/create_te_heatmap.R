library(edgeR)
library(dplyr)
library(pheatmap)

# parse_te_info <- function(name) {...} # my existing parser

create_te_heatmap <- function(
  DGE, 
  grouping = c("subfamily", "family"),
  only_sig = FALSE,
  sig_groups = NULL,
  exclude_groups = NULL,
  do_aggregate = FALSE,
  sample_order = NULL,
  sample_annotation = NULL,
  annotation_colors = NULL,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main_title  = "TE Expression Heatmap",
  value_scale = c("zscore", "raw"),
  expr_transform_fn = NULL,
  # Add these new parameters:
  color_palette = "viridis", # Options: "viridis", "magma", "inferno", "plasma", "cividis", "RdBu", "BrBG"
  color_direction = TRUE,    # Whether to reverse the color scale
  gaps_col        = c(20)

  ...
) {
  # Load required packages if not already loaded
  if (!requireNamespace("viridis", quietly = TRUE)) {
    install.packages("viridis")
    library(viridis)
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
  }


  # 1) Match user choice for grouping
  grouping <- match.arg(grouping)  # "subfamily" or "family"
  
  # 2) Parse rownames to identify TEs, subfamilies/families
  all_rownames <- rownames(DGE$counts)
  parsed_info  <- lapply(all_rownames, parse_te_info)
  
  # group_vec will be either subfamily or family
  group_vec <- sapply(parsed_info, '[[', grouping)
  
  # 3) Compute logCPM (or any other default expression metric you prefer)
  expr_mat <- cpm(DGE, log = TRUE)  # rows=features, cols=samples
  
  # 4) Filter to only significant if requested
  if (only_sig) {
    if (is.null(sig_groups) || length(sig_groups) == 0) {
      stop("only_sig=TRUE but 'sig_groups' is NULL or empty. Provide sig_groups or set only_sig=FALSE.")
    }
    keep_idx <- which(group_vec %in% sig_groups)
    expr_mat <- expr_mat[keep_idx, , drop=FALSE]
    group_vec <- group_vec[keep_idx]
  }
  
  # Exclude user-specified groups if provided
  if (!is.null(exclude_groups) && length(exclude_groups) > 0) {
    exc_idx <- which(group_vec %in% exclude_groups)
    if (length(exc_idx) > 0) {
      keep <- setdiff(seq_along(group_vec), exc_idx)
      expr_mat <- expr_mat[keep, , drop=FALSE]
      group_vec <- group_vec[keep]
    }
  }
  
  # 5) Optionally aggregate by subfamily/family
  if (do_aggregate) {
    df_expr <- as.data.frame(expr_mat)
    df_expr[[grouping]] <- group_vec
    
    agg_expr <- df_expr %>%
      group_by(.data[[grouping]]) %>%
      summarise(across(where(is.numeric), mean))
    
    final_mat <- as.matrix(agg_expr[,-1])
    rownames(final_mat) <- agg_expr[[grouping]]
  } else {
    final_mat <- expr_mat
    rownames(final_mat) <- rownames(expr_mat)
  }
  
  # 6) Reorder columns if sample_order is given
  if (!is.null(sample_order)) {
    sample_order <- sample_order[sample_order %in% colnames(final_mat)]
    final_mat <- final_mat[, sample_order, drop=FALSE]
    
    if (!is.null(sample_annotation)) {
      sample_annotation <- sample_annotation[sample_order, , drop=FALSE]
    }
  }
  
  # 7) If user provided a custom function, apply it row by row
  #    Otherwise, handle "zscore" or "raw"
  if (!is.null(expr_transform_fn)) {
    # The user might supply something like function(x) x - mean(x)
    # We'll apply it to each row of final_mat
    transform_each_row <- function(r) expr_transform_fn(r)
    final_mat <- t(apply(final_mat, 1, transform_each_row))
    
  } else {
    # No custom function => check value_scale
    value_scale <- match.arg(value_scale)
    if (value_scale == "zscore") {
      # row-wise z-score
      final_mat <- t(scale(t(final_mat)))  
      # subtract row mean, divide by row sd
    } else {
      # "raw" => do nothing special
      # i.e. final_mat remains logCPM
    }
  }
  
  # Determine color palette
  if (color_palette %in% c("viridis", "magma", "inferno", "plasma", "cividis")) {
    colors <- viridis::viridis(100, option = color_palette, direction = ifelse(color_direction, 1, -1))
  } else if (color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    # Use RColorBrewer palettes (many are colorblind-friendly)
    max_colors <- RColorBrewer::brewer.pal.info[color_palette, "maxcolors"]
    colors <- RColorBrewer::brewer.pal(max_colors, color_palette)
    colors <- colorRampPalette(colors)(100)
    if (!color_direction) {
      colors <- rev(colors)
    }
  } else {
    # Default to viridis if palette not recognized
    colors <- viridis::viridis(100)
  }
  
  # 8) Plot heatmap with colorblind-friendly palette
  p <- pheatmap(final_mat,
    annotation_col    = sample_annotation,
    annotation_colors = annotation_colors,
    cluster_rows      = cluster_rows,
    cluster_cols      = cluster_cols,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    main              = main_title,
    color             = colors,  # Use our colorblind-friendly palette
    gaps_col          = gaps_col,
    ...
  )
  
  # 9) Return pheatmap object + final matrix
  return(list(
    pheatmap_obj = p,
    data_matrix  = final_mat
  ))
}