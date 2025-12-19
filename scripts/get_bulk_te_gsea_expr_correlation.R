#' Create a single bulk correlation heatmap between TEs and pathways
#'
#' @param te_expr Matrix of TE expression (rows=samples or TEs, cols=TEs or samples)
#' @param pathway_scores Matrix of pathway scores (rows=samples, cols=pathways)
#' @param database_name Name of the database for the title
#' @param output_file File path to save the PDF
#' @param method Correlation method
#' @param color_palette Color palette to use (default: "viridis")
#' @param width PDF width
#' @param height PDF height
#' @param db_type Database type for pathway name cleaning
create_bulk_correlation_heatmap <- function(
  te_expr,
  pathway_scores,
  database_name,
  output_file = NULL,
  method = "spearman",
  color_palette = "viridis",
  width = 12,
  height = 10,
  db_type = NULL
) {
  # Load required packages
  if (!requireNamespace("viridis", quietly = TRUE)) {
    install.packages("viridis")
    library(viridis)
  }
  
  # Check if te_expr needs to be transposed
  # If samples are in rows, we need to transpose
  if (all(rownames(te_expr) %in% rownames(pathway_scores))) {
    message("TE expression matrix has samples in rows - transposing to get TEs in rows")
    te_expr <- t(te_expr)
  }
  
  # Check input dimensions
  message(sprintf("TE expression matrix: %d TEs x %d samples", nrow(te_expr), ncol(te_expr)))
  message(sprintf("Pathway scores matrix: %d samples x %d pathways", nrow(pathway_scores), ncol(pathway_scores)))
  
  # Ensure sample order matches between matrices
  common_samples <- intersect(colnames(te_expr), rownames(pathway_scores))
  message(sprintf("Found %d common samples between matrices", length(common_samples)))
  
  if (length(common_samples) == 0) {
    stop("No common samples between TE expression and pathway scores matrices")
  }
  
  te_expr <- te_expr[, common_samples, drop = FALSE]
  pathway_scores <- pathway_scores[common_samples, , drop = FALSE]
  
  # Check for constant values that would cause correlation to fail
  constant_rows <- which(apply(te_expr, 1, function(x) length(unique(x)) == 1))
  if (length(constant_rows) > 0) {
    message(sprintf("Removing %d TEs with constant values", length(constant_rows)))
    te_expr <- te_expr[-constant_rows, , drop = FALSE]
  }
  
  constant_cols <- which(apply(pathway_scores, 2, function(x) length(unique(x)) == 1))
  if (length(constant_cols) > 0) {
    message(sprintf("Removing %d pathways with constant values", length(constant_cols)))
    pathway_scores <- pathway_scores[, -constant_cols, drop = FALSE]
  }
  
  # Calculate correlation matrix
  message("Calculating correlation matrix...")
  # Transpose the correlation matrix to get pathways as rows and TEs as columns
  cor_matrix <- cor(pathway_scores, t(te_expr), method = method, use = "pairwise.complete.obs")
  message(sprintf("Created correlation matrix: %d pathways x %d TEs", 
                 nrow(cor_matrix), ncol(cor_matrix)))
  
  # Clean pathway names based on database type
  if (!is.null(db_type)) {
    message(sprintf("Cleaning pathway names for database type: %s", db_type))
    clean_names <- rownames(cor_matrix)
    
    if (db_type == "hallmark") {
      clean_names <- gsub("^HALLMARK_", "", clean_names)
    } else if (db_type == "gobp") {
      clean_names <- gsub("^GOBP_", "", clean_names)
    } else if (db_type == "gomf") {
      clean_names <- gsub("^GOMF_", "", clean_names)
    } else if (db_type == "gocc") {
      clean_names <- gsub("^GOCC_", "", clean_names)
    } else if (db_type == "kegg") {
      clean_names <- gsub("^KEGG_MEDICUS_", "", clean_names)
    } else if (db_type == "reactome") {
      clean_names <- gsub("^REACTOME_", "", clean_names)
    } else if (db_type == "biocarta") {
      clean_names <- gsub("^BIOCARTA_", "", clean_names)
    } else if (db_type == "grtd") {
      clean_names <- gsub("^GTRD_", "", clean_names)
    } else if (db_type == "wiki") {
      clean_names <- gsub("^WP_", "", clean_names)
    } else if (db_type == "perturb") {
      clean_names <- gsub("^CGP_", "", clean_names)
    }
    
    # Replace underscores with spaces for better readability
    clean_names <- gsub("_", " ", clean_names)
    rownames(cor_matrix) <- clean_names
  }
  
  # Determine color scale
  max_abs <- max(abs(cor_matrix), na.rm = TRUE)
  
  # Use a proper diverging color palette
  if (color_palette == "viridis") {
    # Create a diverging palette with blue for negative, yellow for positive
    neg_colors <- viridis::viridis(50, option = "viridis", begin = 0, end = 0.5)
    pos_colors <- viridis::viridis(50, option = "viridis", begin = 0.5, end = 1)
    colors <- c(neg_colors, pos_colors)
  } else if (color_palette == "plasma") {
    # Create a diverging palette with blue for negative, orange/red for positive
    neg_colors <- viridis::plasma(50, begin = 0, end = 0.5)
    pos_colors <- viridis::plasma(50, begin = 0.5, end = 1)
    colors <- c(neg_colors, pos_colors)
  } else if (color_palette == "inferno") {
    # Create a diverging palette with purple for negative, yellow for positive
    neg_colors <- viridis::inferno(50, begin = 0, end = 0.5)
    pos_colors <- viridis::inferno(50, begin = 0.5, end = 1)
    colors <- c(neg_colors, pos_colors)
  } else if (color_palette == "magma") {
    # Create a diverging palette with purple for negative, yellow for positive
    neg_colors <- viridis::magma(50, begin = 0, end = 0.5)
    pos_colors <- viridis::magma(50, begin = 0.5, end = 1)
    colors <- c(neg_colors, pos_colors)
  } else if (color_palette == "cividis") {
    # Create a diverging palette with blue for negative, yellow for positive
    neg_colors <- viridis::cividis(50, begin = 0, end = 0.5)
    pos_colors <- viridis::cividis(50, begin = 0.5, end = 1)
    colors <- c(neg_colors, pos_colors)
  } else if (color_palette == "RdBu") {
    # Use RColorBrewer's RdBu palette (red for negative, blue for positive)
    colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
    colors <- colorRampPalette(colors)(100)
  } else if (color_palette == "PiYG") {
    # Use RColorBrewer's PiYG palette (purple for negative, green for positive)
    colors <- rev(RColorBrewer::brewer.pal(11, "PiYG"))
    colors <- colorRampPalette(colors)(100)
  } else if (color_palette == "PRGn") {
    # Use RColorBrewer's PRGn palette (purple for negative, green for positive)
    colors <- rev(RColorBrewer::brewer.pal(11, "PRGn"))
    colors <- colorRampPalette(colors)(100)
  } else {
    # Default to a blue-white-red diverging palette
    colors <- colorRampPalette(c("blue", "white", "red"))(100)
    message(sprintf("Using default blue-white-red palette for '%s'", color_palette))
  }

  
  # Create heatmap with TEs as columns (displayed at bottom with 45-degree angle)
  message("Creating heatmap...")
  p <- pheatmap(
    cor_matrix,
    main = paste("Correlation between", toupper(database_name), "pathways and TE expression"),
    color = colors,
    breaks = seq(-max_abs, max_abs, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,         # Set column names (TEs) to 45-degree angle
    display_numbers = FALSE, # Don't show correlation values in cells
    silent = TRUE
  )
  
  # Save to file if requested
  if (!is.null(output_file)) {
    message(sprintf("Saving heatmap to %s", output_file))
    
    # Adjust width and height based on number of TEs and pathways
    adjusted_width = max(width, ncol(cor_matrix) * 0.15 + 5)
    adjusted_height = max(height, nrow(cor_matrix) * 0.15 + 5)
    
    pdf(output_file, width = adjusted_width, height = adjusted_height)
    print(p)
    dev.off()
  }
  
  return(p)
}