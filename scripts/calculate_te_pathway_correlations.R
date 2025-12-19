#' Calculate correlations between TE expression and pathway scores
#'
#' @param te_expr Matrix of TE expression values (rows=samples, cols=TEs)
#' @param pathway_scores Matrix of pathway scores (rows=samples, cols=pathways)
#' @param metadata Data frame with sample metadata (must include cell type and time columns)
#' @param method Correlation method ("pearson", "spearman", or "kendall")
#' @param cell_type_col Name of the column in metadata containing cell type information
#' @param time_col Name of the column in metadata containing time point information
#' @return List of correlation matrices, one for each celltype-timepoint combination
calculate_te_pathway_correlations <- function(
  te_expr,
  pathway_scores,
  metadata,
  method = "spearman",
  cell_type_col = "CellType",
  time_col = "Time"
) {
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
  
  # Check if metadata contains all common samples
  missing_samples <- common_samples[!common_samples %in% rownames(metadata)]
  if (length(missing_samples) > 0) {
    message(sprintf("Warning: %d samples missing from metadata: %s", 
                   length(missing_samples), 
                   paste(head(missing_samples, 5), collapse=", ")))
    # Only use samples that exist in metadata
    common_samples <- setdiff(common_samples, missing_samples)
    te_expr <- te_expr[, common_samples, drop = FALSE]
    pathway_scores <- pathway_scores[common_samples, , drop = FALSE]
  }
  
  metadata <- metadata[common_samples, , drop = FALSE]
  
  # Check required columns in metadata
  if (!cell_type_col %in% colnames(metadata)) {
    stop(sprintf("Metadata must contain a cell type column named '%s'", cell_type_col))
  }
  if (!time_col %in% colnames(metadata)) {
    stop(sprintf("Metadata must contain a time point column named '%s'", time_col))
  }
  
  # Get unique cell types and timepoints
  cell_types <- unique(metadata[[cell_type_col]])
  timepoints <- unique(metadata[[time_col]])
  
  message(sprintf("Found %d cell types and %d timepoints", length(cell_types), length(timepoints)))
  message(sprintf("Cell types: %s", paste(cell_types, collapse=", ")))
  message(sprintf("Timepoints: %s", paste(timepoints, collapse=", ")))
  
  # Initialize results list
  correlation_matrices <- list()
  
  # Calculate correlations for each cell type and timepoint
  for (cell in cell_types) {
    correlation_matrices[[cell]] <- list()
    
    for (time in timepoints) {
      # Get samples for this cell type and timepoint
      samples <- rownames(metadata)[metadata[[cell_type_col]] == cell & metadata[[time_col]] == time]
      
      if (length(samples) < 3) {
        message(sprintf("Skipping %s-%s: only %d samples (need at least 3 for correlation)", 
                       cell, time, length(samples)))
        next
      }
      
      # Extract expression and scores for these samples
      te_subset <- te_expr[, samples, drop = FALSE]
      pathway_subset <- pathway_scores[samples, , drop = FALSE]
      
      message(sprintf("Calculating correlation for %s-%s: %d TEs x %d pathways using %d samples", 
                     cell, time, nrow(te_subset), ncol(pathway_subset), length(samples)))
      
      # Check for constant values that would cause correlation to fail
      constant_rows <- which(apply(te_subset, 1, function(x) length(unique(x)) == 1))
      if (length(constant_rows) > 0) {
        message(sprintf("Removing %d TEs with constant values", length(constant_rows)))
        te_subset <- te_subset[-constant_rows, , drop = FALSE]
      }
      
      constant_cols <- which(apply(pathway_subset, 2, function(x) length(unique(x)) == 1))
      if (length(constant_cols) > 0) {
        message(sprintf("Removing %d pathways with constant values", length(constant_cols)))
        pathway_subset <- pathway_subset[, -constant_cols, drop = FALSE]
      }
      
      # Skip if no data left
      if (nrow(te_subset) == 0 || ncol(pathway_subset) == 0) {
        message(sprintf("Skipping %s-%s: no valid data after removing constants", cell, time))
        next
      }
      
      # Calculate correlation matrix
      # Pathways in rows, TEs in columns
      tryCatch({
        cor_matrix <- cor(pathway_subset, t(te_subset), method = method, use = "pairwise.complete.obs")
        message(sprintf("Created correlation matrix: %d pathways x %d TEs", 
                       nrow(cor_matrix), ncol(cor_matrix)))
        
        # Store in results
        correlation_matrices[[cell]][[time]] <- cor_matrix
      }, error = function(e) {
        message(sprintf("Error calculating correlation for %s-%s: %s", cell, time, e$message))
      })
    }
  }
  
  return(correlation_matrices)
}