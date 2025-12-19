#' Create correlation heatmaps between TEs and pathways for all databases
#'
#' @param te_expr Matrix of TE expression (rows=samples, cols=TEs)
#' @param gsea_results Results from run_pooled_gsea
#' @param metadata Sample metadata with cell type and time columns
#' @param output_dir Directory to save output files
#' @param correlation_method Correlation method to use
#' @param cluster_rows Whether to cluster rows (pathways)
#' @param cluster_cols Whether to cluster columns (TEs)
#' @param cell_type_col Name of the column in metadata containing cell type information
#' @param time_col Name of the column in metadata containing time point information
create_all_te_pathway_correlations <- function(
  te_expr,
  gsea_results,
  metadata,
  output_dir = "3_Results/TE_Pathway_Correlations",
  correlation_method = "spearman",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  cell_type_col = "CellType",
  time_col = "Time"
) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get list of databases
  databases <- names(gsea_results$scores)
  
  # Process each database
  for (db in databases) {
    message(sprintf("Processing %s database...", toupper(db)))
    
    # Skip if no data
    if (is.null(gsea_results$scores[[db]]) || ncol(gsea_results$scores[[db]]) == 0) {
      message(sprintf("Skipping %s: no data available", db))
      next
    }
    
    # Calculate correlations
    tryCatch({
      correlations <- calculate_te_pathway_correlations(
        te_expr = te_expr,
        pathway_scores = gsea_results$scores[[db]],
        metadata = metadata,
        method = correlation_method,
        cell_type_col = cell_type_col,
        time_col = time_col
      )
      
      # Check if we have any valid correlation matrices
      has_data <- FALSE
      for (cell in names(correlations)) {
        for (time in names(correlations[[cell]])) {
          if (!is.null(correlations[[cell]][[time]]) && 
              nrow(correlations[[cell]][[time]]) > 0 && 
              ncol(correlations[[cell]][[time]]) > 0) {
            has_data <- TRUE
            break
          }
        }
        if (has_data) break
      }
      
      if (!has_data) {
        message(sprintf("Skipping %s: no valid correlation matrices generated", db))
        next
      }
      
      # Create visualization
      output_file <- file.path(output_dir, paste0("te_", db, "_correlations.pdf"))
      visualize_te_pathway_correlations(
        correlation_matrices = correlations,
        database_name = db,
        output_file = output_file,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols
      )
      
      message(sprintf("Saved correlation heatmap to %s", output_file))
    }, error = function(e) {
      message(sprintf("Error processing %s database: %s", db, e$message))
    })
  }
}