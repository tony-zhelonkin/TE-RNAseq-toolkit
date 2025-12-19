#' Create multi-panel heatmap of TE-pathway correlations
#'
#' @param correlation_matrices List of correlation matrices from calculate_te_pathway_correlations
#' @param database_name Name of the pathway database (for title)
#' @param output_file Path to save the PDF output
#' @param cluster_rows Whether to cluster rows (pathways)
#' @param cluster_cols Whether to cluster columns (TEs)
#' @param color_palette Color palette to use
#' @param width PDF width
#' @param height PDF height
visualize_te_pathway_correlations <- function(
  correlation_matrices,
  database_name,
  output_file = NULL,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color_palette = "RdBu",
  width = 14,
  height = 10
) {
  # Load required libraries
  library(pheatmap)
  
  # Check if we have any valid data
  has_data <- FALSE
  for (cell in names(correlation_matrices)) {
    for (time in names(correlation_matrices[[cell]])) {
      if (!is.null(correlation_matrices[[cell]][[time]]) && 
          nrow(correlation_matrices[[cell]][[time]]) > 0 && 
          ncol(correlation_matrices[[cell]][[time]]) > 0) {
        has_data <- TRUE
        break
      }
    }
    if (has_data) break
  }
  
  if (!has_data) {
    message("No valid correlation matrices to visualize")
    return(NULL)
  }
  
  # Get cell types and timepoints
  cell_types <- names(correlation_matrices)
  timepoints <- unique(unlist(lapply(correlation_matrices, names)))
  
  # Sort timepoints in logical order if possible
  timepoint_order <- c("0h", "0.5h", "1h", "24h", "96h")
  timepoints <- timepoints[order(match(timepoints, timepoint_order))]
  
  # Create a common color scale for all heatmaps
  all_values <- c()
  for (cell in names(correlation_matrices)) {
    for (time in names(correlation_matrices[[cell]])) {
      if (!is.null(correlation_matrices[[cell]][[time]]) && 
          nrow(correlation_matrices[[cell]][[time]]) > 0 && 
          ncol(correlation_matrices[[cell]][[time]]) > 0) {
        all_values <- c(all_values, as.vector(correlation_matrices[[cell]][[time]]))
      }
    }
  }
  
  if (length(all_values) == 0) {
    message("No valid correlation values found")
    return(NULL)
  }
  
  min_val <- min(all_values, na.rm = TRUE)
  max_val <- max(all_values, na.rm = TRUE)
  max_abs <- max(abs(c(min_val, max_val)))
  
  # Use a symmetric color scale centered at 0
  if (color_palette == "RdBu") {
    colors <- colorRampPalette(c("blue", "white", "red"))(100)
  } else {
    colors <- viridis::viridis(100, option = color_palette)
  }
  
  # Set up the plotting device
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
  }
  
  # Calculate layout dimensions
  n_rows <- length(cell_types)
  n_cols <- length(timepoints)
  
  # Set up layout
  layout_matrix <- matrix(1:(n_rows * n_cols), nrow = n_rows, byrow = TRUE)
  layout(layout_matrix)
  
  # Create each heatmap
  for (i in 1:length(cell_types)) {
    cell <- cell_types[i]
    
    for (j in 1:length(timepoints)) {
      time <- timepoints[j]
      
      if (!time %in% names(correlation_matrices[[cell]]) || 
          is.null(correlation_matrices[[cell]][[time]]) || 
          nrow(correlation_matrices[[cell]][[time]]) == 0 || 
          ncol(correlation_matrices[[cell]][[time]]) == 0) {
        # Create an empty plot if no data
        plot.new()
        title(main = paste(cell, time))
        next
      }
      
      cor_matrix <- correlation_matrices[[cell]][[time]]
      
      # Determine which labels to show
      show_rownames <- (j == length(timepoints))  # Show row names only in rightmost column
      show_colnames <- (i == length(cell_types))  # Show column names only in bottom row
      
      # Create the heatmap
      heatmap_title <- paste(cell, time)
      
      pheatmap(
        cor_matrix,
        main = heatmap_title,
        color = colors,
        breaks = seq(-max_abs, max_abs, length.out = 101),
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = show_rownames,
        show_colnames = show_colnames,
        fontsize_row = 8,
        fontsize_col = 8,
        silent = TRUE
      )
    }
  }
  
  # Add an overall title
  title(main = paste("Correlation between", toupper(database_name), "pathways and TE expression"),
        outer = TRUE, line = -1)
  
  if (!is.null(output_file)) {
    dev.off()
  }
}