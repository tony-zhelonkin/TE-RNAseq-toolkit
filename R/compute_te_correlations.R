# compute_te_correlations.R - Compute correlations between TEs and gene pathways
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Computes correlations between aggregated TE expression and 
#              gene pathway activity (e.g., GSEA scores or mean pathway expression).
#
# Dependencies: stats, dplyr, tibble

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

#' Compute TE-Pathway Correlations
#'
#' Calculates Pearson or Spearman correlations between aggregated TE expression 
#' and a matrix of pathway activities (e.g., GSVA scores, Z-scores, or mean expression).
#'
#' @param te_aggregates Output from compute_te_aggregates() or a numeric matrix (groups x samples).
#' @param pathway_matrix Numeric matrix of pathway activities (pathways x samples).
#' @param method Correlation method: "pearson" (default) or "spearman".
#' @param p_adjust_method Method for p-value adjustment (default: "BH").
#'
#' @return A list containing:
#'   \item{cor_matrix}{Matrix of correlation coefficients (TE groups x Pathways)}
#'   \item{p_matrix}{Matrix of raw p-values}
#'   \item{padj_matrix}{Matrix of adjusted p-values}
#'   \item{results_long}{Tibble in long format with all statistics}
#'
#' @export
compute_te_correlations <- function(te_aggregates, 
                                    pathway_matrix, 
                                    method = c("pearson", "spearman"),
                                    p_adjust_method = "BH") {
    
    method <- match.arg(method)
    
    # Extract matrix if needed
    if (is.list(te_aggregates) && "matrix" %in% names(te_aggregates)) {
        te_mat <- te_aggregates$matrix
    } else {
        te_mat <- te_aggregates
    }
    
    # Ensure samples match
    common_samples <- intersect(colnames(te_mat), colnames(pathway_matrix))
    if (length(common_samples) < 3) {
        stop("Insufficient common samples for correlation (minimum 3 needed).")
    }
    
    te_mat <- te_mat[, common_samples, drop=FALSE]
    path_mat <- pathway_matrix[, common_samples, drop=FALSE]
    
    n_te <- nrow(te_mat)
    n_path <- nrow(path_mat)
    
    message("[COR] Computing correlations for ", n_te, " TE groups and ", n_path, " pathways...")
    
    # Initialize matrices
    cor_mat <- matrix(NA, nrow = n_te, ncol = n_path, dimnames = list(rownames(te_mat), rownames(path_mat)))
    p_mat <- matrix(NA, nrow = n_te, ncol = n_path, dimnames = list(rownames(te_mat), rownames(path_mat)))
    
    # Compute correlations
    for (i in 1:n_te) {
        for (j in 1:n_path) {
            test <- cor.test(te_mat[i, ], path_mat[j, ], method = method)
            cor_mat[i, j] <- test$estimate
            p_mat[i, j] <- test$p.value
        }
    }
    
    # Adjust p-values across all tests
    padj_mat <- matrix(p.adjust(p_mat, method = p_adjust_method), 
                       nrow = n_te, ncol = n_path, 
                       dimnames = dimnames(p_mat))
    
    # Convert to long format
    results_long <- expand.grid(te_group = rownames(te_mat), pathway = rownames(path_mat), stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        mutate(
            correlation = as.vector(cor_mat),
            pvalue = as.vector(p_mat),
            padj = as.vector(padj_mat)
        )
    
    list(
        cor_matrix = cor_mat,
        p_matrix = p_mat,
        padj_matrix = padj_mat,
        results_long = results_long,
        metadata = list(
            method = method,
            p_adjust_method = p_adjust_method,
            n_te = n_te,
            n_path = n_path,
            n_samples = length(common_samples)
        )
    )
}
