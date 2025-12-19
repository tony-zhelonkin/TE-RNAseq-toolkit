# compute_te_de.R - Compute TE Differential Expression
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Wrapper function to perform differential expression analysis
#              specifically for TEs (or aggregated TEs). This function handles
#              linear modeling and empirical Bayes moderation.
#
# Dependencies: limma, edgeR

suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
})

#' Compute TE Differential Expression
#'
#' Fits a linear model and performs empirical Bayes moderation to identify
#' differentially expressed TEs. This function is a wrapper around limma::voom
#' (if not already voomed) and limma::lmFit/eBayes.
#'
#' @param dge A DGEList object (either element-level or aggregated).
#' @param design Design matrix for the linear model.
#' @param contrast_matrix Contrast matrix defining comparisons of interest.
#' @param robust Logical, whether to use robust empirical Bayes moderation (default: TRUE).
#' @param voom_normalize Logical, whether to apply voom normalization (default: TRUE).
#'                       Set to FALSE if dge input is already an EList (voom output).
#' @param quality_weights Logical, whether to use array quality weights in voom (default: FALSE).
#'
#' @return A limma MArrayLM fit object (after eBayes).
#'
#' @details
#' This function standardizes the DE workflow for TEs:
#' 1. Checks input object type (DGEList vs EList).
#' 2. Applies voom normalization if needed (transforming counts to logCPM with weights).
#' 3. Fits linear model (lmFit).
#' 4. Applies contrasts (contrasts.fit).
#' 5. Computes statistics (eBayes).
#'
#' This can be used on:
#' - Combined gene+TE DGEList
#' - Separate TE-only DGEList
#' - Aggregated TE DGEList (created by create_aggregated_dge)
#'
#' @examples
#' \dontrun{
#' # Design and contrasts
#' design <- model.matrix(~ 0 + group, data = dge$samples)
#' contrast_matrix <- makeContrasts(Treat-Ctrl, levels=design)
#'
#' # Run DE
#' fit <- compute_te_de(dge, design, contrast_matrix)
#' }
#'
#' @export
compute_te_de <- function(dge, design, contrast_matrix,
                          robust = TRUE,
                          voom_normalize = TRUE,
                          quality_weights = FALSE) {

    message("=== Computing TE Differential Expression ===")

    # 1. Voom transformation (if needed)
    if (voom_normalize) {
        if (inherits(dge, "DGEList")) {
            message("[NORM] Applying voom transformation...")
            if (quality_weights) {
                message("[NORM] Estimating sample quality weights...")
                v <- limma::voomWithQualityWeights(dge, design, plot = FALSE)
            } else {
                v <- limma::voom(dge, design, plot = FALSE)
            }
        } else if (inherits(dge, "EList")) {
            message("[INFO] Input is already an EList (voom output). Skipping voom.")
            v <- dge
        } else {
            stop("Input 'dge' must be a DGEList or EList (voom) object.")
        }
    } else {
        # Assume input is suitable for lmFit directly (e.g. logCPM matrix or EList)
        v <- dge
    }

    # 2. Linear Model Fit
    message("[FIT] Fitting linear model...")
    fit <- limma::lmFit(v, design)

    # 3. Contrasts
    message("[FIT] Applying contrasts...")
    fit <- limma::contrasts.fit(fit, contrast_matrix)

    # 4. Empirical Bayes
    message("[FIT] Computing empirical Bayes statistics (robust=", robust, ")...")
    fit <- limma::eBayes(fit, robust = robust)

    message("=== DE Analysis Complete ===")
    return(fit)
}
