# format_te_de.R - Format TE differential expression results
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Extracts and formats TE-specific DE results from a limma fit object.
#
# Dependencies: limma, dplyr, tibble, te_utils.R

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

#' Format TE Differential Expression Results
#'
#' Extracts differential expression statistics for TEs from a limma fit object
#' and formats them into a standardized master table.
#'
#' @param fit A limma MArrayLM fit object (from eBayes)
#' @param dge A DGEList object with TE annotation in $genes slot
#' @param contrasts Character vector of contrasts to extract (NULL = all)
#' @param p_adjust_method Method for p-value adjustment (default: "BH")
#'
#' @return A tibble with columns:
#'   - te_id: Unique TE identifier
#'   - subfamily: TE subfamily
#'   - family: TE family
#'   - te_class: TE class
#'   - contrast: Contrast name
#'   - logfc: Log2 fold change
#'   - pvalue: Raw p-value
#'   - padj: Adjusted p-value
#'   - t_stat: t-statistic
#'   - ave_expr: Average expression (logCPM)
#'   - analysis_mode: combined or separate
#'
#' @export
format_te_de <- function(fit, dge, contrasts = NULL, p_adjust_method = "BH") {

    # Validate inputs
    if (!inherits(fit, "MArrayLM")) {
        stop("fit must be a limma MArrayLM object")
    }

    # Determine contrasts
    available_contrasts <- colnames(fit$coefficients)
    if (is.null(contrasts)) {
        contrasts <- available_contrasts
    } else {
        missing <- setdiff(contrasts, available_contrasts)
        if (length(missing) > 0) {
            stop("Contrasts not found: ", paste(missing, collapse = ", "))
        }
    }

    # Check for TE annotation
    if (!"feature_type" %in% colnames(dge$genes)) {
        stop("DGEList missing feature_type annotation")
    }

    # Identify TEs
    te_idx <- which(dge$genes$feature_type == "te")
    if (length(te_idx) == 0) {
        stop("No TE features found in DGEList")
    }

    results_list <- list()

    message("[EXTRACT] Extracting DE results for ", length(te_idx), " TEs x ",
            length(contrasts), " contrasts...")

    for (contrast in contrasts) {
        # Extract topTable for TEs only
        # We process all TEs, not just top ones
        tt <- limma::topTable(
            fit,
            coef = contrast,
            number = Inf,
            sort.by = "none",
            adjust.method = p_adjust_method
        )

        # Filter to TEs
        # topTable returns row names as feature IDs
        te_tt <- tt[rownames(dge)[te_idx], ]

        # Format
        res <- tibble::tibble(
            te_id = rownames(te_tt),
            subfamily = dge$genes$subfamily[te_idx],
            family = dge$genes$family[te_idx],
            te_class = dge$genes$class[te_idx],
            contrast = contrast,
            logfc = te_tt$logFC,
            pvalue = te_tt$P.Value,
            padj = te_tt$adj.P.Val,
            t_stat = te_tt$t,
            ave_expr = te_tt$AveExpr,
            analysis_mode = ifelse(is.null(dge$analysis_mode), "unknown", dge$analysis_mode)
        )

        results_list[[contrast]] <- res
    }

    # Combine
    combined <- bind_rows(results_list)

    return(combined)
}
