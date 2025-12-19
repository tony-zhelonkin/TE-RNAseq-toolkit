# create_separate_dge.R - Factory function for separate TE-only DGEList
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates a TE-only DGEList for separate analysis.
#              Supports library size borrowing from gene data to prevent
#              normalization artifacts from global TE derepression.
#
# Dependencies: edgeR, te_utils.R, validate_te_input.R

#' Create Separate TE-Only DGEList
#'
#' Factory function that creates a DGEList containing only TEs for separate
#' analysis. This approach is recommended when:
#' - Gene and TE annotations may overlap
#' - You want genes and TEs as separate hypothesis universes
#' - Different preprocessing was used for genes vs TEs
#'
#' @param te_counts TE count matrix (TEs x samples). Can also be a file path.
#' @param samples Sample metadata data.frame. Must have rownames matching sample IDs,
#'                or a 'sample_id' column.
#' @param te_annotation Optional pre-computed TE annotation from build_te_annotation()
#' @param gene_dge_reference Optional gene DGEList for library size borrowing.
#'                           Strongly recommended to prevent normalization artifacts.
#' @param validate If TRUE (default), run validation checks
#' @param normalize If TRUE (default), run calcNormFactors with TMM method
#' @param filter_low_expr If TRUE, apply filterByExpr. Requires design parameter.
#' @param design Design matrix for filterByExpr (optional)
#'
#' @return A DGEList object with:
#'   - TE count matrix
#'   - $genes slot with TE annotation
#'   - $analysis_mode = "separate"
#'   - Library sizes (borrowed from genes if reference provided)
#'   - TMM normalization factors (if normalize=TRUE)
#'
#' @details
#' The function performs these steps:
#' 1. Loads TE count matrix (from file if path provided)
#' 2. Validates input structure
#' 3. Creates DGEList with TE annotation
#' 4. Optionally borrows library sizes from gene reference
#' 5. Optionally normalizes and filters
#'
#' @section Library Size Borrowing:
#' When \code{gene_dge_reference} is provided, library sizes are copied from
#' the gene DGEList. This is strongly recommended because:
#' - TEs are often globally de-repressed (volatile) or low-count (noisy)
#' - Genes provide a stable reference for sequencing depth
#' - Without borrowing, global TE changes can skew normalization
#'
#' @section When to Use:
#' Use this function when:
#' - Annotations may overlap (conservative approach)
#' - You want genes and TEs as separate hypothesis universes
#' - Different preprocessing was used for genes vs TEs
#' - You cannot validate that exonic TEs were removed
#'
#' Use \code{\link{create_combined_dge}} instead when:
#' - Your TE SAF annotation excludes exonic TEs
#' - You want a unified statistical framework
#'
#' @examples
#' \dontrun{
#' # Recommended: with library size borrowing
#' dge_genes <- readRDS("checkpoints/1.1_dge_normalized.rds")
#' dge_te <- create_separate_te_dge(
#'     te_counts = te_matrix,
#'     samples = sample_metadata,
#'     gene_dge_reference = dge_genes
#' )
#'
#' # Without reference (will show warning)
#' dge_te <- create_separate_te_dge(
#'     te_counts = te_matrix,
#'     samples = sample_metadata
#' )
#'
#' # From file
#' dge_te <- create_separate_te_dge(
#'     te_counts = "path/to/te_counts.txt",
#'     samples = sample_metadata,
#'     gene_dge_reference = dge_genes
#' )
#' }
#'
#' @seealso \code{\link{create_combined_dge}} for combined analysis mode
#' @seealso \code{\link{validate_separate_input}} for validation details
#'
#' @export
create_separate_te_dge <- function(te_counts,
                                   samples,
                                   te_annotation = NULL,
                                   gene_dge_reference = NULL,
                                   validate = TRUE,
                                   normalize = TRUE,
                                   filter_low_expr = FALSE,
                                   design = NULL) {

    # Load required packages
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR package is required. Install with BiocManager::install('edgeR')")
    }

    message("=== Creating Separate TE-Only DGEList ===")

    # -------------------------------------------------------------------------
    # Load TE count matrix if file path provided
    # -------------------------------------------------------------------------
    if (is.character(te_counts) && length(te_counts) == 1) {
        message("[LOAD] Reading TE counts from: ", te_counts)
        te_counts <- load_count_matrix(te_counts)
    }

    # Ensure matrix
    te_counts <- as.matrix(te_counts)

    # -------------------------------------------------------------------------
    # Validation
    # -------------------------------------------------------------------------
    if (validate) {
        validate_separate_input(te_counts, samples, stop_on_error = TRUE)
    }

    # -------------------------------------------------------------------------
    # Prepare sample metadata
    # -------------------------------------------------------------------------
    if (!is.null(samples)) {
        if (is.null(rownames(samples))) {
            if ("sample_id" %in% colnames(samples)) {
                rownames(samples) <- samples$sample_id
            } else if ("ID" %in% colnames(samples)) {
                rownames(samples) <- samples$ID
            }
        }
        # Reorder to match count matrix columns
        samples <- samples[colnames(te_counts), , drop = FALSE]
    }

    # -------------------------------------------------------------------------
    # Build TE annotation
    # -------------------------------------------------------------------------
    message("[ANNOTATE] Building TE annotation...")

    feature_ids <- rownames(te_counts)

    if (is.null(te_annotation)) {
        te_annotation <- build_te_annotation(feature_ids)
    }

    # Create feature annotation data.frame
    genes_df <- data.frame(
        feature_id = feature_ids,
        feature_type = "te",
        is_te = TRUE,
        subfamily = te_annotation$subfamily,
        family = te_annotation$family,
        class = te_annotation$class,
        replication_type = te_annotation$replication_type,
        row.names = feature_ids,
        stringsAsFactors = FALSE
    )

    # -------------------------------------------------------------------------
    # Create DGEList
    # -------------------------------------------------------------------------
    message("[CREATE] Building DGEList...")

    dge <- edgeR::DGEList(
        counts = te_counts,
        samples = samples,
        genes = genes_df
    )

    # Add analysis mode metadata
    dge$analysis_mode <- "separate"
    dge$te_toolkit_version <- "2.0.0"

    # -------------------------------------------------------------------------
    # Library size borrowing (if reference provided)
    # -------------------------------------------------------------------------
    if (!is.null(gene_dge_reference)) {
        message("[LIBSIZE] Borrowing library sizes from gene reference...")

        # Validate reference is a DGEList
        if (!inherits(gene_dge_reference, "DGEList")) {
            stop("gene_dge_reference must be a DGEList object", call. = FALSE)
        }

        # Find common samples
        te_samples <- colnames(dge)
        gene_samples <- colnames(gene_dge_reference)
        common_samples <- intersect(te_samples, gene_samples)

        if (length(common_samples) == 0) {
            stop("No common samples between TE and gene reference DGELists",
                 call. = FALSE)
        }

        if (length(common_samples) < length(te_samples)) {
            missing <- setdiff(te_samples, common_samples)
            warning("[LIBSIZE] ", length(missing), " TE samples not found in gene reference: ",
                    paste(head(missing, 3), collapse = ", "),
                    if (length(missing) > 3) "..." else "",
                    "\nThese samples will use their own library sizes.",
                    call. = FALSE)
        }

        # Borrow library sizes for common samples
        for (s in common_samples) {
            dge$samples[s, "lib.size"] <- gene_dge_reference$samples[s, "lib.size"]
        }

        message("[LIBSIZE] Borrowed library sizes for ", length(common_samples),
                " of ", length(te_samples), " samples")
        dge$lib_size_borrowed <- TRUE
        dge$lib_size_source <- "gene_reference"

    } else {
        warning(
            "[CAUTION] No gene_dge_reference provided.\n",
            "Library sizes will be calculated from TE counts only.\n",
            "This may cause normalization artifacts if TEs are globally ",
            "de-repressed or have low counts.\n",
            "Consider providing a gene DGEList for library size borrowing.",
            call. = FALSE
        )
        dge$lib_size_borrowed <- FALSE
        dge$lib_size_source <- "te_counts"
    }

    # -------------------------------------------------------------------------
    # Normalize
    # -------------------------------------------------------------------------
    if (normalize) {
        message("[NORMALIZE] Calculating TMM normalization factors...")
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        message("[NORMALIZE] Normalization factors range: ",
                round(min(dge$samples$norm.factors), 3), " - ",
                round(max(dge$samples$norm.factors), 3))
    }

    # -------------------------------------------------------------------------
    # Filter low expression (optional)
    # -------------------------------------------------------------------------
    if (filter_low_expr) {
        if (is.null(design)) {
            warning("filter_low_expr=TRUE but no design provided. Skipping filtering.")
        } else {
            message("[FILTER] Applying filterByExpr...")
            keep <- edgeR::filterByExpr(dge, design = design)
            n_before <- nrow(dge)
            dge <- dge[keep, , keep.lib.sizes = FALSE]
            n_after <- nrow(dge)
            message("[FILTER] Kept ", n_after, " of ", n_before, " TEs")
        }
    }

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    message("=== Separate TE DGEList created ===")
    message("  TEs: ", nrow(dge))
    message("  Samples: ", ncol(dge))
    message("  Analysis mode: separate")
    message("  Library sizes: ", ifelse(dge$lib_size_borrowed,
                                       "borrowed from gene reference",
                                       "from TE counts (CAUTION)"))

    # Summary by class
    if (!is.null(dge$genes$class)) {
        te_classes <- table(dge$genes$class)
        message("  TE classes:")
        for (cls in names(te_classes)) {
            message("    - ", cls, ": ", te_classes[[cls]])
        }
    }

    return(dge)
}


#' Create Gene-Only DGEList (for reference)
#'
#' Creates a gene-only DGEList that can be used as a reference for
#' library size borrowing in separate TE analysis.
#'
#' @param gene_counts Gene count matrix (genes x samples)
#' @param samples Sample metadata data.frame
#' @param normalize If TRUE (default), run calcNormFactors with TMM
#' @return A DGEList object suitable for use as gene_dge_reference
#'
#' @examples
#' \dontrun{
#' # Create gene reference, then use for TE analysis
#' dge_genes <- create_gene_dge_reference(gene_counts, samples)
#' dge_te <- create_separate_te_dge(te_counts, samples,
#'                                  gene_dge_reference = dge_genes)
#' }
#'
#' @export
create_gene_dge_reference <- function(gene_counts, samples, normalize = TRUE) {

    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR package is required. Install with BiocManager::install('edgeR')")
    }

    message("[CREATE] Building gene reference DGEList...")

    # Ensure matrix
    gene_counts <- as.matrix(gene_counts)

    # Prepare samples
    if (!is.null(samples)) {
        if (is.null(rownames(samples))) {
            if ("sample_id" %in% colnames(samples)) {
                rownames(samples) <- samples$sample_id
            } else if ("ID" %in% colnames(samples)) {
                rownames(samples) <- samples$ID
            }
        }
        samples <- samples[colnames(gene_counts), , drop = FALSE]
    }

    # Create DGEList
    dge <- edgeR::DGEList(counts = gene_counts, samples = samples)

    # Add feature type annotation
    dge$genes <- data.frame(
        feature_id = rownames(gene_counts),
        feature_type = "gene",
        row.names = rownames(gene_counts)
    )

    # Normalize
    if (normalize) {
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
    }

    message("[DONE] Gene reference: ", nrow(dge), " genes x ", ncol(dge), " samples")

    return(dge)
}
