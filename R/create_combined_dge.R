# create_combined_dge.R - Factory function for combined gene+TE DGEList
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates a combined DGEList containing both genes and TEs
#              for unified analysis. Use when annotations are non-overlapping.
#
# Dependencies: edgeR, te_utils.R, validate_te_input.R

#' Create Combined Gene+TE DGEList
#'
#' Factory function that creates a DGEList containing both genes and TEs
#' for combined analysis. This is the recommended approach when:
#' - Gene and TE annotations are non-overlapping (exonic TEs removed)
#' - You want a unified statistical framework
#' - You want TMM normalization on the combined library
#'
#' @param gene_counts Gene count matrix (genes x samples). Can also be a file path.
#' @param te_counts TE count matrix (TEs x samples). Can also be a file path.
#' @param samples Sample metadata data.frame. Must have rownames matching sample IDs,
#'                or a 'sample_id' column.
#' @param gene_annotation Optional pre-computed gene annotation data.frame
#' @param te_annotation Optional pre-computed TE annotation from build_te_annotation()
#' @param validate If TRUE (default), run all validation checks before creating DGE
#' @param normalize If TRUE (default), run calcNormFactors with TMM method
#' @param filter_low_expr If TRUE, apply filterByExpr. Requires design parameter.
#' @param design Design matrix for filterByExpr (optional)
#'
#' @return A DGEList object with:
#'   - Combined count matrix (genes + TEs)
#'   - $genes slot with feature_type ("gene" or "te") and TE annotation
#'   - $analysis_mode = "combined"
#'   - TMM normalization factors (if normalize=TRUE)
#'
#' @details
#' The function performs these steps:
#' 1. Loads count matrices (from files if paths provided)
#' 2. Validates no ID collision (critical for combined analysis)
#' 3. Validates sample consistency
#' 4. Validates library size ratios
#' 5. Combines matrices (rbind)
#' 6. Creates DGEList with feature annotation
#' 7. Optionally normalizes and filters
#'
#' @section When to Use:
#' Use this function when:
#' - Your TE SAF annotation excludes exonic TEs (verified with bedtools intersect)
#' - Gene and TE matrices have identical samples
#' - You want genes and TEs in a single statistical framework
#'
#' Use \code{\link{create_separate_te_dge}} instead when:
#' - Annotations may overlap
#' - You want separate hypothesis universes for genes vs TEs
#' - Different preprocessing was used for genes vs TEs
#'
#' @section Normalization:
#' TMM normalization is performed on the combined matrix by default.
#' With ~78k genes and ~1.2k TEs, genes dominate the normalization.
#' This is a feature, not a bug: TE changes are measured relative to
#' the stable gene background.
#'
#' @examples
#' \dontrun{
#' # From matrices
#' dge <- create_combined_dge(
#'     gene_counts = gene_matrix,
#'     te_counts = te_matrix,
#'     samples = sample_metadata
#' )
#'
#' # From files
#' dge <- create_combined_dge(
#'     gene_counts = "path/to/gene_counts.txt",
#'     te_counts = "path/to/te_counts.txt",
#'     samples = sample_metadata
#' )
#'
#' # With filtering
#' design <- model.matrix(~ 0 + group, data = sample_metadata)
#' dge <- create_combined_dge(
#'     gene_counts = gene_matrix,
#'     te_counts = te_matrix,
#'     samples = sample_metadata,
#'     filter_low_expr = TRUE,
#'     design = design
#' )
#' }
#'
#' @seealso \code{\link{create_separate_te_dge}} for separate analysis mode
#' @seealso \code{\link{validate_combined_input}} for validation details
#'
#' @export
create_combined_dge <- function(gene_counts,
                                te_counts,
                                samples,
                                gene_annotation = NULL,
                                te_annotation = NULL,
                                validate = TRUE,
                                normalize = TRUE,
                                filter_low_expr = FALSE,
                                design = NULL) {

    # Load required packages
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR package is required. Install with BiocManager::install('edgeR')")
    }

    message("=== Creating Combined Gene+TE DGEList ===")

    # -------------------------------------------------------------------------
    # Load count matrices if file paths provided
    # -------------------------------------------------------------------------
    if (is.character(gene_counts) && length(gene_counts) == 1) {
        message("[LOAD] Reading gene counts from: ", gene_counts)
        gene_counts <- load_count_matrix(gene_counts)
    }

    if (is.character(te_counts) && length(te_counts) == 1) {
        message("[LOAD] Reading TE counts from: ", te_counts)
        te_counts <- load_count_matrix(te_counts)
    }

    # Ensure matrices
    gene_counts <- as.matrix(gene_counts)
    te_counts <- as.matrix(te_counts)

    # -------------------------------------------------------------------------
    # Validation
    # -------------------------------------------------------------------------
    if (validate) {
        validate_combined_input(gene_counts, te_counts, samples, stop_on_error = TRUE)
    }

    # -------------------------------------------------------------------------
    # Prepare sample metadata
    # -------------------------------------------------------------------------
    # Ensure samples is a data.frame with proper rownames
    if (!is.null(samples)) {
        if (is.null(rownames(samples))) {
            if ("sample_id" %in% colnames(samples)) {
                rownames(samples) <- samples$sample_id
            } else if ("ID" %in% colnames(samples)) {
                rownames(samples) <- samples$ID
            }
        }
        # Reorder to match count matrix columns
        samples <- samples[colnames(gene_counts), , drop = FALSE]
    }

    # -------------------------------------------------------------------------
    # Combine count matrices
    # -------------------------------------------------------------------------
    message("[COMBINE] Merging gene (", nrow(gene_counts), ") and TE (",
            nrow(te_counts), ") features")

    combined_counts <- rbind(gene_counts, te_counts)

    message("[COMBINE] Combined matrix: ", nrow(combined_counts), " features x ",
            ncol(combined_counts), " samples")

    # -------------------------------------------------------------------------
    # Build feature annotation
    # -------------------------------------------------------------------------
    message("[ANNOTATE] Building feature annotation...")

    # Detect feature types
    feature_ids <- rownames(combined_counts)
    feature_types <- c(rep("gene", nrow(gene_counts)), rep("te", nrow(te_counts)))

    # Build TE annotation
    if (is.null(te_annotation)) {
        te_annotation <- build_te_annotation(rownames(te_counts))
    }

    # Create base annotation data.frame
    genes_df <- data.frame(
        feature_id = feature_ids,
        feature_type = feature_types,
        row.names = feature_ids,
        stringsAsFactors = FALSE
    )

    # Add TE-specific columns (NA for genes)
    genes_df$is_te <- feature_types == "te"
    genes_df$subfamily <- NA_character_
    genes_df$family <- NA_character_
    genes_df$class <- NA_character_
    genes_df$replication_type <- NA_character_

    # Fill in TE annotation
    te_rows <- which(genes_df$feature_type == "te")
    genes_df$subfamily[te_rows] <- te_annotation$subfamily
    genes_df$family[te_rows] <- te_annotation$family
    genes_df$class[te_rows] <- te_annotation$class
    genes_df$replication_type[te_rows] <- te_annotation$replication_type

    # Add gene annotation if provided
    if (!is.null(gene_annotation)) {
        gene_rows <- which(genes_df$feature_type == "gene")
        # Merge gene-specific columns
        for (col in setdiff(colnames(gene_annotation), colnames(genes_df))) {
            genes_df[[col]] <- NA
            if (col %in% colnames(gene_annotation)) {
                genes_df[[col]][gene_rows] <- gene_annotation[[col]]
            }
        }
    }

    # -------------------------------------------------------------------------
    # Create DGEList
    # -------------------------------------------------------------------------
    message("[CREATE] Building DGEList...")

    dge <- edgeR::DGEList(
        counts = combined_counts,
        samples = samples,
        genes = genes_df
    )

    # Add analysis mode metadata
    dge$analysis_mode <- "combined"
    dge$te_toolkit_version <- "2.0.0"

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
            message("[FILTER] Kept ", n_after, " of ", n_before, " features (",
                    sum(dge$genes$feature_type == "gene"), " genes, ",
                    sum(dge$genes$feature_type == "te"), " TEs)")
        }
    }

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    message("=== Combined DGEList created ===")
    message("  Features: ", nrow(dge), " (", sum(dge$genes$feature_type == "gene"),
            " genes, ", sum(dge$genes$feature_type == "te"), " TEs)")
    message("  Samples: ", ncol(dge))
    message("  Analysis mode: combined")

    return(dge)
}


#' Load count matrix from file
#'
#' Helper function to load count matrices from various file formats.
#'
#' @param file_path Path to count matrix file
#' @return Matrix with rownames (feature IDs) and colnames (sample IDs)
#'
#' @keywords internal
load_count_matrix <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("File not found: ", file_path, call. = FALSE)
    }

    # Detect format from extension
    ext <- tolower(tools::file_ext(file_path))

    if (ext %in% c("tsv", "txt")) {
        df <- read.delim(file_path, row.names = 1, check.names = FALSE)
    } else if (ext == "csv") {
        df <- read.csv(file_path, row.names = 1, check.names = FALSE)
    } else if (ext == "rds") {
        df <- readRDS(file_path)
    } else {
        # Try tab-delimited as default
        df <- read.delim(file_path, row.names = 1, check.names = FALSE)
    }

    as.matrix(df)
}


#' Quick summary of a combined DGEList
#'
#' @param dge A DGEList object created by create_combined_dge
#' @return Invisible NULL (prints summary)
#'
#' @export
print_dge_summary <- function(dge) {
    cat("=== DGEList Summary ===\n")
    cat("Analysis mode:", ifelse(is.null(dge$analysis_mode), "unknown", dge$analysis_mode), "\n")
    cat("Features:", nrow(dge), "\n")
    cat("  - Genes:", sum(dge$genes$feature_type == "gene", na.rm = TRUE), "\n")
    cat("  - TEs:", sum(dge$genes$feature_type == "te", na.rm = TRUE), "\n")
    cat("Samples:", ncol(dge), "\n")

    if (!is.null(dge$samples$norm.factors)) {
        cat("Normalized: Yes (TMM)\n")
        cat("  Norm factor range:", round(min(dge$samples$norm.factors), 3), "-",
            round(max(dge$samples$norm.factors), 3), "\n")
    } else {
        cat("Normalized: No\n")
    }

    if (!is.null(dge$genes$class)) {
        te_classes <- table(dge$genes$class[dge$genes$feature_type == "te"])
        cat("TE classes:\n")
        for (cls in names(te_classes)) {
            cat("  -", cls, ":", te_classes[[cls]], "\n")
        }
    }

    invisible(NULL)
}
