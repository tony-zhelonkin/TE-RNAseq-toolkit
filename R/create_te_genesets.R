# create_te_genesets.R - Generate TERM2GENE format for TE GSEA
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Creates clusterProfiler-compatible gene set collections
#              for running classical GSEA on TE data.
#              - Individual TEs (subfamilies) are treated as "genes"
#              - TE families or classes become "gene sets"
#
# Dependencies: dplyr, tibble

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

# =============================================================================
# MAIN FUNCTION: CREATE TE GENE SETS
# =============================================================================

#' Create TERM2GENE/TERM2NAME for TE GSEA
#'
#' Generates clusterProfiler-compatible gene set collections from TE annotations.
#' Individual TEs (identified by their full ID: subfamily:family:class) serve as
#' "genes", while families or classes serve as "gene sets".
#'
#' @param te_annotation Tibble with TE annotations. Must have columns:
#'   - feature_id: Full TE ID (e.g., "L1MdA_VI:L1:LINE")
#'   - subfamily: TE subfamily name
#'   - family: TE family name
#'   - class: TE class name
#'   - is_te: Logical indicating if feature is a TE
#' @param level Character, aggregation level: "family" or "class"
#' @param use_subfamily_as_gene Logical. If TRUE (default), uses subfamily as
#'   the gene identifier. If FALSE, uses the full feature_id.
#'
#' @return List with components:
#'   - T2G: TERM2GENE data frame (gs_name, gene_symbol)
#'   - T2N: TERM2NAME data frame (gs_name, description)
#'   - stats: Summary statistics
#'
#' @examples
#' # From DGEList
#' genesets <- create_te_genesets(dge_te$genes, level = "family")
#'
#' # Use with clusterProfiler::GSEA()
#' gsea_result <- GSEA(ranked_genes, TERM2GENE = genesets$T2G, TERM2NAME = genesets$T2N)
#'
#' @export
create_te_genesets <- function(te_annotation,
                                level = c("family", "class"),
                                use_subfamily_as_gene = TRUE) {

    level <- match.arg(level)

    # Validate input
    required_cols <- c("is_te", "subfamily", "family", "class")
    missing_cols <- setdiff(required_cols, colnames(te_annotation))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
             "\nUse build_te_annotation() to create proper annotation.")
    }

    # Filter to TEs only
    te_only <- te_annotation %>%
        filter(is_te == TRUE) %>%
        filter(!is.na(subfamily), !is.na(family), !is.na(class))

    if (nrow(te_only) == 0) {
        stop("No valid TE annotations found. Check that is_te column is properly set.")
    }

    # Determine gene identifier
    if (use_subfamily_as_gene) {
        gene_col <- "subfamily"
    } else if ("feature_id" %in% colnames(te_annotation)) {
        gene_col <- "feature_id"
    } else {
        gene_col <- "subfamily"
    }

    # Create TERM2GENE based on level
    if (level == "family") {
        T2G <- te_only %>%
            select(gs_name = family, gene_symbol = !!sym(gene_col)) %>%
            distinct() %>%
            arrange(gs_name, gene_symbol)

        # Create descriptions with class info
        family_to_class <- te_only %>%
            select(family, class) %>%
            distinct() %>%
            group_by(family) %>%
            summarise(classes = paste(unique(class), collapse = "/"), .groups = "drop")

        T2N <- data.frame(
            gs_name = family_to_class$family,
            description = paste0(family_to_class$family, " (", family_to_class$classes, ")"),
            stringsAsFactors = FALSE
        )

    } else if (level == "class") {
        T2G <- te_only %>%
            select(gs_name = class, gene_symbol = !!sym(gene_col)) %>%
            distinct() %>%
            arrange(gs_name, gene_symbol)

        # Create descriptions with replication type
        source("01_modules/TE-RNAseq-toolkit/R/te_utils.R")
        T2N <- data.frame(
            gs_name = unique(T2G$gs_name),
            stringsAsFactors = FALSE
        ) %>%
            mutate(
                replication = case_when(
                    gs_name %in% TE_CYTOPLASMIC_CLASSES ~ "retrotransposon",
                    gs_name %in% TE_NUCLEAR_CLASSES ~ "DNA transposon",
                    TRUE ~ "other"
                ),
                description = paste0(gs_name, " (", replication, ")")
            ) %>%
            select(gs_name, description)
    }

    # Calculate statistics
    set_sizes <- T2G %>%
        group_by(gs_name) %>%
        summarise(size = n(), .groups = "drop")

    stats <- list(
        level = level,
        n_sets = nrow(T2N),
        n_genes = length(unique(T2G$gene_symbol)),
        min_size = min(set_sizes$size),
        max_size = max(set_sizes$size),
        median_size = median(set_sizes$size),
        set_sizes = set_sizes
    )

    message(sprintf("[TE GSEA] Created %s-level gene sets:", level))
    message(sprintf("  Gene sets: %d", stats$n_sets))
    message(sprintf("  Total TEs: %d", stats$n_genes))
    message(sprintf("  Size range: %d - %d (median: %.0f)",
                    stats$min_size, stats$max_size, stats$median_size))

    return(list(
        T2G = as.data.frame(T2G),
        T2N = as.data.frame(T2N),
        stats = stats
    ))
}


#' Create TE gene sets from a DGEList
#'
#' Convenience wrapper that extracts annotations from a DGEList object.
#'
#' @param dge DGEList with TE annotations in $genes slot
#' @param level Character, aggregation level: "family" or "class"
#'
#' @return List with T2G, T2N, and stats (see create_te_genesets)
#'
#' @export
create_te_genesets_from_dge <- function(dge, level = c("family", "class")) {

    if (is.null(dge$genes)) {
        stop("DGEList must have $genes annotation slot. ",
             "Run annotate_dge_features() first.")
    }

    # Check for required columns
    if (!"subfamily" %in% colnames(dge$genes)) {
        # Try to build annotation from rownames
        message("Building TE annotation from rownames...")
        source("01_modules/TE-RNAseq-toolkit/R/te_utils.R")
        annotation <- build_te_annotation(rownames(dge$counts))
    } else {
        annotation <- dge$genes
    }

    create_te_genesets(annotation, level = level)
}


#' Filter TE gene sets by size
#'
#' Removes gene sets that are too small or too large for reliable GSEA.
#'
#' @param te_genesets List from create_te_genesets()
#' @param min_size Minimum number of TEs per set (default: 5)
#' @param max_size Maximum number of TEs per set (default: 500)
#'
#' @return Filtered list with T2G, T2N, and updated stats
#'
#' @export
filter_te_genesets <- function(te_genesets,
                                min_size = 5,
                                max_size = 500) {

    T2G <- te_genesets$T2G
    T2N <- te_genesets$T2N

    # Calculate set sizes
    set_sizes <- T2G %>%
        group_by(gs_name) %>%
        summarise(size = n(), .groups = "drop")

    # Find valid sets
    valid_sets <- set_sizes %>%
        filter(size >= min_size, size <= max_size) %>%
        pull(gs_name)

    n_original <- nrow(T2N)
    n_removed_small <- sum(set_sizes$size < min_size)
    n_removed_large <- sum(set_sizes$size > max_size)

    # Filter
    T2G_filtered <- T2G %>% filter(gs_name %in% valid_sets)
    T2N_filtered <- T2N %>% filter(gs_name %in% valid_sets)

    message(sprintf("[TE GSEA] Filtering gene sets [%d, %d]:", min_size, max_size))
    message(sprintf("  Original: %d sets", n_original))
    message(sprintf("  Too small (<%d): %d removed", min_size, n_removed_small))
    message(sprintf("  Too large (>%d): %d removed", max_size, n_removed_large))
    message(sprintf("  Retained: %d sets (%.1f%%)",
                    length(valid_sets), 100 * length(valid_sets) / n_original))

    # Update stats
    new_sizes <- T2G_filtered %>%
        group_by(gs_name) %>%
        summarise(size = n(), .groups = "drop")

    stats <- list(
        level = te_genesets$stats$level,
        n_sets = length(valid_sets),
        n_genes = length(unique(T2G_filtered$gene_symbol)),
        min_size = min(new_sizes$size),
        max_size = max(new_sizes$size),
        median_size = median(new_sizes$size),
        set_sizes = new_sizes,
        filtering = list(
            min_size = min_size,
            max_size = max_size,
            n_removed_small = n_removed_small,
            n_removed_large = n_removed_large
        )
    )

    return(list(
        T2G = as.data.frame(T2G_filtered),
        T2N = as.data.frame(T2N_filtered),
        stats = stats
    ))
}


#' Export TE gene sets to GMT format
#'
#' Creates a GMT file compatible with GSEA desktop and other tools.
#'
#' @param te_genesets List from create_te_genesets()
#' @param output_file Path to output GMT file
#'
#' @return Invisibly returns the output file path
#'
#' @export
export_te_genesets_gmt <- function(te_genesets, output_file) {

    T2G <- te_genesets$T2G
    T2N <- te_genesets$T2N

    # Get unique gene sets
    gene_sets <- split(T2G$gene_symbol, T2G$gs_name)

    # Build GMT lines
    gmt_lines <- sapply(names(gene_sets), function(gs_name) {
        genes <- gene_sets[[gs_name]]
        desc <- T2N$description[T2N$gs_name == gs_name]
        if (length(desc) == 0) desc <- gs_name
        paste(c(gs_name, desc, genes), collapse = "\t")
    })

    # Write to file
    writeLines(gmt_lines, output_file)

    message(sprintf("Exported %d gene sets to: %s",
                    length(gene_sets), output_file))

    invisible(output_file)
}


message("[TE-RNAseq-toolkit] Loaded: create_te_genesets.R")
