# te_utils.R - Core TE parsing and annotation utilities
# TE-RNAseq-toolkit
# Version: 2.0.0
#
# Description: Functions for parsing TE identifiers and building annotation tables.
#              Supports TEtranscripts-style labels: "subfamily:family:class"
#
# Dependencies: tibble, dplyr, stringr

suppressPackageStartupMessages({
    library(tibble)
    library(dplyr)
    library(stringr)
})

# =============================================================================
# TE CLASSIFICATION CONSTANTS
# =============================================================================

#' TE classes that replicate via cytoplasmic intermediates (retrotransposons)
TE_CYTOPLASMIC_CLASSES <- c("LINE", "SINE", "LTR", "Retroposon")

#' TE classes that replicate in the nucleus (DNA transposons)
TE_NUCLEAR_CLASSES <- c("DNA", "RC")

#' All recognized TE classes
TE_ALL_CLASSES <- c(TE_CYTOPLASMIC_CLASSES, TE_NUCLEAR_CLASSES)

# =============================================================================
# CORE PARSING FUNCTIONS
# =============================================================================

#' Parse a single TE identifier
#'
#' Parses TEtranscripts-style TE labels in the format "subfamily:family:class"
#' (e.g., "HAL1:L1:LINE", "B1_Mus1:Alu:SINE")
#'
#' @param te_id Character string, a single TE identifier
#' @return Named list with components: subfamily, family, class, is_te
#'
#' @examples
#' parse_te_id("HAL1:L1:LINE")
#' # list(subfamily = "HAL1", family = "L1", class = "LINE", is_te = TRUE)
#'
#' parse_te_id("ENSMUSG00000000001")
#' # list(subfamily = NA, family = NA, class = NA, is_te = FALSE)
#'
#' @export
parse_te_id <- function(te_id) {
    parts <- strsplit(te_id, ":", fixed = TRUE)[[1]]

    if (length(parts) >= 3) {
        list(
            subfamily = parts[1],
            family = parts[2],
            class = parts[3],
            is_te = TRUE
        )
    } else {
        list(
            subfamily = NA_character_,
            family = NA_character_,
            class = NA_character_,
            is_te = FALSE
        )
    }
}

#' Alias for parse_te_id (backwards compatibility)
#' @rdname parse_te_id
#' @export
parse_te_info <- parse_te_id

#' Alias for parse_te_id (backwards compatibility)
#' @rdname parse_te_id
#' @export
parse_te_label <- parse_te_id


#' Build TE annotation table from a vector of IDs
#'
#' Parses a vector of feature IDs (genes and/or TEs) and returns a structured
#' annotation table. Non-TE IDs will have NA for TE-specific columns.
#'
#' @param feature_ids Character vector of feature IDs
#' @return Tibble with columns: feature_id, is_te, subfamily, family, class, replication_type
#'
#' @examples
#' ids <- c("HAL1:L1:LINE", "B1_Mus1:Alu:SINE", "ENSMUSG00000000001")
#' build_te_annotation(ids)
#'
#' @export
build_te_annotation <- function(feature_ids) {
    # Parse all IDs
    parsed <- lapply(feature_ids, parse_te_id)

    # Extract components
    tibble(
        feature_id = feature_ids,
        is_te = vapply(parsed, `[[`, logical(1), "is_te"),
        subfamily = vapply(parsed, `[[`, character(1), "subfamily"),
        family = vapply(parsed, `[[`, character(1), "family"),
        class = vapply(parsed, `[[`, character(1), "class")
    ) %>%
        mutate(
            replication_type = case_when(
                class %in% TE_CYTOPLASMIC_CLASSES ~ "cytoplasmic",
                class %in% TE_NUCLEAR_CLASSES ~ "nuclear",
                TRUE ~ NA_character_
            )
        )
}


#' Detect feature type from ID
#'
#' Determines whether a feature ID represents a TE or a gene based on
#' the presence of the "subfamily:family:class" pattern.
#'
#' @param feature_ids Character vector of feature IDs
#' @return Character vector with "te" or "gene" for each ID
#'
#' @examples
#' detect_feature_type(c("HAL1:L1:LINE", "ENSMUSG00000000001"))
#' # c("te", "gene")
#'
#' @export
detect_feature_type <- function(feature_ids) {
    # TE IDs have at least 2 colons (subfamily:family:class)
    n_colons <- str_count(feature_ids, ":")
    ifelse(n_colons >= 2, "te", "gene")
}


#' Check if a feature ID is a TE
#'
#' @param feature_ids Character vector of feature IDs
#' @return Logical vector
#'
#' @export
is_te <- function(feature_ids) {
    detect_feature_type(feature_ids) == "te"
}


#' Get unique TE groups at a specified level
#'
#' @param annotation Annotation tibble from build_te_annotation()
#' @param level One of "subfamily", "family", or "class"
#' @param exclude_groups Character vector of groups to exclude
#' @return Character vector of unique group names
#'
#' @export
get_te_groups <- function(annotation, level = c("subfamily", "family", "class"),
                          exclude_groups = NULL) {
    level <- match.arg(level)

    groups <- annotation %>%
        filter(is_te) %>%
        pull(!!sym(level)) %>%
        unique() %>%
        na.omit()

    if (!is.null(exclude_groups)) {
        groups <- setdiff(groups, exclude_groups)
    }

    sort(groups)
}


#' Create index list for TE groups
#'
#' Creates a named list where each element contains the row indices
#' for features belonging to that TE group. Used for gene set testing.
#'
#' @param annotation Annotation tibble from build_te_annotation()
#' @param level One of "subfamily", "family", or "class"
#' @param exclude_groups Character vector of groups to exclude
#' @return Named list of integer indices
#'
#' @examples
#' # Create index for geneSetTest
#' idx <- create_te_index(annotation, level = "family")
#' geneSetTest(index = idx[["L1"]], statistics = t_stats, alternative = "either")
#'
#' @export
create_te_index <- function(annotation, level = c("subfamily", "family", "class"),
                            exclude_groups = NULL) {
    level <- match.arg(level)

    # Get groups
    groups <- get_te_groups(annotation, level = level, exclude_groups = exclude_groups)

    # Build index list
    index_list <- lapply(groups, function(grp) {
        which(annotation[[level]] == grp)
    })
    names(index_list) <- groups

    index_list
}


#' Map TE groups to their parent class
#'
#' Creates a lookup table from subfamily/family to class.
#'
#' @param annotation Annotation tibble from build_te_annotation()
#' @param from_level Source level ("subfamily", "family", or "class")
#' @return Named character vector (names = groups, values = classes)
#'
#' @export
map_te_to_class <- function(annotation, from_level = c("subfamily", "family", "class")) {
    from_level <- match.arg(from_level)

    if (from_level == "class") {
        # Identity mapping for class level
        classes <- annotation %>%
            filter(is_te) %>%
            pull(class) %>%
            unique() %>%
            na.omit()
        return(setNames(classes, classes))
    }

    annotation %>%
        filter(is_te) %>%
        select(all_of(c(from_level, "class"))) %>%
        distinct() %>%
        {setNames(.[[2]], .[[1]])}
}


#' Map TE groups to replication type
#'
#' @param annotation Annotation tibble from build_te_annotation()
#' @param from_level Source level ("subfamily", "family", or "class")
#' @return Named character vector (names = groups, values = replication_type)
#'
#' @export
map_te_to_replication <- function(annotation, from_level = c("subfamily", "family", "class")) {
    from_level <- match.arg(from_level)

    annotation %>%
        filter(is_te) %>%
        select(all_of(c(from_level, "replication_type"))) %>%
        distinct() %>%
        {setNames(.[[2]], .[[1]])}
}


# =============================================================================
# HELPER FUNCTIONS FOR DGE INTEGRATION
# =============================================================================

#' Add TE annotation to a DGEList's $genes slot
#'
#' @param dge A DGEList object
#' @param annotation Optional pre-computed annotation. If NULL, will be computed.
#' @return DGEList with annotation merged into $genes
#'
#' @export
annotate_dge_features <- function(dge, annotation = NULL) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR package is required for DGEList operations")
    }

    feature_ids <- rownames(dge$counts)

    if (is.null(annotation)) {
        annotation <- build_te_annotation(feature_ids)
    }

    # Initialize or merge with existing $genes
    if (is.null(dge$genes)) {
        dge$genes <- annotation
    } else {
        # Merge, avoiding duplicate columns
        existing_cols <- colnames(dge$genes)
        new_cols <- setdiff(colnames(annotation), c("feature_id", existing_cols))

        if (length(new_cols) > 0) {
            dge$genes <- cbind(dge$genes, annotation[, new_cols, drop = FALSE])
        }
    }

    # Ensure feature_type column exists
    if (!"feature_type" %in% colnames(dge$genes)) {
        dge$genes$feature_type <- ifelse(annotation$is_te, "te", "gene")
    }

    dge
}


#' Subset DGEList to TEs only
#'
#' @param dge A DGEList object with feature_type annotation
#' @return DGEList containing only TE features
#'
#' @export
subset_te <- function(dge) {
    if (is.null(dge$genes$feature_type)) {
        stop("DGEList must have feature_type annotation. Run annotate_dge_features() first.")
    }

    te_idx <- dge$genes$feature_type == "te"
    dge[te_idx, ]
}


#' Subset DGEList to genes only
#'
#' @param dge A DGEList object with feature_type annotation
#' @return DGEList containing only gene features
#'
#' @export
subset_genes <- function(dge) {
    if (is.null(dge$genes$feature_type)) {
        stop("DGEList must have feature_type annotation. Run annotate_dge_features() first.")
    }

    gene_idx <- dge$genes$feature_type == "gene"
    dge[gene_idx, ]
}


# =============================================================================
# SUMMARY FUNCTIONS
# =============================================================================
#' Summarize TE content in a count matrix or DGEList
#'
#' @param x A count matrix or DGEList
#' @return Data frame with summary statistics
#'
#' @export
summarize_te_content <- function(x) {
    if (inherits(x, "DGEList")) {
        counts <- x$counts
        feature_type <- x$genes$feature_type
        if (is.null(feature_type)) {
            feature_type <- detect_feature_type(rownames(counts))
        }
    } else {
        counts <- x
        feature_type <- detect_feature_type(rownames(counts))
    }

    is_te_vec <- feature_type == "te"

    gene_counts <- colSums(counts[!is_te_vec, , drop = FALSE])
    te_counts <- colSums(counts[is_te_vec, , drop = FALSE])
    total_counts <- gene_counts + te_counts
    te_proportion <- te_counts / total_counts

    data.frame(
        sample = colnames(counts),
        gene_counts = gene_counts,
        te_counts = te_counts,
        total_counts = total_counts,
        te_proportion = te_proportion,
        row.names = NULL
    )
}
