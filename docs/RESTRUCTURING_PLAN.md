# TE-RNAseq-toolkit Restructuring Plan

**Version:** 1.0.0
**Date:** 2025-12-18
**Prepared by:** Claude Code with Gemini Pro 3 preview consultation

---

## Executive Summary

This document outlines the comprehensive restructuring of the TE-RNAseq-toolkit to align with the established RNAseq analysis philosophy. The key changes are:

1. **Separate computation from visualization** (Rule of Three)
2. **Introduce checkpoint caching** via `load_or_compute()` pattern
3. **Remove hardcoded project-specific logic** (no more "Th1"/"Th17")
4. **Create master tables** for cross-language analysis
5. **Support both combined AND separate gene/TE analysis** via Dual-Mode Factory Architecture
6. **Follow phase-based execution model**

---

## Table of Contents

1. [Current Problems](#1-current-problems)
2. [Proposed Directory Structure](#2-proposed-directory-structure)
3. [Dual-Mode Architecture: Combined vs Separate Analysis](#3-dual-mode-architecture-combined-vs-separate-analysis)
4. [The Rule of Three: Function Separation](#4-the-rule-of-three-function-separation)
5. [Master Table Schemas](#5-master-table-schemas)
6. [Checkpoint Strategy](#6-checkpoint-strategy)
7. [Parameterization Strategy](#7-parameterization-strategy)
8. [Implementation Roadmap](#8-implementation-roadmap)
9. [Script Migration Map](#9-script-migration-map)

---

## 1. Current Problems

### 1.1 Mixed Computation and Visualization

**`te_enrichment_multi_contrast.R`** (`create_te_enrichment_plot()`):
- Lines 50-74: Computes enrichment statistics inline
- Lines 145-173: Creates visualization
- **Problem**: Cannot reuse enrichment computation without generating a plot

**`create_te_heatmap.R`**:
- Line 51: Computes logCPM inside visualization function
- **Problem**: Recomputes expression every time instead of using pre-computed data

### 1.2 Hardcoded Project-Specific Logic

```r
# CURRENT (te_enrichment_multi_contrast.R:97-101) - BAD
celltype = case_when(
    str_starts(contrast, "Th1_")  ~ "Th1",
    str_starts(contrast, "Th17_") ~ "Th17",
    TRUE                          ~ "Other"
)
```

This makes the function **unusable for other datasets**.

### 1.3 No Checkpoint Caching

The 1100+ line example workflow (`te_example_workflow.R`) performs expensive operations:
- DE analysis with voom (lines 215-226)
- TE enrichment across many contrasts (lines 250-280)
- Multiple heatmap generations (lines 462-765)

**None are cached** - every run starts from scratch.

### 1.4 No Standardized Output

Results saved as ad-hoc CSVs with no schema:
```r
write.csv(results_sub_all$enrichment_data, "/workspaces/Yasmine-retroT/3_Results/...")
```

### 1.5 No Flexibility for Different Preprocessing Scenarios

The example workflow assumes a single combined matrix, but:
- Some projects may have separate gene/TE matrices
- Strandedness may differ between genes and TEs
- Annotation overlap handling varies by project

---

## 2. Proposed Directory Structure

```
TE-RNAseq-toolkit/
├── README.md                           # Module documentation
├── RESTRUCTURING_PLAN.md               # This document
│
├── R/                                  # Core R functions (installable as package)
│   ├── # --- DATA INGESTION (Factory Functions) ---
│   ├── create_combined_dge.R           # Combined gene+TE DGEList
│   ├── create_separate_dge.R           # Separate TE-only DGEList
│   ├── validate_te_input.R             # Input validation functions
│   │
│   ├── # --- UTILITIES ---
│   ├── te_utils.R                      # Core parsing utilities
│   │
│   ├── # --- PHASE 1: Processing (Compute) ---
│   ├── compute_te_enrichment.R         # TE family/subfamily enrichment
│   ├── compute_te_aggregates.R         # Aggregate expression by group
│   ├── compute_te_correlations.R       # TE-pathway correlations
│   ├── compute_te_proportions.R        # TE proportion calculations
│   ├── compute_te_de.R                 # TE differential expression
│   │
│   ├── # --- PHASE 2: Master Tables (Format) ---
│   ├── format_te_enrichment.R          # Normalize to master table schema
│   ├── format_te_de.R                  # TE DE results to master table
│   ├── format_te_aggregates.R          # Expression aggregates to master table
│   │
│   └── # --- PHASE 3: Visualization (Plot) ---
│       ├── plot_te_enrichment_dotplot.R    # Multi-contrast enrichment dotplot
│       ├── plot_te_enrichment_single.R     # Single contrast paired viz
│       ├── plot_te_heatmap.R               # Expression heatmaps
│       ├── plot_te_correlation.R           # TE-pathway correlation heatmaps
│       ├── plot_te_volcano.R               # Volcano with TE highlighting
│       ├── plot_te_scatter.R               # Scatter plots
│       ├── plot_te_proportions.R           # Stacked bar/proportion plots
│       └── plot_utils.R                    # Shared plotting utilities
│
├── scripts/                            # Shell scripts (unchanged)
│   ├── runFeatureCounts.sh
│   └── runFeatureCounts_TE_and_genes.sh
│
├── examples/                           # Template pipeline scripts
│   ├── 1.5.te_processing_combined.R    # Phase 1: Combined gene+TE
│   ├── 1.5.te_processing_separate.R    # Phase 1: Separate TE-only
│   ├── 1.6.te_master_tables.R          # Phase 2: Master table generation
│   └── 2.5.te_visualization.R          # Phase 3: TE visualizations
│
├── inst/
│   └── config/
│       └── te_config_template.yaml     # Config template for projects
│
└── tests/
    ├── test_te_enrichment.R
    ├── test_te_parsing.R
    └── testdata/
```

---

## 3. Dual-Mode Architecture: Combined vs Separate Analysis

### 3.1 The Scientific Reality

**Why This Matters:**

Different preprocessing pipelines produce different types of count matrices:

| Aspect | Combined Analysis | Separate Analysis |
|--------|-------------------|-------------------|
| Input | Single `combined_gene_TE_counts.tsv` | Separate `gene_counts.txt` + `te_counts.txt` |
| When Valid | Non-overlapping annotations, consistent library prep | Overlapping annotations, or methodological separation desired |
| Normalization | Single TMM on combined matrix | Separate (or borrowed library sizes) |
| Hypothesis Universe | Genes + TEs tested together | Genes and TEs as separate universes |

### 3.2 Factory Function Architecture

Create two entry points that produce structurally identical outputs:

#### `create_combined_dge.R` - The Gold Standard

```r
#' Create Combined Gene+TE DGEList
#'
#' For projects where gene and TE annotations are non-overlapping and
#' combined analysis is scientifically valid.
#'
#' @param gene_counts Gene count matrix (genes x samples)
#' @param te_counts TE count matrix (TEs x samples)
#' @param samples Sample metadata data.frame
#' @param te_annotation TE annotation from build_te_annotation() or NULL (auto-parse)
#' @param validate Check for ID collisions and sample consistency (recommended)
#' @return DGEList with feature_type column ("gene" or "te")
#'
#' @details
#' This function:
#' 1. Validates no duplicate IDs between gene and TE matrices
#' 2. Validates identical sample order
#' 3. Row-binds the matrices
#' 4. Adds feature_type metadata
#' 5. Parses TE annotation (subfamily, family, class)
#'
#' Use this when:
#' - Your TE annotation excludes exonic TEs (no overlap with genes)
#' - You want genes and TEs in a single statistical framework
#' - You want TMM normalization on the combined library
create_combined_dge <- function(
    gene_counts,
    te_counts,
    samples,
    te_annotation = NULL,
    validate = TRUE
) {
    # Validation
    if (validate) {
        validate_no_id_collision(gene_counts, te_counts)
        validate_sample_consistency(gene_counts, te_counts)
        validate_library_size_ratio(gene_counts, te_counts)
    }

    # Combine matrices
    combined <- rbind(gene_counts, te_counts)

    # Create DGEList
    dge <- DGEList(counts = combined, samples = samples)

    # Add feature metadata
    dge$genes <- data.frame(
        feature_id = rownames(combined),
        feature_type = c(
            rep("gene", nrow(gene_counts)),
            rep("te", nrow(te_counts))
        ),
        row.names = rownames(combined)
    )

    # Add TE annotation
    te_ids <- rownames(te_counts)
    if (is.null(te_annotation)) {
        te_annotation <- build_te_annotation(te_ids)
    }

    # Merge TE annotation (genes will have NA for TE-specific columns)
    dge$genes <- merge_te_annotation(dge$genes, te_annotation)

    # Add metadata about analysis mode
    dge$analysis_mode <- "combined"

    return(dge)
}
```

#### `create_separate_dge.R` - For Separate Analysis

```r
#' Create TE-Only DGEList for Separate Analysis
#'
#' For projects where separate gene/TE analysis is required.
#'
#' @param te_counts TE count matrix (TEs x samples)
#' @param samples Sample metadata data.frame
#' @param te_annotation TE annotation or NULL (auto-parse)
#' @param gene_dge_reference Optional: Gene DGEList for library size borrowing
#' @return DGEList with TE features only
#'
#' @details
#' If gene_dge_reference is provided, library sizes are borrowed from the
#' gene matrix. This is recommended to prevent scaling artifacts from
#' global TE derepression.
#'
#' Use this when:
#' - Annotations may overlap (conservative approach)
#' - You want genes and TEs as separate hypothesis universes
#' - Different preprocessing was used for genes vs TEs
create_separate_te_dge <- function(
    te_counts,
    samples,
    te_annotation = NULL,
    gene_dge_reference = NULL
) {
    # Create DGEList
    dge <- DGEList(counts = te_counts, samples = samples)

    # Add feature metadata
    dge$genes <- data.frame(
        feature_id = rownames(te_counts),
        feature_type = "te",
        row.names = rownames(te_counts)
    )

    # Add TE annotation
    if (is.null(te_annotation)) {
        te_annotation <- build_te_annotation(rownames(te_counts))
    }
    dge$genes <- cbind(dge$genes, te_annotation[, -1])  # Exclude Symbol column

    # Library size handling
    if (!is.null(gene_dge_reference)) {
        message("[INFO] Borrowing library sizes from gene reference for normalization")
        # Use gene library sizes to normalize TEs
        # This ensures TE expression changes are measured relative to total biological library
        common_samples <- intersect(colnames(dge), colnames(gene_dge_reference))
        if (length(common_samples) < ncol(dge)) {
            warning("Not all TE samples found in gene reference")
        }
        dge$samples$lib.size[match(common_samples, colnames(dge))] <-
            gene_dge_reference$samples$lib.size[match(common_samples, colnames(gene_dge_reference))]
    } else {
        warning("[CAUTION] Normalizing TEs without gene reference. ",
                "Global TE derepression may skew results.")
    }

    # Add metadata about analysis mode
    dge$analysis_mode <- "separate"

    return(dge)
}
```

### 3.3 Validation Functions

#### `validate_te_input.R`

```r
#' Validate No ID Collision Between Gene and TE Matrices
#'
#' @param gene_counts Gene count matrix
#' @param te_counts TE count matrix
#' @return TRUE if valid, error if collision detected
validate_no_id_collision <- function(gene_counts, te_counts) {
    gene_ids <- rownames(gene_counts)
    te_ids <- rownames(te_counts)

    collision <- intersect(gene_ids, te_ids)

    if (length(collision) > 0) {
        stop(
            "ID collision detected between gene and TE matrices!\n",
            "Found ", length(collision), " duplicate IDs:\n",
            paste(head(collision, 10), collapse = ", "),
            if (length(collision) > 10) paste0(" ... and ", length(collision) - 10, " more"),
            "\n\nThis usually means TE annotation overlaps with gene annotation.\n",
            "Solutions:\n",
            "  1. Use separate analysis mode: create_separate_te_dge()\n",
            "  2. Filter overlapping TEs from your SAF annotation\n",
            "  3. Review preprocessing pipeline"
        )
    }

    message("[OK] No ID collision between gene and TE matrices")
    return(TRUE)
}

#' Validate Sample Consistency Between Matrices
validate_sample_consistency <- function(gene_counts, te_counts) {
    gene_samples <- colnames(gene_counts)
    te_samples <- colnames(te_counts)

    if (!identical(gene_samples, te_samples)) {
        missing_in_te <- setdiff(gene_samples, te_samples)
        missing_in_gene <- setdiff(te_samples, gene_samples)

        stop(
            "Sample mismatch between gene and TE matrices!\n",
            if (length(missing_in_te) > 0) paste0("Missing in TE: ", paste(missing_in_te, collapse = ", "), "\n"),
            if (length(missing_in_gene) > 0) paste0("Missing in Gene: ", paste(missing_in_gene, collapse = ", "), "\n"),
            "\nEnsure both matrices have identical samples in the same order."
        )
    }

    message("[OK] Sample consistency validated")
    return(TRUE)
}

#' Validate Library Size Ratio (Warning Only)
validate_library_size_ratio <- function(gene_counts, te_counts, warn_threshold = 0.3) {
    gene_lib_sizes <- colSums(gene_counts)
    te_lib_sizes <- colSums(te_counts)
    te_proportion <- te_lib_sizes / (gene_lib_sizes + te_lib_sizes)

    if (any(te_proportion > warn_threshold)) {
        warning(
            "[CAUTION] High TE proportion detected in some samples:\n",
            "Max TE proportion: ", round(max(te_proportion) * 100, 1), "%\n",
            "Samples above ", warn_threshold * 100, "% threshold: ",
            sum(te_proportion > warn_threshold), "\n",
            "This may indicate preprocessing issues or extreme TE derepression.\n",
            "Review your annotation and counting pipeline."
        )
    } else {
        message("[OK] TE proportions within normal range (max: ",
                round(max(te_proportion) * 100, 1), "%)")
    }

    return(invisible(te_proportion))
}
```

### 3.4 When to Use Which Mode

| Scenario | Recommended Mode | Rationale |
|----------|-----------------|-----------|
| Exonic TEs removed from SAF | Combined | No overlap, unified statistical framework |
| Used TEtranscripts pipeline | Combined | Produces intrinsically consistent matrices |
| Annotations may overlap | Separate | Conservative; avoids double-counting |
| Want separate FDR thresholds | Separate | Different hypothesis universes |
| Different library prep for genes/TEs | Separate | Cannot assume comparable library sizes |
| Exploratory cross-analysis | Combined | Allows TE-gene correlation analysis |

### 3.5 Normalization Philosophy

**For Combined Mode:**
- Use TMM on the combined DGEList
- Genes (~78k features) dominate the normalization
- This is a **feature, not a bug**: TE changes are measured relative to the stable gene background

**For Separate Mode:**
- **Recommended**: Borrow library sizes from the gene matrix
- **Alternative**: Independent TMM (with warning about potential artifacts)
- This prevents global TE derepression from skewing the normalization

```r
# Combined mode normalization
dge_combined <- calcNormFactors(dge_combined, method = "TMM")

# Separate mode with library size borrowing
dge_te <- create_separate_te_dge(
    te_counts,
    samples,
    gene_dge_reference = dge_genes  # Borrow library sizes
)
dge_te <- calcNormFactors(dge_te, method = "TMM")
```

---

## 4. The Rule of Three: Function Separation

Each analysis capability is split into **three distinct functions** across phases:

### 4.1 TE Enrichment Example

#### Phase 1: Compute (`compute_te_enrichment.R`)

```r
#' Compute TE Family/Subfamily Enrichment Statistics
#'
#' @param fit Limma fit object with t-statistics
#' @param dge DGEList with TE annotation in $genes slot
#' @param te_level "subfamily", "family", or "class"
#' @param contrasts Character vector of contrast names (NULL = all)
#' @param p_adjust_method P-value adjustment method
#' @param adjust_across "all" or "by_contrast"
#' @param exclude_groups Groups to exclude (e.g., c("Unspecified", "Unknown"))
#' @return List with raw enrichment statistics (NO visualization)
compute_te_enrichment <- function(
    fit,
    dge,
    te_level = c("subfamily", "family", "class"),
    contrasts = NULL,
    p_adjust_method = "BH",
    adjust_across = c("all", "by_contrast"),
    exclude_groups = NULL
) {
    # Returns: list(
    #   results = tibble(te_group, contrast, pval_raw, t_mean, avg_logFC, n_elements),
    #   metadata = list(te_level, adjust_method, analysis_mode, ...)
    # )
}
```

#### Phase 2: Format (`format_te_enrichment.R`)

```r
#' Format TE Enrichment Results to Master Table Schema
#'
#' @param enrichment_results Output from compute_te_enrichment()
#' @param alpha Significance threshold for direction annotation
#' @return Tibble conforming to master_te_enrichment.csv schema
format_te_enrichment <- function(
    enrichment_results,
    alpha = 0.05
) {
    # Returns: tibble with standardized schema
}
```

#### Phase 3: Plot (`plot_te_enrichment_dotplot.R`)

```r
#' Create TE Enrichment Dotplot
#'
#' @param enrichment_table Master table from format_te_enrichment() or CSV
#' @param contrast_groups Named list for faceting (e.g., list(GroupA = c("c1","c2")))
#' @param contrast_labels Named vector for display names
#' @param pval_threshold Threshold for significance outline
#' @param color_scheme Color palette name
#' @return ggplot object (NO computation)
plot_te_enrichment_dotplot <- function(
    enrichment_table,
    contrast_groups = NULL,
    contrast_labels = NULL,
    pval_threshold = 0.05,
    color_scheme = "plasma",
    te_levels_to_show = NULL
) {
    # Pure visualization - reads pre-computed data
}
```

### 4.2 Function Responsibility Summary

| Function Type | Phase | Input | Output | Side Effects |
|--------------|-------|-------|--------|--------------|
| `create_*_dge` | 0 | Count matrices | DGEList | None |
| `compute_*` | 1 | DGE/fit objects | Raw results list | None (cacheable) |
| `format_*` | 2 | Raw results | Master table tibble | Writes CSV |
| `plot_*` | 3 | Master table | ggplot object | None |

---

## 5. Master Table Schemas

### 5.1 `master_te_enrichment.csv`

| Column | Type | Description |
|--------|------|-------------|
| `te_group` | character | TE subfamily/family/class name |
| `te_level` | character | "subfamily", "family", or "class" |
| `te_class` | character | TE class (LINE, SINE, LTR, DNA, RC) |
| `replication_type` | character | "cytoplasmic" or "nuclear" |
| `contrast` | character | Contrast name |
| `pvalue` | numeric | Raw p-value from geneSetTest |
| `padj` | numeric | Adjusted p-value |
| `enrichment_score` | numeric | Mean t-statistic (for magnitude visualization) |
| `avg_logfc` | numeric | Mean logFC of elements in group |
| `direction` | character | "up", "down", or "ns" |
| `n_elements` | integer | Number of TEs in group |
| `prop_significant` | numeric | Proportion individually significant (optional) |
| `analysis_mode` | character | "combined" or "separate" |

### 5.2 `master_te_de.csv`

| Column | Type | Description |
|--------|------|-------------|
| `te_id` | character | Full TE identifier (e.g., "HAL1:L1:LINE") |
| `subfamily` | character | TE subfamily |
| `family` | character | TE family |
| `te_class` | character | TE class |
| `contrast` | character | Contrast name |
| `logfc` | numeric | Log fold change |
| `pvalue` | numeric | Raw p-value |
| `padj` | numeric | Adjusted p-value |
| `t_stat` | numeric | t-statistic |
| `ave_expr` | numeric | Average expression |
| `analysis_mode` | character | "combined" or "separate" |

### 5.3 `master_te_aggregates.csv`

| Column | Type | Description |
|--------|------|-------------|
| `te_group` | character | TE subfamily or family |
| `te_level` | character | "subfamily" or "family" |
| `sample_id` | character | Sample ID |
| `expression` | numeric | Aggregated expression value |
| `expression_type` | character | "logCPM", "zscore", or "cpm" |

---

## 6. Checkpoint Strategy

### 6.1 Naming Convention

Checkpoint files differentiate mode via **filename**, not structure:

```
03_results/checkpoints/
├── # Main pipeline (Phase 1.1)
├── 1.1_dge_normalized.rds              # May be combined or genes-only
├── 1.1_fit_object.rds
├── 1.1_de_results.rds
│
├── # TE-specific checkpoints (Phase 1.5)
├── 1.5_te_dge_combined.rds             # Combined mode DGE
├── 1.5_te_dge_separate.rds             # Separate mode DGE
├── 1.5_te_fit_combined.rds             # Fit from combined analysis
├── 1.5_te_fit_separate.rds             # Fit from separate analysis
├── 1.5_te_enrichment_subfamily.rds     # All contrasts, subfamily level
├── 1.5_te_enrichment_family.rds        # All contrasts, family level
├── 1.5_te_aggregates.rds               # Aggregated expression
└── 1.5_te_correlations.rds             # TE-pathway correlations
```

### 6.2 Checkpoint Structure

Each checkpoint includes provenance metadata:

```r
# Structure of saved checkpoint
result_object <- list(
    data = computed_results,
    metadata = list(
        analysis_mode = "combined",  # or "separate"
        te_level = "subfamily",
        p_adjust_method = "BH",
        excluded_groups = c("Unspecified", "Unknown"),
        created_at = Sys.time(),
        toolkit_version = "1.0.0"
    )
)
saveRDS(result_object, checkpoint_path)
```

---

## 7. Parameterization Strategy

### 7.1 Remove All Hardcoded Project Logic

**Before (BAD)**:
```r
celltype = case_when(
    str_starts(contrast, "Th1_")  ~ "Th1",
    str_starts(contrast, "Th17_") ~ "Th17",
    TRUE                          ~ "Other"
)
```

**After (GOOD)**:
```r
# User provides grouping via config/parameter
plot_te_enrichment_dotplot <- function(
    enrichment_table,
    contrast_groups = NULL,  # User-provided grouping
    ...
) {
    if (!is.null(contrast_groups)) {
        enrichment_table <- enrichment_table %>%
            mutate(group = map_contrast_to_group(contrast, contrast_groups))
    }
}

# Usage:
plot_te_enrichment_dotplot(
    enrichment_table,
    contrast_groups = list(
        Th1 = c("Th1_0h_vs_ctrl", "Th1_24h_vs_ctrl"),
        Th17 = c("Th17_0h_vs_ctrl", "Th17_24h_vs_ctrl")
    )
)
```

### 7.2 Project Configuration Template

```yaml
# inst/config/te_config_template.yaml

project:
  id: "PROJECT-ID"
  species: "Mus musculus"

te_analysis:
  # Analysis mode: "combined" or "separate"
  mode: "combined"

  # TE grouping levels to analyze
  levels: ["subfamily", "family"]

  # Groups to exclude from enrichment
  exclude_groups: ["Unspecified", "Unknown"]

  # Replication location classification
  cytoplasmic_classes: ["LINE", "SINE", "LTR", "Retroposon"]
  nuclear_classes: ["DNA", "RC"]

  # Statistical parameters
  p_adjust_method: "BH"
  adjust_across: "all"  # or "by_contrast"

visualization:
  # Contrast grouping for faceted plots
  contrast_groups:
    Group1: ["contrast_A", "contrast_B"]
    Group2: ["contrast_C", "contrast_D"]

  # Display labels
  contrast_labels:
    contrast_A: "0h"
    contrast_B: "24h"

  # Direction labels (for single contrast plots)
  direction_labels:
    positive: "Higher in Treatment"
    negative: "Higher in Control"

  # Color schemes
  colors:
    enrichment_palette: "plasma"
    heatmap_palette: "viridis"

  # Plot dimensions
  plot_width: 10
  plot_height: 8
  plot_dpi: 300
```

---

## 8. Implementation Roadmap

### Phase 1: Data Ingestion (Priority: CRITICAL)

| Task | Files | Dependencies |
|------|-------|--------------|
| Create factory functions | `create_combined_dge.R`, `create_separate_dge.R` | None |
| Create validation functions | `validate_te_input.R` | None |
| Update `te_utils.R` | `te_utils.R` | None |

### Phase 2: Core Processing Functions (Priority: HIGH)

| Task | Files | Dependencies |
|------|-------|--------------|
| Create `compute_te_enrichment.R` | `R/compute_te_enrichment.R` | Factory functions |
| Create `compute_te_aggregates.R` | `R/compute_te_aggregates.R` | Factory functions |
| Create `compute_te_de.R` | `R/compute_te_de.R` | Factory functions |

### Phase 3: Master Table Functions (Priority: HIGH)

| Task | Files | Dependencies |
|------|-------|--------------|
| Create `format_te_enrichment.R` | `R/format_te_enrichment.R` | compute_te_enrichment |
| Create `format_te_de.R` | `R/format_te_de.R` | None |
| Create `format_te_aggregates.R` | `R/format_te_aggregates.R` | compute_te_aggregates |

### Phase 4: Visualization Functions (Priority: MEDIUM)

| Task | Files | Dependencies |
|------|-------|--------------|
| Create `plot_te_enrichment_dotplot.R` | `R/plot_te_enrichment_dotplot.R` | format_te_enrichment |
| Create `plot_te_heatmap.R` | `R/plot_te_heatmap.R` | format_te_aggregates |
| Create `plot_te_enrichment_single.R` | `R/plot_te_enrichment_single.R` | format_te_enrichment |
| Create `plot_utils.R` | `R/plot_utils.R` | None |

### Phase 5: Extended Features (Priority: LOW)

| Task | Files | Dependencies |
|------|-------|--------------|
| Create `compute_te_correlations.R` | `R/compute_te_correlations.R` | compute_te_aggregates |
| Create `plot_te_correlation.R` | `R/plot_te_correlation.R` | compute_te_correlations |
| Create `plot_te_volcano.R` | `R/plot_te_volcano.R` | format_te_de |
| Create `plot_te_proportions.R` | `R/plot_te_proportions.R` | compute_te_proportions |

### Phase 6: Examples & Documentation (Priority: MEDIUM)

| Task | Files | Dependencies |
|------|-------|--------------|
| Create combined mode example | `examples/1.5.te_processing_combined.R` | All compute functions |
| Create separate mode example | `examples/1.5.te_processing_separate.R` | All compute functions |
| Create visualization example | `examples/2.5.te_visualization.R` | All plot functions |
| Update README.md | `README.md` | All |

---

## 9. Script Migration Map

| Old Script | New Location(s) | Notes |
|------------|-----------------|-------|
| `te_utils.R` | `R/te_utils.R` | Keep + enhance with feature_type |
| `te_enrichment_multi_contrast.R` | Split into: | |
| - Computation | `R/compute_te_enrichment.R` | Pure stats |
| - Visualization | `R/plot_te_enrichment_dotplot.R` | Pure viz |
| `te_enrichment_single_contrast.R` | Split into: | |
| - Computation | `R/compute_te_enrichment.R` | Reuse same function |
| - Visualization | `R/plot_te_enrichment_single.R` | Different viz style |
| `create_te_heatmap.R` | Split into: | |
| - Aggregation | `R/compute_te_aggregates.R` | Pre-compute |
| - Visualization | `R/plot_te_heatmap.R` | Read aggregates |
| `calculate_te_pathway_correlations.R` | `R/compute_te_correlations.R` | Add caching |
| `visualise_te_pathway_correlations.R` | `R/plot_te_correlation.R` | Keep logic |
| `create_all_te_pathway_correlations.R` | **DELETE** | Move to examples |
| `te_example_workflow.R` | **REPLACE** with `examples/*.R` | Phase-based |

---

## Appendix A: Integration with Main Pipeline

```
Main RNAseq Pipeline
│
├── Phase 1.1: Core Analysis
│   ├── Load counts matrix
│   │   ├── Option A: Combined gene+TE matrix → use directly
│   │   └── Option B: Separate matrices → decide analysis mode
│   ├── Create DGEList (combined or genes-only)
│   ├── TMM normalization
│   ├── DE analysis (limma-voom)
│   └── Checkpoints: 1.1_dge_normalized.rds, 1.1_fit_object.rds
│
├── Phase 1.3: GSEA (regular genes)
│   └── Checkpoint: 1.3_gsea_results.rds
│
├── Phase 1.5: TE Analysis ← TE-RNASEQ-TOOLKIT INTEGRATION
│   ├── Choose analysis mode (combined vs separate)
│   ├── compute_te_enrichment()
│   ├── compute_te_aggregates()
│   └── Checkpoints: 1.5_te_*.rds
│
├── Phase 2.x: Master Tables
│   ├── master_de_table.csv (genes)
│   ├── master_gsea_table.csv
│   ├── master_te_enrichment.csv ← NEW
│   ├── master_te_de.csv ← NEW
│   └── master_te_aggregates.csv ← NEW
│
└── Phase 3.x: Visualization
    ├── Gene visualizations (2.1-2.4)
    └── TE visualizations (2.5) ← NEW
```

---

## Appendix B: Example Pipeline Script (Combined Mode)

```r
#!/usr/bin/env Rscript
# 1.5.te_analysis_combined.R - TE-specific analysis (combined mode)
# Project: {PROJECT_ID}
# Phase: 1.5
# Dependencies: gene_counts.txt, te_counts.txt, or combined_gene_TE_counts.tsv
# Outputs: 1.5_te_*.rds, master_te_*.csv

# ============================================================================
# SETUP
# ============================================================================
message("=================================================================")
message("Phase 1.5: TE Analysis (Combined Mode)")
message("=================================================================")

source("02_analysis/config/config.R")
source_toolkit("TE-RNAseq-toolkit")

# ============================================================================
# LOAD AND VALIDATE INPUT
# ============================================================================

# Option 1: Load pre-combined matrix
combined_counts <- read.delim(
    file.path(DIR_DATA, "combined_gene_TE_counts.tsv"),
    row.names = 1
)

# Parse which rows are TEs vs genes
te_annotation <- build_te_annotation(rownames(combined_counts))
is_te <- !is.na(te_annotation$subfamily)

gene_counts <- combined_counts[!is_te, ]
te_counts <- combined_counts[is_te, ]

# ============================================================================
# CREATE COMBINED DGEList
# ============================================================================
dge <- load_or_compute(
    checkpoint_file = "1.5_te_dge_combined.rds",
    description = "Combined gene+TE DGEList",
    compute_fn = function() {
        dge <- create_combined_dge(
            gene_counts = gene_counts,
            te_counts = te_counts,
            samples = metadata,
            validate = TRUE
        )
        calcNormFactors(dge, method = "TMM")
    }
)

# ============================================================================
# FIT MODEL (if not using main pipeline fit)
# ============================================================================
fit <- load_or_compute(
    checkpoint_file = "1.5_te_fit_combined.rds",
    description = "limma-voom fit on combined data",
    compute_fn = function() {
        design <- model.matrix(~ 0 + group, data = dge$samples)
        v <- voom(dge, design)
        fit <- lmFit(v, design)
        fit <- contrasts.fit(fit, contrast_matrix)
        eBayes(fit)
    }
)

# ============================================================================
# COMPUTE TE ENRICHMENT
# ============================================================================
te_enrichment_sub <- load_or_compute(
    checkpoint_file = "1.5_te_enrichment_subfamily.rds",
    description = "TE subfamily enrichment",
    compute_fn = function() {
        compute_te_enrichment(
            fit = fit,
            dge = dge,
            te_level = "subfamily",
            exclude_groups = CONFIG$te_analysis$exclude_groups
        )
    }
)

te_enrichment_fam <- load_or_compute(
    checkpoint_file = "1.5_te_enrichment_family.rds",
    description = "TE family enrichment",
    compute_fn = function() {
        compute_te_enrichment(
            fit = fit,
            dge = dge,
            te_level = "family",
            exclude_groups = CONFIG$te_analysis$exclude_groups
        )
    }
)

# ============================================================================
# CREATE MASTER TABLES
# ============================================================================
message("[EXPORT] Creating master tables...")

te_master_sub <- format_te_enrichment(te_enrichment_sub)
te_master_fam <- format_te_enrichment(te_enrichment_fam)
te_master_all <- bind_rows(te_master_sub, te_master_fam)
write_csv(te_master_all, file.path(DIR_TABLES, "master_te_enrichment.csv"))

message("[DONE] Phase 1.5 complete")
```

---

## Appendix C: Example Pipeline Script (Separate Mode)

```r
#!/usr/bin/env Rscript
# 1.5.te_analysis_separate.R - TE-specific analysis (separate mode)
# Project: {PROJECT_ID}
# Phase: 1.5
# Dependencies: 1.1_dge_normalized.rds (genes), te_counts.txt
# Outputs: 1.5_te_*.rds, master_te_*.csv

# ============================================================================
# SETUP
# ============================================================================
message("=================================================================")
message("Phase 1.5: TE Analysis (Separate Mode)")
message("=================================================================")

source("02_analysis/config/config.R")
source_toolkit("TE-RNAseq-toolkit")

# ============================================================================
# LOAD DEPENDENCIES
# ============================================================================

# Load gene DGE from main pipeline (for library size reference)
dge_genes <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_dge_normalized.rds"))

# Load TE counts
te_counts <- read.delim(
    file.path(DIR_DATA, "te_counts_matrix.txt"),
    row.names = 1
)

# ============================================================================
# CREATE SEPARATE TE DGEList (with library size borrowing)
# ============================================================================
dge_te <- load_or_compute(
    checkpoint_file = "1.5_te_dge_separate.rds",
    description = "Separate TE DGEList",
    compute_fn = function() {
        dge <- create_separate_te_dge(
            te_counts = te_counts,
            samples = metadata,
            gene_dge_reference = dge_genes  # Borrow library sizes
        )
        calcNormFactors(dge, method = "TMM")
    }
)

# ============================================================================
# FIT MODEL ON TE DATA
# ============================================================================
fit_te <- load_or_compute(
    checkpoint_file = "1.5_te_fit_separate.rds",
    description = "limma-voom fit on TE data",
    compute_fn = function() {
        design <- model.matrix(~ 0 + group, data = dge_te$samples)
        v <- voom(dge_te, design)
        fit <- lmFit(v, design)
        fit <- contrasts.fit(fit, contrast_matrix)
        eBayes(fit)
    }
)

# ... rest of analysis same as combined mode
```

---

**End of Restructuring Plan**
