#!/usr/bin/env Rscript
# =============================================================================
# 2.5.te_visualization.R - TE Visualization Pipeline
# =============================================================================
# TE-RNAseq-toolkit Example Script
# Phase: 2.5
#
# Description:
#   Template script for creating TE visualizations from master tables.
#   This script reads pre-computed results from Phase 1.5 and generates
#   publication-ready figures.
#
# Prerequisites:
#   - Phase 1.5 complete (checkpoints and master tables exist)
#   - master_te_enrichment.csv
#   - master_te_de.csv
#   - master_te_aggregates.csv
#
# Outputs:
#   - plots/te_enrichment_dotplot.pdf
#   - plots/te_enrichment_heatmap.pdf
#   - plots/te_volcano.pdf
#   - plots/te_proportions.pdf
#
# Usage:
#   1. Copy this file to your project's 02_analysis/scripts/ directory
#   2. Update paths for your project
#   3. Run: Rscript 2.5.te_visualization.R
#
# =============================================================================

message("=================================================================")
message("Phase 2.5: TE Visualizations")
message("=================================================================")
message("Started: ", Sys.time())

# =============================================================================
# SETUP
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(edgeR)
    library(tibble)
})

# -----------------------------------------------------------------------------
# PROJECT CONFIGURATION - UPDATE THESE PATHS
# -----------------------------------------------------------------------------

PROJECT_ROOT <- "/workspaces/YOUR_PROJECT"  # <-- UPDATE THIS

DIR_CHECKPOINTS <- file.path(PROJECT_ROOT, "03_results/checkpoints")
DIR_TABLES      <- file.path(PROJECT_ROOT, "03_results/tables")
DIR_PLOTS       <- file.path(PROJECT_ROOT, "03_results/plots")

TOOLKIT_PATH <- "/workspaces/AdaW_eWAT_WL_2025/01_modules/TE-RNAseq-toolkit/R"

# -----------------------------------------------------------------------------
# SOURCE TOOLKIT FUNCTIONS
# -----------------------------------------------------------------------------

message("[SETUP] Loading TE-RNAseq-toolkit visualization functions...")

source(file.path(TOOLKIT_PATH, "te_utils.R"))
source(file.path(TOOLKIT_PATH, "plot_utils.R"))
source(file.path(TOOLKIT_PATH, "plot_te_enrichment_dotplot.R"))
source(file.path(TOOLKIT_PATH, "plot_te_enrichment_single.R"))
source(file.path(TOOLKIT_PATH, "plot_te_heatmap.R"))
source(file.path(TOOLKIT_PATH, "plot_te_volcano.R"))
source(file.path(TOOLKIT_PATH, "plot_te_proportions.R"))
source(file.path(TOOLKIT_PATH, "plot_te_correlation.R"))
source(file.path(TOOLKIT_PATH, "compute_te_aggregates.R"))
source(file.path(TOOLKIT_PATH, "compute_te_proportions.R"))

message("[SETUP] Toolkit loaded successfully")

# Create output directory
dir.create(DIR_PLOTS, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD MASTER TABLES AND CHECKPOINTS
# =============================================================================

message("\n[STEP 1] Loading data...")

# NOTE: In a real project, these files will be loaded from your results directory.
# For demonstration, we'll create synthetic data if files don't exist.

# Check if master tables exist
if (file.exists(file.path(DIR_TABLES, "master_te_enrichment.csv"))) {
    master_te_enrichment <- read_csv(file.path(DIR_TABLES, "master_te_enrichment.csv"),
                                     show_col_types = FALSE)
    message("  Loaded master_te_enrichment.csv")
} else {
    # For demonstration, create synthetic master table
    message("  Master tables not found - creating synthetic data for demo...")

    # Create synthetic enrichment data
    te_groups <- c("L1Md_A", "L1Md_T", "IAP1", "IAP2", "MMERGLN", "MERVL", "B1_Mur1", "B1_Mur2")
    contrasts <- c("Treatment_vs_Control")
    te_levels <- c("subfamily")

    master_te_enrichment <- expand.grid(
        te_group = te_groups,
        contrast = contrasts,
        te_level = te_levels,
        stringsAsFactors = FALSE
    ) %>%
        mutate(
            te_class = case_when(
                grepl("^L1", te_group) ~ "LINE",
                grepl("^IAP|MMERGLN|MERVL", te_group) ~ "LTR",
                grepl("^B1", te_group) ~ "SINE",
                TRUE ~ "Other"
            ),
            replication_type = ifelse(te_class %in% c("LINE", "SINE", "LTR"), "cytoplasmic", "nuclear"),
            pvalue = runif(n()),
            padj = p.adjust(pvalue, "BH"),
            enrichment_score = rnorm(n(), mean = 0, sd = 1.5),
            avg_logfc = enrichment_score * 0.3 + rnorm(n(), 0, 0.2),
            direction = case_when(
                padj >= 0.05 ~ "ns",
                avg_logfc > 0 ~ "up",
                avg_logfc < 0 ~ "down"
            ),
            n_elements = sample(5:50, n(), replace = TRUE),
            prop_significant = runif(n(), 0, 0.3),
            analysis_mode = "combined"
        )
}

# Load TE DGE and aggregates for heatmap and proportions
if (file.exists(file.path(DIR_CHECKPOINTS, "1.5_te_dge_combined.rds"))) {
    dge <- readRDS(file.path(DIR_CHECKPOINTS, "1.5_te_dge_combined.rds"))
    message("  Loaded DGE from checkpoint")
} else {
    # Create synthetic DGE for demo
    message("  Creating synthetic DGE for demo...")

    set.seed(42)
    n_tes <- 20
    n_samples <- 6
    sample_names <- paste0("Sample", 1:n_samples)

    te_ids <- c(
        "L1Md_A:L1:LINE", "L1Md_T:L1:LINE", "L1Md_F:L1:LINE", "L1Md_Gf:L1:LINE",
        "IAP1:ERVK:LTR", "IAP2:ERVK:LTR", "MMERGLN:ERVL:LTR", "MERVL:ERVL:LTR",
        "MMERVK10C:ERVK:LTR", "ETnERV3:ERVK:LTR", "MuRRS:ERVK:LTR", "RLTR1B:ERV1:LTR",
        "B1_Mur1:Alu:SINE", "B1_Mur2:Alu:SINE", "B2_Mm1a:B2:SINE", "B2_Mm2:B2:SINE",
        "MER20:hAT-Charlie:DNA", "Charlie1:hAT-Charlie:DNA", "RSINE1:SINE:SINE", "ID_B1:Alu:SINE"
    )

    te_counts <- matrix(
        rpois(n_tes * n_samples, lambda = 50),
        nrow = n_tes,
        dimnames = list(te_ids, sample_names)
    )

    metadata <- data.frame(
        sample_id = sample_names,
        group = factor(rep(c("Control", "Treatment"), each = 3)),
        row.names = sample_names
    )

    te_annotation <- build_te_annotation(te_ids)

    dge <- edgeR::DGEList(counts = te_counts, samples = metadata)
    dge$genes <- data.frame(
        feature_id = te_ids,
        feature_type = "te",
        subfamily = te_annotation$subfamily,
        family = te_annotation$family,
        class = te_annotation$class,
        replication_type = te_annotation$replication_type,
        row.names = te_ids
    )
    dge$analysis_mode <- "combined"
    dge <- edgeR::calcNormFactors(dge)
}

# Load or compute aggregates
if (file.exists(file.path(DIR_CHECKPOINTS, "1.5_te_aggregates.rds"))) {
    te_aggregates <- readRDS(file.path(DIR_CHECKPOINTS, "1.5_te_aggregates.rds"))
    message("  Loaded aggregates from checkpoint")
} else {
    message("  Computing TE aggregates for demo...")
    te_aggregates <- compute_te_aggregates(
        dge = dge,
        te_level = "family",
        method = "mean",
        value_type = "logcpm"
    )
}

message("  Data loaded successfully")

# =============================================================================
# VISUALIZATION CONFIGURATION
# =============================================================================

# Contrast grouping (adjust for your project)
# Example: Group contrasts by cell type, treatment, or timepoint
contrast_groups <- list(
    Comparisons = c("Treatment_vs_Control")
    # Add more groups:
    # CellTypeA = c("A_Treatment_vs_Control"),
    # CellTypeB = c("B_Treatment_vs_Control")
)

# Contrast display labels (adjust for your project)
contrast_labels <- c(
    "Treatment_vs_Control" = "Treatment vs. Control"
    # Add more:
    # "A_Treatment_vs_Control" = "A: Treat/Ctrl"
)

# =============================================================================
# PLOT 1: TE ENRICHMENT DOTPLOT
# =============================================================================

message("\n[STEP 2] Creating TE enrichment dotplot...")

p_dotplot <- plot_te_enrichment_dotplot(
    enrichment_table = master_te_enrichment,
    contrast_groups = NULL,  # Use NULL for automatic grouping
    contrast_labels = contrast_labels,
    te_levels_to_show = c("subfamily", "family"),
    pval_threshold = 0.05,
    color_by = "padj",
    size_by = "enrichment_score",
    facet_by = c("te_class")
)

# Save plot
output_file <- file.path(DIR_PLOTS, "te_enrichment_dotplot.pdf")
ggsave(output_file, p_dotplot, width = 10, height = 8)
message("  Saved: ", output_file)

# =============================================================================
# PLOT 2: TE ENRICHMENT SINGLE CONTRAST
# =============================================================================

message("\n[STEP 3] Creating single-contrast enrichment plot...")

# For single contrast visualization
single_contrast <- unique(master_te_enrichment$contrast)[1]

p_single <- plot_te_enrichment_single(
    enrichment_table = master_te_enrichment,
    contrast_name = single_contrast,
    direction_labels = list(
        positive = "Up in Treatment",
        negative = "Down in Treatment"
    )
)

output_file <- file.path(DIR_PLOTS, "te_enrichment_single.pdf")
ggsave(output_file, p_single, width = 10, height = 6)
message("  Saved: ", output_file)

# =============================================================================
# PLOT 3: TE EXPRESSION HEATMAP
# =============================================================================

message("\n[STEP 4] Creating TE expression heatmap...")

# Heatmap of aggregated TE expression by family
p_heatmap <- plot_te_heatmap(
    aggregates = te_aggregates,
    samples = dge$samples,
    annotation_col = c("group"),
    annotation_row = c("class"),
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color_palette = "viridis",
    main_title = "TE Family Expression (z-score)"
)

# pheatmap saves to device, so we need to handle differently
output_file <- file.path(DIR_PLOTS, "te_expression_heatmap.pdf")
pdf(output_file, width = 8, height = 10)
plot_te_heatmap(
    aggregates = te_aggregates,
    samples = dge$samples,
    annotation_col = c("group"),
    annotation_row = c("class"),
    scale = "row",
    main_title = "TE Family Expression (z-score)"
)
dev.off()
message("  Saved: ", output_file)

# =============================================================================
# PLOT 4: TE PROPORTIONS
# =============================================================================

message("\n[STEP 5] Creating TE proportion plots...")

# Compute TE proportions if DGE available
te_props <- compute_te_proportions(dge, group_by = "group")

# Total TE proportion by sample
p_total_prop <- plot_te_total_prop(te_props, group_by = "group")
output_file <- file.path(DIR_PLOTS, "te_total_proportion.pdf")
ggsave(output_file, p_total_prop, width = 8, height = 5)
message("  Saved: ", output_file)

# TE composition by class
p_composition <- plot_te_composition(te_props, level = "class", group_by = "group")
output_file <- file.path(DIR_PLOTS, "te_composition_class.pdf")
ggsave(output_file, p_composition, width = 8, height = 6)
message("  Saved: ", output_file)

# =============================================================================
# PLOT 5: TE VOLCANO PLOT (if DE results available)
# =============================================================================

message("\n[STEP 6] Creating TE volcano plot...")

# Load or create DE results
if (file.exists(file.path(DIR_TABLES, "master_te_de.csv"))) {
    master_te_de <- read_csv(file.path(DIR_TABLES, "master_te_de.csv"), show_col_types = FALSE)
} else {
    # Create synthetic DE results for demo
    master_te_de <- data.frame(
        te_id = rownames(dge),
        contrast = single_contrast,
        logfc = rnorm(nrow(dge)),
        pvalue = runif(nrow(dge)),
        stringsAsFactors = FALSE
    ) %>%
        mutate(
            padj = p.adjust(pvalue, "BH"),
            feature_type = "te"
        )
}

# Volcano plot
p_volcano <- plot_te_volcano(
    de_results = master_te_de,
    highlight_te = TRUE,
    padj_threshold = 0.05,
    logfc_threshold = 1
)

output_file <- file.path(DIR_PLOTS, "te_volcano.pdf")
ggsave(output_file, p_volcano, width = 8, height = 7)
message("  Saved: ", output_file)

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=================================================================")
message("Phase 2.5 Complete: TE Visualizations")
message("=================================================================")
message("Completed: ", Sys.time())
message("")
message("Plots saved to: ", DIR_PLOTS)
message("")
message("Generated plots:")
message("  1. te_enrichment_dotplot.pdf - Multi-contrast enrichment overview")
message("  2. te_enrichment_single.pdf  - Single contrast detail view")
message("  3. te_expression_heatmap.pdf - Expression heatmap by family")
message("  4. te_total_proportion.pdf   - TE proportion by sample")
message("  5. te_composition_class.pdf  - TE class composition")
message("  6. te_volcano.pdf            - Volcano plot with TE highlighting")
message("")
message("Customize by editing:")
message("  - contrast_groups: Group contrasts for faceting")
message("  - contrast_labels: Human-readable contrast names")
message("  - color_palette: viridis, magma, rdbu")
