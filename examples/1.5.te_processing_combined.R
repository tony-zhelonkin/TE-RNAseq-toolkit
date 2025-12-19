#!/usr/bin/env Rscript
# =============================================================================
# 1.5.te_processing_combined.R - TE Analysis Pipeline (Combined Mode)
# =============================================================================
# TE-RNAseq-toolkit Example Script
# Phase: 1.5
#
# Description:
#   Template script for TE-specific analysis using COMBINED mode.
#   In combined mode, gene and TE counts are analyzed together in a single
#   DGEList with unified TMM normalization.
#
# Prerequisites:
#   - Gene count matrix (genes x samples)
#   - TE count matrix (TEs x samples)
#   - Sample metadata
#   - Non-overlapping gene/TE annotations (exonic TEs removed from SAF)
#
# Outputs:
#   - Checkpoints: 1.5_te_dge_combined.rds, 1.5_te_fit.rds, 1.5_te_enrichment_*.rds
#   - Master tables: master_te_enrichment.csv, master_te_de.csv
#
# Usage:
#   1. Copy this file to your project's 02_analysis/scripts/ directory
#   2. Update paths and parameters for your project
#   3. Run: Rscript 1.5.te_processing_combined.R
#
# =============================================================================

message("=================================================================")
message("Phase 1.5: TE Analysis (Combined Mode)")
message("=================================================================")
message("Started: ", Sys.time())

# =============================================================================
# SETUP
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(dplyr)
    library(tibble)
    library(readr)
})

# -----------------------------------------------------------------------------
# PROJECT CONFIGURATION - UPDATE THESE PATHS FOR YOUR PROJECT
# -----------------------------------------------------------------------------

# Option 1: Set paths directly
PROJECT_ROOT <- "/workspaces/YOUR_PROJECT"  # <-- UPDATE THIS

# Data directories
DIR_DATA        <- file.path(PROJECT_ROOT, "00_data/processed")
DIR_CHECKPOINTS <- file.path(PROJECT_ROOT, "03_results/checkpoints")
DIR_TABLES      <- file.path(PROJECT_ROOT, "03_results/tables")

# Input files
GENE_COUNTS_FILE <- file.path(DIR_DATA, "gene_counts.txt")       # <-- UPDATE
TE_COUNTS_FILE   <- file.path(DIR_DATA, "te_counts.txt")         # <-- UPDATE
METADATA_FILE    <- file.path(DIR_DATA, "sample_metadata.csv")   # <-- UPDATE

# TE-RNAseq-toolkit location
TOOLKIT_PATH <- "/workspaces/AdaW_eWAT_WL_2025/01_modules/TE-RNAseq-toolkit/R"

# Option 2: Or load from project config (recommended)
# CONFIG_FILE <- file.path(PROJECT_ROOT, "02_analysis/config/analysis_config.yaml")
# config <- yaml::read_yaml(CONFIG_FILE)
# DIR_DATA <- config$paths$data_dir
# etc.

# -----------------------------------------------------------------------------
# SOURCE TOOLKIT FUNCTIONS
# -----------------------------------------------------------------------------

message("[SETUP] Loading TE-RNAseq-toolkit functions...")

source(file.path(TOOLKIT_PATH, "te_utils.R"))
source(file.path(TOOLKIT_PATH, "validate_te_input.R"))
source(file.path(TOOLKIT_PATH, "create_combined_dge.R"))
source(file.path(TOOLKIT_PATH, "compute_te_enrichment.R"))
source(file.path(TOOLKIT_PATH, "compute_te_aggregates.R"))
source(file.path(TOOLKIT_PATH, "compute_te_de.R"))
source(file.path(TOOLKIT_PATH, "format_te_enrichment.R"))
source(file.path(TOOLKIT_PATH, "format_te_de.R"))
source(file.path(TOOLKIT_PATH, "format_te_aggregates.R"))

message("[SETUP] Toolkit loaded successfully")

# -----------------------------------------------------------------------------
# HELPER: load_or_compute pattern
# -----------------------------------------------------------------------------

#' Load from checkpoint or compute and save
#'
#' @param checkpoint_file Path to checkpoint RDS file
#' @param description Description for logging
#' @param compute_fn Function that computes the result
#' @param force_recompute If TRUE, always recompute
load_or_compute <- function(checkpoint_file, description, compute_fn, force_recompute = FALSE) {

    checkpoint_path <- file.path(DIR_CHECKPOINTS, checkpoint_file)

    if (file.exists(checkpoint_path) && !force_recompute) {
        message("[CACHE] Loading ", description, " from: ", checkpoint_file)
        return(readRDS(checkpoint_path))
    }

    message("[COMPUTE] Computing ", description, "...")
    result <- compute_fn()

    # Ensure directory exists
    dir.create(dirname(checkpoint_path), recursive = TRUE, showWarnings = FALSE)

    message("[SAVE] Saving to: ", checkpoint_file)
    saveRDS(result, checkpoint_path)

    return(result)
}

# =============================================================================
# LOAD INPUT DATA
# =============================================================================

message("\n[STEP 1] Loading input data...")

# NOTE: In a real project, replace this section with your data loading code
# This example shows the expected format

# Load sample metadata
# metadata <- read.csv(METADATA_FILE, row.names = 1)
# Or: metadata <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_sample_metadata.rds"))

# For demonstration, we'll create synthetic data
# REMOVE THIS SECTION and load your real data
set.seed(42)
n_genes <- 100
n_tes <- 20
n_samples <- 6

gene_ids <- paste0("Gene", 1:n_genes)
sample_names <- paste0("Sample", 1:n_samples)

# Synthetic gene counts
gene_counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes,
    dimnames = list(gene_ids, sample_names)
)

# Synthetic TE counts (TEtranscripts format: subfamily:family:class)
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

# Sample metadata
metadata <- data.frame(
    sample_id = sample_names,
    group = factor(rep(c("Control", "Treatment"), each = 3)),
    batch = factor(rep(c("A", "B"), 3)),
    row.names = sample_names
)

message("  Gene counts: ", nrow(gene_counts), " genes x ", ncol(gene_counts), " samples")
message("  TE counts: ", nrow(te_counts), " TEs x ", ncol(te_counts), " samples")
message("  Metadata: ", nrow(metadata), " samples, groups: ",
        paste(levels(metadata$group), collapse = ", "))

# =============================================================================
# CREATE COMBINED DGEList
# =============================================================================

message("\n[STEP 2] Creating combined gene+TE DGEList...")

dge <- load_or_compute(
    checkpoint_file = "1.5_te_dge_combined.rds",
    description = "combined gene+TE DGEList",
    compute_fn = function() {
        create_combined_dge(
            gene_counts = gene_counts,
            te_counts = te_counts,
            samples = metadata,
            validate = TRUE,      # Run all validation checks
            normalize = TRUE      # Apply TMM normalization
        )
    }
)

# Print summary
print_dge_summary(dge)

# =============================================================================
# SET UP DESIGN AND CONTRASTS
# =============================================================================

message("\n[STEP 3] Setting up experimental design...")

# Design matrix (cell-means parameterization)
design <- model.matrix(~ 0 + group, data = dge$samples)
colnames(design) <- gsub("^group", "", colnames(design))

message("  Design matrix columns: ", paste(colnames(design), collapse = ", "))

# Define contrasts
contrast_matrix <- limma::makeContrasts(
    Treatment_vs_Control = Treatment - Control,
    levels = design
)

message("  Contrasts: ", paste(colnames(contrast_matrix), collapse = ", "))

# =============================================================================
# FIT LINEAR MODEL
# =============================================================================

message("\n[STEP 4] Fitting linear model (limma-voom)...")

fit <- load_or_compute(
    checkpoint_file = "1.5_te_fit_combined.rds",
    description = "limma-voom fit on combined data",
    compute_fn = function() {
        compute_te_de(
            dge = dge,
            design = design,
            contrast_matrix = contrast_matrix,
            robust = TRUE
        )
    }
)

message("  Fit object: ", nrow(fit), " features x ", ncol(fit$coefficients), " contrasts")

# =============================================================================
# COMPUTE TE ENRICHMENT
# =============================================================================

message("\n[STEP 5] Computing TE enrichment...")

# Subfamily level
te_enrichment_sub <- load_or_compute(
    checkpoint_file = "1.5_te_enrichment_subfamily.rds",
    description = "TE subfamily enrichment",
    compute_fn = function() {
        compute_te_enrichment(
            fit = fit,
            dge = dge,
            te_level = "subfamily",
            exclude_groups = c("Unspecified", "Unknown"),
            min_size = 3
        )
    }
)

# Family level
te_enrichment_fam <- load_or_compute(
    checkpoint_file = "1.5_te_enrichment_family.rds",
    description = "TE family enrichment",
    compute_fn = function() {
        compute_te_enrichment(
            fit = fit,
            dge = dge,
            te_level = "family",
            exclude_groups = c("Unspecified", "Unknown"),
            min_size = 3
        )
    }
)

# Class level
te_enrichment_class <- load_or_compute(
    checkpoint_file = "1.5_te_enrichment_class.rds",
    description = "TE class enrichment",
    compute_fn = function() {
        compute_te_enrichment(
            fit = fit,
            dge = dge,
            te_level = "class",
            exclude_groups = c("Unspecified", "Unknown"),
            min_size = 3
        )
    }
)

# Print summary
summarize_te_enrichment(te_enrichment_class)

# =============================================================================
# COMPUTE TE AGGREGATES
# =============================================================================

message("\n[STEP 6] Computing TE aggregates...")

te_aggregates <- load_or_compute(
    checkpoint_file = "1.5_te_aggregates.rds",
    description = "TE expression aggregates",
    compute_fn = function() {
        compute_te_aggregates(
            dge = dge,
            te_level = "family",
            method = "mean",
            value_type = "logcpm"
        )
    }
)

message("  Aggregated groups: ", nrow(te_aggregates$matrix))

# =============================================================================
# CREATE MASTER TABLES
# =============================================================================

message("\n[STEP 7] Creating master tables...")

# Ensure output directory exists
dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)

# TE Enrichment master table
message("  Creating master_te_enrichment.csv...")
te_enrich_sub_fmt <- format_te_enrichment(te_enrichment_sub)
te_enrich_fam_fmt <- format_te_enrichment(te_enrichment_fam)
te_enrich_class_fmt <- format_te_enrichment(te_enrichment_class)

master_te_enrichment <- bind_rows(
    te_enrich_sub_fmt,
    te_enrich_fam_fmt,
    te_enrich_class_fmt
)

write_csv(master_te_enrichment, file.path(DIR_TABLES, "master_te_enrichment.csv"))
message("    Rows: ", nrow(master_te_enrichment))

# TE DE master table
message("  Creating master_te_de.csv...")
master_te_de <- format_te_de(fit, dge)
write_csv(master_te_de, file.path(DIR_TABLES, "master_te_de.csv"))
message("    Rows: ", nrow(master_te_de))

# TE Aggregates master table
message("  Creating master_te_aggregates.csv...")
master_te_aggregates <- format_te_aggregates(te_aggregates)
write_csv(master_te_aggregates, file.path(DIR_TABLES, "master_te_aggregates.csv"))
message("    Rows: ", nrow(master_te_aggregates))

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=================================================================")
message("Phase 1.5 Complete: TE Analysis (Combined Mode)")
message("=================================================================")
message("Completed: ", Sys.time())
message("")
message("Checkpoints saved:")
message("  - ", file.path(DIR_CHECKPOINTS, "1.5_te_dge_combined.rds"))
message("  - ", file.path(DIR_CHECKPOINTS, "1.5_te_fit_combined.rds"))
message("  - ", file.path(DIR_CHECKPOINTS, "1.5_te_enrichment_*.rds"))
message("  - ", file.path(DIR_CHECKPOINTS, "1.5_te_aggregates.rds"))
message("")
message("Master tables saved:")
message("  - ", file.path(DIR_TABLES, "master_te_enrichment.csv"))
message("  - ", file.path(DIR_TABLES, "master_te_de.csv"))
message("  - ", file.path(DIR_TABLES, "master_te_aggregates.csv"))
message("")
message("Next step: Run 2.5.te_visualization.R for plots")
