#!/usr/bin/env Rscript
# =============================================================================
# 1.5.te_processing_separate.R - TE Analysis Pipeline (Separate Mode)
# =============================================================================
# TE-RNAseq-toolkit Example Script
# Phase: 1.5
#
# Description:
#   Template script for TE-specific analysis using SEPARATE mode.
#   In separate mode, TEs are analyzed independently from genes, with
#   optional library size borrowing from the gene reference to prevent
#   normalization artifacts.
#
# Prerequisites:
#   - TE count matrix (TEs x samples)
#   - Gene DGEList from main pipeline (for library size reference)
#   - Sample metadata
#
# When to use Separate Mode:
#   - Gene and TE annotations may overlap
#   - You want genes and TEs as separate hypothesis universes
#   - Different preprocessing was used for genes vs TEs
#   - You cannot verify that exonic TEs were removed from SAF
#
# Outputs:
#   - Checkpoints: 1.5_te_dge_separate.rds, 1.5_te_fit_separate.rds, etc.
#   - Master tables: master_te_enrichment.csv, master_te_de.csv
#
# Usage:
#   1. Copy this file to your project's 02_analysis/scripts/ directory
#   2. Update paths and parameters for your project
#   3. Run: Rscript 1.5.te_processing_separate.R
#
# =============================================================================

message("=================================================================")
message("Phase 1.5: TE Analysis (Separate Mode)")
message("=================================================================")
message("Started: ", Sys.time())

# =============================================================================
# SETUP
# =============================================================================

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

PROJECT_ROOT <- "/workspaces/YOUR_PROJECT"  # <-- UPDATE THIS

# Data directories
DIR_DATA        <- file.path(PROJECT_ROOT, "00_data/processed")
DIR_CHECKPOINTS <- file.path(PROJECT_ROOT, "03_results/checkpoints")
DIR_TABLES      <- file.path(PROJECT_ROOT, "03_results/tables")

# Input files
TE_COUNTS_FILE  <- file.path(DIR_DATA, "te_counts.txt")          # <-- UPDATE
GENE_DGE_FILE   <- file.path(DIR_CHECKPOINTS, "1.1_dge_normalized.rds")  # From main pipeline
METADATA_FILE   <- file.path(DIR_DATA, "sample_metadata.csv")    # <-- UPDATE

# TE-RNAseq-toolkit location
TOOLKIT_PATH <- "/workspaces/AdaW_eWAT_WL_2025/01_modules/TE-RNAseq-toolkit/R"

# -----------------------------------------------------------------------------
# SOURCE TOOLKIT FUNCTIONS
# -----------------------------------------------------------------------------

message("[SETUP] Loading TE-RNAseq-toolkit functions...")

source(file.path(TOOLKIT_PATH, "te_utils.R"))
source(file.path(TOOLKIT_PATH, "validate_te_input.R"))
source(file.path(TOOLKIT_PATH, "create_separate_dge.R"))
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

load_or_compute <- function(checkpoint_file, description, compute_fn, force_recompute = FALSE) {

    checkpoint_path <- file.path(DIR_CHECKPOINTS, checkpoint_file)

    if (file.exists(checkpoint_path) && !force_recompute) {
        message("[CACHE] Loading ", description, " from: ", checkpoint_file)
        return(readRDS(checkpoint_path))
    }

    message("[COMPUTE] Computing ", description, "...")
    result <- compute_fn()

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

# For demonstration, create synthetic data
# REMOVE THIS SECTION and load your real data
set.seed(42)
n_genes <- 100
n_tes <- 20
n_samples <- 6

gene_ids <- paste0("Gene", 1:n_genes)
sample_names <- paste0("Sample", 1:n_samples)

# Synthetic gene counts (for reference)
gene_counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes,
    dimnames = list(gene_ids, sample_names)
)

# Synthetic TE counts
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

message("  TE counts: ", nrow(te_counts), " TEs x ", ncol(te_counts), " samples")
message("  Metadata: ", nrow(metadata), " samples, groups: ",
        paste(levels(metadata$group), collapse = ", "))

# =============================================================================
# CREATE GENE REFERENCE DGEList (for library size borrowing)
# =============================================================================

message("\n[STEP 2] Creating gene reference DGEList...")

# In a real project, you would load this from the main pipeline:
# dge_genes <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_dge_normalized.rds"))

# For demonstration, create from synthetic gene counts
dge_genes <- create_gene_dge_reference(
    gene_counts = gene_counts,
    samples = metadata,
    normalize = TRUE
)

message("  Gene reference: ", nrow(dge_genes), " genes x ", ncol(dge_genes), " samples")
message("  Library size range: ", round(min(dge_genes$samples$lib.size)),
        " - ", round(max(dge_genes$samples$lib.size)))

# =============================================================================
# CREATE SEPARATE TE DGEList (with library size borrowing)
# =============================================================================

message("\n[STEP 3] Creating separate TE DGEList with library size borrowing...")

dge_te <- load_or_compute(
    checkpoint_file = "1.5_te_dge_separate.rds",
    description = "separate TE DGEList",
    compute_fn = function() {
        create_separate_te_dge(
            te_counts = te_counts,
            samples = metadata,
            gene_dge_reference = dge_genes,  # Borrow library sizes
            validate = TRUE,
            normalize = TRUE
        )
    }
)

# Verify library size borrowing
message("  TE DGEList: ", nrow(dge_te), " TEs x ", ncol(dge_te), " samples")
message("  Library sizes borrowed: ", dge_te$lib_size_borrowed)
message("  Analysis mode: ", dge_te$analysis_mode)

# =============================================================================
# SET UP DESIGN AND CONTRASTS
# =============================================================================

message("\n[STEP 4] Setting up experimental design...")

design <- model.matrix(~ 0 + group, data = dge_te$samples)
colnames(design) <- gsub("^group", "", colnames(design))

message("  Design matrix columns: ", paste(colnames(design), collapse = ", "))

contrast_matrix <- limma::makeContrasts(
    Treatment_vs_Control = Treatment - Control,
    levels = design
)

message("  Contrasts: ", paste(colnames(contrast_matrix), collapse = ", "))

# =============================================================================
# FIT LINEAR MODEL (on TE data only)
# =============================================================================

message("\n[STEP 5] Fitting linear model (limma-voom) on TE data...")

fit_te <- load_or_compute(
    checkpoint_file = "1.5_te_fit_separate.rds",
    description = "limma-voom fit on TE data",
    compute_fn = function() {
        compute_te_de(
            dge = dge_te,
            design = design,
            contrast_matrix = contrast_matrix,
            robust = TRUE
        )
    }
)

message("  Fit object: ", nrow(fit_te), " features x ", ncol(fit_te$coefficients), " contrasts")

# =============================================================================
# COMPUTE TE ENRICHMENT
# =============================================================================

message("\n[STEP 6] Computing TE enrichment...")

# Family level (most informative for separate mode)
te_enrichment_fam <- load_or_compute(
    checkpoint_file = "1.5_te_enrichment_family_separate.rds",
    description = "TE family enrichment (separate mode)",
    compute_fn = function() {
        compute_te_enrichment(
            fit = fit_te,
            dge = dge_te,
            te_level = "family",
            exclude_groups = c("Unspecified", "Unknown"),
            min_size = 3
        )
    }
)

# Class level
te_enrichment_class <- load_or_compute(
    checkpoint_file = "1.5_te_enrichment_class_separate.rds",
    description = "TE class enrichment (separate mode)",
    compute_fn = function() {
        compute_te_enrichment(
            fit = fit_te,
            dge = dge_te,
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

message("\n[STEP 7] Computing TE aggregates...")

te_aggregates <- load_or_compute(
    checkpoint_file = "1.5_te_aggregates_separate.rds",
    description = "TE expression aggregates (separate mode)",
    compute_fn = function() {
        compute_te_aggregates(
            dge = dge_te,
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

message("\n[STEP 8] Creating master tables...")

dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)

# TE Enrichment master table
message("  Creating master_te_enrichment.csv...")
te_enrich_fam_fmt <- format_te_enrichment(te_enrichment_fam)
te_enrich_class_fmt <- format_te_enrichment(te_enrichment_class)

master_te_enrichment <- bind_rows(
    te_enrich_fam_fmt,
    te_enrich_class_fmt
)

write_csv(master_te_enrichment, file.path(DIR_TABLES, "master_te_enrichment.csv"))
message("    Rows: ", nrow(master_te_enrichment))
message("    Analysis mode: ", unique(master_te_enrichment$analysis_mode))

# TE DE master table
message("  Creating master_te_de.csv...")
master_te_de <- format_te_de(fit_te, dge_te)
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
message("Phase 1.5 Complete: TE Analysis (Separate Mode)")
message("=================================================================")
message("Completed: ", Sys.time())
message("")
message("KEY DIFFERENCES FROM COMBINED MODE:")
message("  - TEs analyzed as separate hypothesis universe")
message("  - Library sizes borrowed from gene reference")
message("  - FDR controlled within TE population only")
message("")
message("Checkpoints saved:")
message("  - ", file.path(DIR_CHECKPOINTS, "1.5_te_dge_separate.rds"))
message("  - ", file.path(DIR_CHECKPOINTS, "1.5_te_fit_separate.rds"))
message("  - ", file.path(DIR_CHECKPOINTS, "1.5_te_enrichment_*_separate.rds"))
message("")
message("Master tables saved:")
message("  - ", file.path(DIR_TABLES, "master_te_enrichment.csv"))
message("  - ", file.path(DIR_TABLES, "master_te_de.csv"))
message("  - ", file.path(DIR_TABLES, "master_te_aggregates.csv"))
message("")
message("Next step: Run 2.5.te_visualization.R for plots")
