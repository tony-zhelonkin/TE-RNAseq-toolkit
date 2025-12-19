# Test Script for TE-RNAseq-toolkit
# Validates core functionality using synthetic data

# Set working directory to project root if needed
# setwd("/workspaces/AdaW_eWAT_WL_2025")

message("=================================================================")
message("STARTING TE-RNAseq-toolkit VALIDATION")
message("=================================================================")

# Source required files
base_dir <- "/workspaces/AdaW_eWAT_WL_2025/01_modules/TE-RNAseq-toolkit/R"
source(file.path(base_dir, "te_utils.R"))
source(file.path(base_dir, "validate_te_input.R"))
source(file.path(base_dir, "create_combined_dge.R"))
source(file.path(base_dir, "create_separate_dge.R"))
source(file.path(base_dir, "compute_te_enrichment.R"))
source(file.path(base_dir, "compute_te_aggregates.R"))

library(edgeR)
library(limma)
library(dplyr)
library(tibble)

# =============================================================================
# 1. GENERATE SYNTHETIC DATA
# =============================================================================
message("\n[TEST] Generating synthetic data...")

n_genes <- 100
n_tes <- 20
n_samples <- 6

# Gene IDs
gene_ids <- paste0("Gene", 1:n_genes)

# TE IDs (TEtranscripts format)
# Structure: subfamily:family:class
te_subfamilies <- c("L1Md_A", "L1Md_T", "IAP_1", "IAP_2", "MuLV")
te_families <- c("L1", "L1", "ERVK", "ERVK", "ERVL")
te_classes <- c("LINE", "LINE", "LTR", "LTR", "LTR")

te_ids <- paste(te_subfamilies, te_families, te_classes, sep=":")
# Replicate to get 20 TEs
te_ids <- rep(te_ids, length.out=n_tes)
# Make unique by adding suffix if needed, but here we want unique IDs
te_ids <- paste0(te_ids, "_", 1:n_tes)
# Actually, let's make them look real
te_real_ids <- c()
for(i in 1:n_tes) {
  fam_idx <- (i - 1) %% 5 + 1
  te_real_ids <- c(te_real_ids, paste(te_subfamilies[fam_idx], te_families[fam_idx], te_classes[fam_idx], sep=":"))
}
# TE names must be unique usually? No, TEtranscripts aggregates by subfamily usually.
# But here we are simulating locus-specific or unique IDs if we have a count matrix.
# Let's assume unique IDs like "L1Md_A:L1:LINE_dup1" is not how it looks.
# Usually it's "subfamily:family:class" and rows are unique subfamilies if using TEtranscripts count table.
# OR if locus specific: "chr1:100-200|subfamily:family:class"
# The parsing function handles "subfamily:family:class".
# Let's stick to unique identifiers that contain the info.
te_ids <- paste0("TE", 1:n_tes, "|", te_real_ids) 
# Wait, parse_te_id expects "subfamily:family:class" or it fails?
# Let's check parse_te_id in te_utils.R
# It splits by ":" and expects >= 3 parts.
# So "TE1|L1Md_A:L1:LINE" works?
# The code says: `parts <- strsplit(te_id, ":", fixed = TRUE)[[1]]`
# If I have "TE1|L1Md_A:L1:LINE", parts will be "TE1|L1Md_A", "L1", "LINE".
# This matches subfamily="TE1|L1Md_A", family="L1", class="LINE". This is valid.

sample_names <- paste0("Sample", 1:n_samples)
groups <- rep(c("Ctrl", "Treat"), each=3)
metadata <- data.frame(
    sample_id = sample_names,
    group = factor(groups),
    row.names = sample_names
)

# Counts
set.seed(123)
gene_counts <- matrix(rpois(n_genes * n_samples, lambda=100), nrow=n_genes)
rownames(gene_counts) <- gene_ids
colnames(gene_counts) <- sample_names

te_counts <- matrix(rpois(n_tes * n_samples, lambda=50), nrow=n_tes)
rownames(te_counts) <- te_ids
colnames(te_counts) <- sample_names

message("  Gene matrix: ", nrow(gene_counts), "x", ncol(gene_counts))
message("  TE matrix: ", nrow(te_counts), "x", ncol(te_counts))

# =============================================================================
# 2. TEST TE UTILS
# =============================================================================
message("\n[TEST] Testing te_utils.R...")

parsed <- parse_te_id(te_ids[1])
if(parsed$is_te && parsed$family == "L1" && parsed$class == "LINE") {
    message("[PASS] parse_te_id")
} else {
    stop("[FAIL] parse_te_id failed")
}

anno <- build_te_annotation(te_ids)
if(nrow(anno) == n_tes && all(anno$is_te)) {
    message("[PASS] build_te_annotation")
} else {
    stop("[FAIL] build_te_annotation failed")
}

# =============================================================================
# 3. TEST VALIDATION
# =============================================================================
message("\n[TEST] Testing validate_te_input.R...")

# Test 1: No collision
if(validate_no_id_collision(gene_counts, te_counts, stop_on_error=FALSE)) {
    message("[PASS] No collision check")
} else {
    stop("[FAIL] No collision check failed on valid data")
}

# Test 2: Sample consistency
if(validate_sample_consistency(gene_counts, te_counts, stop_on_error=FALSE)) {
    message("[PASS] Sample consistency check")
} else {
    stop("[FAIL] Sample consistency check failed on valid data")
}

# Test 3: Create collision
bad_gene_counts <- gene_counts
rownames(bad_gene_counts)[1] <- te_ids[1]
if(!validate_no_id_collision(bad_gene_counts, te_counts, stop_on_error=FALSE)) {
    message("[PASS] Collision detection works")
} else {
    stop("[FAIL] Collision detection missed a duplicate")
}

# =============================================================================
# 4. TEST COMBINED DGE
# =============================================================================
message("\n[TEST] Testing create_combined_dge.R...")

dge_combined <- create_combined_dge(
    gene_counts = gene_counts,
    te_counts = te_counts,
    samples = metadata,
    validate = TRUE
)

if(nrow(dge_combined) == n_genes + n_tes && dge_combined$analysis_mode == "combined") {
    message("[PASS] Combined DGE creation")
} else {
    stop("[FAIL] Combined DGE creation failed")
}

# =============================================================================
# 5. TEST SEPARATE DGE
# =============================================================================
message("\n[TEST] Testing create_separate_dge.R...")

# Create reference gene DGE
dge_ref <- create_gene_dge_reference(gene_counts, metadata)

dge_separate <- create_separate_te_dge(
    te_counts = te_counts,
    samples = metadata,
    gene_dge_reference = dge_ref
)

if(nrow(dge_separate) == n_tes && 
   dge_separate$analysis_mode == "separate" &&
   dge_separate$lib_size_borrowed == TRUE) {
    message("[PASS] Separate DGE creation with borrowing")
} else {
    stop("[FAIL] Separate DGE creation failed")
}

# =============================================================================
# 6. TEST ENRICHMENT
# =============================================================================
message("\n[TEST] Testing compute_te_enrichment.R...")

# Mock a fit object
# We need dge_combined to have a design and be processed
design <- model.matrix(~ group, data = metadata)
v <- voom(dge_combined, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Compute enrichment
# We only have one coef "groupTreat" (vs Ctrl)
enrichment <- compute_te_enrichment(
    fit = fit,
    dge = dge_combined,
    te_level = "class",
    contrasts = "groupTreat" # coef name in design
)

if(!is.null(enrichment$results) && nrow(enrichment$results) > 0) {
    message("[PASS] TE enrichment computation")
    print(head(enrichment$results))
} else {
    stop("[FAIL] TE enrichment returned empty results")
}

# =============================================================================
# 7. TEST AGGREGATES
# =============================================================================
message("\n[TEST] Testing compute_te_aggregates.R...")

agg <- compute_te_aggregates(
    dge = dge_combined,
    te_level = "family",
    method = "mean",
    value_type = "cpm"
)

if(!is.null(agg$matrix) && nrow(agg$matrix) > 0) {
    message("[PASS] TE aggregation")
    print(head(agg$matrix[,1:3]))
} else {
    stop("[FAIL] TE aggregation failed")
}

dge_agg <- create_aggregated_dge(
    dge = dge_combined,
    te_level = "class"
)

if(inherits(dge_agg, "DGEList") && dge_agg$analysis_mode == "aggregated") {
    message("[PASS] Aggregated DGE creation")
} else {
    stop("[FAIL] Aggregated DGE creation failed")
}

# =============================================================================
# 8. TEST COMPUTE TE DE
# =============================================================================
message("\n[TEST] Testing compute_te_de.R...")
source(file.path(base_dir, "compute_te_de.R"))

contrast_matrix <- makeContrasts(groupTreat, levels=design)
fit_de <- compute_te_de(dge_combined, design, contrast_matrix)

if(inherits(fit_de, "MArrayLM")) {
    message("[PASS] compute_te_de")
} else {
    stop("[FAIL] compute_te_de failed")
}

# =============================================================================
# 9. TEST FORMATTING FUNCTIONS
# =============================================================================
message("\n[TEST] Testing formatting functions...")
source(file.path(base_dir, "format_te_enrichment.R"))
source(file.path(base_dir, "format_te_aggregates.R"))
source(file.path(base_dir, "format_te_de.R"))

# Test format_te_enrichment
fmt_enrich <- format_te_enrichment(enrichment)
if(nrow(fmt_enrich) > 0 && "enrichment_score" %in% colnames(fmt_enrich)) {
    message("[PASS] format_te_enrichment")
} else {
    stop("[FAIL] format_te_enrichment failed")
}

# Test format_te_aggregates
fmt_agg <- format_te_aggregates(agg)
if(nrow(fmt_agg) > 0 && "expression" %in% colnames(fmt_agg)) {
    message("[PASS] format_te_aggregates")
} else {
    stop("[FAIL] format_te_aggregates failed")
}

# Test format_te_de
fmt_de <- format_te_de(fit_de, dge_combined)
if(nrow(fmt_de) > 0 && "logfc" %in% colnames(fmt_de)) {
    message("[PASS] format_te_de")
} else {
    stop("[FAIL] format_te_de failed")
}

# =============================================================================
# 10. TEST VISUALIZATION FUNCTIONS
# =============================================================================
message("\n[TEST] Testing visualization functions...")
source(file.path(base_dir, "plot_utils.R"))
source(file.path(base_dir, "plot_te_enrichment_dotplot.R"))
source(file.path(base_dir, "plot_te_heatmap.R"))
source(file.path(base_dir, "plot_te_volcano.R"))

# Test dotplot
p_dot <- plot_te_enrichment_dotplot(fmt_enrich)
if(inherits(p_dot, "ggplot")) {
    message("[PASS] plot_te_enrichment_dotplot")
} else {
    stop("[FAIL] plot_te_enrichment_dotplot failed")
}

# Test heatmap
# Need to ensure pheatmap doesn't try to open a device that fails
# But pheatmap returns the object invisibly
tryCatch({
    hm <- plot_te_heatmap(agg)
    message("[PASS] plot_te_heatmap")
}, error = function(e) {
    stop("[FAIL] plot_te_heatmap failed: ", e$message)
})

# Test volcano
# Create synthetic DE results with feature_type
de_res <- fmt_de
de_res$feature_type <- "te"
p_vol <- plot_te_volcano(de_res, highlight_te=TRUE)
if(inherits(p_vol, "ggplot")) {
    message("[PASS] plot_te_volcano")
} else {
    stop("[FAIL] plot_te_volcano failed")
}

# =============================================================================
# 11. TEST PROPORTIONS
# =============================================================================
message("\n[TEST] Testing compute_te_proportions.R...")
source(file.path(base_dir, "compute_te_proportions.R"))
source(file.path(base_dir, "plot_te_proportions.R"))

props <- compute_te_proportions(dge_combined, group_by="group")
if(!is.null(props$samples) && !is.null(props$te_composition)) {
    message("[PASS] compute_te_proportions")
} else {
    stop("[FAIL] compute_te_proportions failed")
}

p_prop <- plot_te_total_prop(props, group_by="group")
if(inherits(p_prop, "ggplot")) {
    message("[PASS] plot_te_total_prop")
} else {
    stop("[FAIL] plot_te_total_prop failed")
}

p_comp <- plot_te_composition(props, level="class", group_by="group")
if(inherits(p_comp, "ggplot")) {
    message("[PASS] plot_te_composition")
} else {
    stop("[FAIL] plot_te_composition failed")
}

# =============================================================================
# 12. TEST SINGLE ENRICHMENT PLOT
# =============================================================================
message("\n[TEST] Testing plot_te_enrichment_single.R...")
source(file.path(base_dir, "plot_te_enrichment_single.R"))

# Use existing fmt_enrich from Section 9
p_single <- plot_te_enrichment_single(fmt_enrich, contrast_name="groupTreat")
if(inherits(p_single, "ggplot")) {
    message("[PASS] plot_te_enrichment_single")
} else {
    stop("[FAIL] plot_te_enrichment_single failed")
}


# =============================================================================
# 13. TEST CORRELATIONS
# =============================================================================
message("\n[TEST] Testing compute_te_correlations.R...")
source(file.path(base_dir, "compute_te_correlations.R"))
source(file.path(base_dir, "plot_te_correlation.R"))

# Create synthetic pathway matrix
pathways <- c("Glycolysis", "TCA_Cycle", "Hypoxia")
path_mat <- matrix(rnorm(length(pathways) * n_samples), nrow=length(pathways), dimnames=list(pathways, sample_names))

cor_res <- compute_te_correlations(agg, path_mat)
if(!is.null(cor_res$cor_matrix) && nrow(cor_res$cor_matrix) == nrow(agg$matrix)) {
    message("[PASS] compute_te_correlations")
} else {
    stop("[FAIL] compute_te_correlations failed")
}

tryCatch({
    p_cor <- plot_te_correlation(cor_res)
    message("[PASS] plot_te_correlation")
}, error = function(e) {
    stop("[FAIL] plot_te_correlation failed: ", e$message)
})

message("\n=================================================================")
message("ALL TESTS PASSED SUCCESSFULLY")
message("=================================================================")
