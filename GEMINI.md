# GEMINI.md - TE-RNAseq-toolkit Context

**Project:** TE-RNAseq-toolkit
**Version:** 2.0.0
**Type:** R Analysis Toolkit (Modular Scripts)
**Primary Goal:** Robust differential expression and enrichment analysis of Transposable Elements (TEs) alongside genes in bulk RNA-seq.

---

## 1. Project Overview

This toolkit provides a standardized, statistically rigorous framework for integrating TE analysis into RNA-seq pipelines. It addresses the specific challenges of TE quantification (multi-mapping reads) and differential expression (low counts, high variance).

### Core Philosophies
1.  **Multi-Mapping Handling:** Supports "Random-One" (integer counts) and "Fractional" (weighted counts) strategies.
2.  **Unified Statistics:** "Combined Mode" (Genes + TEs) allows TEs to borrow strength from gene dispersion estimates, providing more stable FDR control.
3.  **Phased Execution:** Separates expensive computation (Phase 1) from data export (Phase 2) and visualization (Phase 3).
4.  **Master Tables:** All visualizations are derived from standardized CSV exports, ensuring consistency between analysis and reporting.

---

## 2. Directory Structure

| Path | Purpose |
|------|---------|
| `R/` | Core function library. **Do not modify** unless developing the toolkit. |
| `examples/` | Template scripts for user analysis (`1.5.te_processing_*.R`). |
| `inst/config/` | Configuration templates (`te_config_template.yaml`). |
| `docs/` | Methodological background and restructuring plans. |
| `tests/` | Unit and integration tests. |

---

## 3. Analysis Workflows

### Mode A: Combined (Recommended)
Analyze genes and TEs in a single `DGEList`.
*   **Prerequisite:** Gene and TE annotations must **not** overlap (exonic TEs removed).
*   **Method:** Unified TMM normalization and `limma-voom` modeling.
*   **Benefit:** Shared dispersion estimation, single FDR correction.

### Mode B: Separate
Analyze TEs in a separate `DGEList`.
*   **Prerequisite:** None (can use overlapping annotations).
*   **Method:** Can optionally borrow library sizes from gene counts for normalization.
*   **Benefit:** Simpler to add to existing pipelines; conservative.

---

## 4. Operational Workflow (The "Phased" Approach)

Agents should follow this pattern when creating or modifying analysis scripts.

### Phase 1: Compute (Cached)
*   **Input:** Raw count matrices (Genes, TEs), Metadata.
*   **Action:** Create DGEList, Normalize, Fit Model, Compute Enrichment.
*   **Output:** `.rds` checkpoints (managed by `load_or_compute` helper).
*   **Key Functions:** `create_combined_dge`, `compute_te_de`, `compute_te_enrichment`.

### Phase 2: Tables (Export)
*   **Input:** Phase 1 RDS objects.
*   **Action:** Format results into human-readable tables.
*   **Output:** `master_te_de.csv`, `master_te_enrichment.csv`.
*   **Key Functions:** `format_te_de`, `format_te_enrichment`.

### Phase 3: Visualization (Read-Only)
*   **Input:** Master CSVs (preferred) or RDS objects.
*   **Action:** Generate Plots (Volcano, Dotplot, Heatmap).
*   **Output:** PDF/PNG figures.
*   **Key Functions:** `plot_te_enrichment_dotplot`, `plot_te_volcano`.

---

## 5. Development Conventions

### Script Setup
All analysis scripts must source the toolkit functions manually (until this is a proper package):
```r
TOOLKIT_PATH <- "path/to/TE-RNAseq-toolkit/R"
source(file.path(TOOLKIT_PATH, "te_utils.R"))
source(file.path(TOOLKIT_PATH, "create_combined_dge.R"))
# ... source other modules as needed
```

### The `load_or_compute` Pattern
Expensive steps **must** be wrapped in this local helper function to prevent unnecessary re-computation during iterative development.
```r
result <- load_or_compute(
    checkpoint_file = "1.5_te_result.rds",
    description = "description",
    compute_fn = function() { ... }
)
```

### TE Identifiers
The toolkit expects TE IDs in the format: `Subfamily:Family:Class`
*   Example: `L1Md_A:L1:LINE`
*   Parsing: Use `te_utils.R` functions (`parse_te_id`, `build_te_annotation`).

### Configuration
Use `inst/config/te_config_template.yaml` as the source of truth for analysis parameters (grouping levels, exclusions, p-value adjustments).

---

## 6. Key Files to Know

*   `R/te_utils.R`: Core parsing logic (`parse_te_id`).
*   `R/create_combined_dge.R`: Data ingestion and validation.
*   `R/compute_te_de.R`: The statistical heart (limma-voom wrapper).
*   `examples/1.5.te_processing_combined.R`: The **golden standard** implementation of the workflow.

## 7. Common Tasks & Commands

**Q: How do I start a new analysis?**
A: Copy `examples/1.5.te_processing_combined.R` to your analysis folder and update paths.

**Q: How do I change the statistical model?**
A: `compute_te_de.R` uses `limma::voom`. Modifications to the model formula (e.g., adding batch effects) happen in the script where `design` is defined, passed to `compute_te_de`.

**Q: What if I have fractional counts?**
A: The toolkit uses `limma-voom`, which handles fractional counts natively. No special flag is needed, but ensure you are using the correct count matrix input.
