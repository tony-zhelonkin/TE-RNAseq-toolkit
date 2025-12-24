# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

TE-RNAseq-toolkit (v2.0.0) is an R package for analyzing Transposable Elements (TEs) in bulk RNA-seq data. It supports both **combined** (genes + TEs in one matrix) and **separate** (TE-only) analysis modes, following a 3-phase workflow: Compute → Tables → Visualization.

## Running Tests

```bash
# Run the synthetic data validation suite
Rscript tests/test_validation_synthetic.R

# Run integration test with real data (requires project data)
Rscript tests/test_integration_real_data.R
```

## Running Example Pipelines

```bash
# Copy example to your project and configure paths, then run:
Rscript examples/1.5.te_processing_combined.R
Rscript examples/1.5.te_processing_separate.R
Rscript examples/2.5.te_visualization.R
```

## Architecture

### Three-Phase Design
1. **Phase 1 (Compute)**: `R/compute_*.R` - Expensive calculations cached via `load_or_compute()`
2. **Phase 2 (Tables)**: `R/format_*.R` - Export to standardized "Master Tables" (CSV)
3. **Phase 3 (Viz)**: `R/plot_*.R` - Pure visualization, reads from checkpoints/tables only

### Core Function Categories

| Category | Files | Purpose |
|----------|-------|---------|
| Factory | `create_combined_dge.R`, `create_separate_dge.R` | Build DGEList objects |
| Utils | `te_utils.R`, `validate_te_input.R` | Parsing, annotation, validation |
| Compute | `compute_te_de.R`, `compute_te_enrichment.R`, `compute_te_aggregates.R`, `compute_te_proportions.R`, `compute_te_correlations.R` | Analysis |
| Format | `format_te_de.R`, `format_te_enrichment.R`, `format_te_aggregates.R` | Master table export |
| Plot | `plot_te_*.R`, `plot_utils.R` | Visualization |

### Key Design Decisions

- **TE ID Format**: Expects TEtranscripts-style `subfamily:family:class` (e.g., `L1Md_A:L1:LINE`)
- **Combined Mode**: Unified TMM normalization and FDR control; requires non-overlapping annotations
- **Separate Mode**: TE-only analysis with library size borrowing from genes
- **Statistical Framework**: Uses limma-voom (handles both integer and fractional counts)
- **Enrichment Testing**: Uses `limma::geneSetTest` for competitive gene set testing

## Sourcing the Toolkit

```r
# Source all core functions
base_dir <- "path/to/TE-RNAseq-toolkit/R"
source(file.path(base_dir, "te_utils.R"))
source(file.path(base_dir, "validate_te_input.R"))
source(file.path(base_dir, "create_combined_dge.R"))
source(file.path(base_dir, "compute_te_enrichment.R"))
source(file.path(base_dir, "compute_te_aggregates.R"))
source(file.path(base_dir, "compute_te_de.R"))
source(file.path(base_dir, "format_te_enrichment.R"))
source(file.path(base_dir, "format_te_de.R"))
source(file.path(base_dir, "format_te_aggregates.R"))
```

## Common Workflow Pattern

```r
# 1. Create combined DGEList
dge <- create_combined_dge(gene_counts, te_counts, samples)

# 2. Set up design and contrasts
design <- model.matrix(~ 0 + group, data = dge$samples)
contrast_matrix <- makeContrasts(Treatment_vs_Control = Treatment - Control, levels = design)

# 3. Fit model
fit <- compute_te_de(dge, design, contrast_matrix)

# 4. Compute enrichment
enrichment <- compute_te_enrichment(fit, dge, te_level = "family")

# 5. Format and export
write_csv(format_te_enrichment(enrichment), "master_te_enrichment.csv")
```

## Dependencies

**Required**: edgeR, limma, dplyr, tibble, stringr, readr
**Visualization**: ggplot2, pheatmap, RColorBrewer, scales

## Configuration

Copy `inst/config/te_config_template.yaml` to your project and customize analysis parameters.

## Key Constants (te_utils.R)

```r
TE_CYTOPLASMIC_CLASSES <- c("LINE", "SINE", "LTR", "Retroposon")  # Retrotransposons
TE_NUCLEAR_CLASSES <- c("DNA", "RC")                               # DNA transposons
```

## Preprocessing Requirements

Before using this toolkit, ensure:
1. **STAR alignment** used `--outFilterMultimapNmax 100 --outSAMmultNmax 1 --outMultimapperOrder Random` for Random-One strategy
2. **TE annotation** has exonic regions removed (bedtools subtract) for Combined Mode
3. **featureCounts** used `-M` for TEs and `-s 0` (unstranded) for TE counting
