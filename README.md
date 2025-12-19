# TE-RNAseq-toolkit

**Version:** 2.0.0
**Date:** 2025-12-18

A standardized R toolkit for Transposable Element (TE) analysis in bulk RNA-seq data. This module supports both **combined** (genes + TEs) and **separate** (TE-only) analysis modes, integrating seamlessly with the project's standard 5-phase workflow.

---

## Key Features

- **Dual-Mode Analysis**:
  - **Combined Mode**: Analyze genes and TEs in a single statistical framework (Recommended when annotations are non-overlapping).
  - **Separate Mode**: Analyze TEs independently with library size borrowing from gene counts.
- **Phase-Based Workflow**:
  - **Phase 1 (Compute)**: Expensive calculations (enrichment, aggregation) are cached via `load_or_compute()`.
  - **Phase 2 (Tables)**: Standardized "Master Tables" for cross-language compatibility.
  - **Phase 3 (Viz)**: Pure visualization functions that read from master tables.
- **Rich Visualization**:
  - Multi-contrast enrichment dotplots.
  - TE-specific heatmaps and volcano plots.
  - Correlation analysis with gene pathways.

## Documentation

For a deep dive into the theoretical background, mapping strategies (Random-One vs Fractional), and statistical frameworks (edgeR vs DESeq2 vs limma-voom), please read:

👉 **[Detailed Methodology & Best Practices](docs/METHODOLOGY.md)**

---

## Directory Structure

```
TE-RNAseq-toolkit/
├── R/                                  # Core R functions
│   ├── create_combined_dge.R           # Factory for combined mode
│   ├── create_separate_dge.R           # Factory for separate mode
│   ├── compute_te_*.R                  # Computational functions (Phase 1)
│   ├── format_te_*.R                   # Formatting functions (Phase 2)
│   └── plot_te_*.R                     # Visualization functions (Phase 3)
│
├── examples/                           # Template scripts
│   ├── 1.5.te_processing_combined.R    # Main pipeline (Combined)
│   ├── 1.5.te_processing_separate.R    # Main pipeline (Separate)
│   └── 2.5.te_visualization.R          # Visualization pipeline
│
└── inst/config/                        # Configuration templates
    └── te_config_template.yaml
```

---

## Preprocessing & Prerequisites

**Crucial:** This toolkit is designed for data processed with specific strategies to ensure valid quantification of repetitive elements.

### 1. Alignment (STAR Parameters)
For valid integer counting of TE subfamilies/families, the aligner must be configured to handle multi-mappers by retaining them but emitting only one alignment per read ("Random-One" strategy).

**Key STAR Flags:**
- `--outFilterMultimapNmax 100`: Allow reads to map to up to 100 loci (essential for TEs).
- `--outSAMmultNmax 1`: Emit only **one** alignment per read to the BAM.
- `--outMultimapperOrder Random`: Randomly select the reported alignment.
- `--runRNGseed 777`: Ensure reproducibility.

### 2. Annotation Preparation
To use **Combined Mode**, gene and TE annotations must not overlap. Overlapping regions lead to double-counting and inflated library sizes.

**Recommended Workflow:**
1.  **Grouped SAF**: Create a SAF file where `GeneID` is the TE group (e.g., `HAL1:L1:LINE`) rather than the unique locus ID.
2.  **Filter Exons**: Remove TE loci that overlap with annotated gene exons using `bedtools subtract`.

### 3. Feature Counting
This toolkit expects counts generated with specific handling for TEs vs Genes:

| Feature | Strandedness | Multi-mappers | Rationale |
|---------|--------------|---------------|-----------|
| **Genes** | Library-specific | Excluded | Standard RNA-seq practice. |
| **TEs** | Unstranded (`-s 0`) | Included (`-M`) | TEs are repetitive and often bidirectional. |

**Note on Combined Matrix:**
Even with asymmetric counting parameters, a combined matrix is valid **IF** exonic TEs were removed from the annotation. A multi-mapping read might be excluded from genes (due to ambiguity) but counted for TEs (due to `-M`), but it will never be counted for both.

---

## Quick Start

### 1. Source the Toolkit

In your project's `config.R`:

```r
source_toolkit("TE-RNAseq-toolkit")
```

### 2. Choose Your Mode

| Scenario | Recommended Mode | Rationale |
|----------|-----------------|-----------|
| Exonic TEs removed from SAF | **Combined** | Unified FDR, robust TMM normalization |
| Annotations overlap | **Separate** | Conservative, avoids double-counting |
| "TEtranscripts" pipeline | **Combined** | Designed for combined analysis |

### 3. Run the Pipeline

Copy the appropriate example script to your analysis folder:

```bash
cp 01_modules/TE-RNAseq-toolkit/examples/1.5.te_processing_combined.R 02_analysis/1.5.te_analysis.R
```

Run it:

```bash
Rscript 02_analysis/1.5.te_analysis.R
```

---

## Core Functions

### Factory Functions (Data Ingestion)
- `create_combined_dge(gene_counts, te_counts, ...)`: Creates a unified DGEList. Validates no ID collisions.
- `create_separate_te_dge(te_counts, gene_dge_reference, ...)`: Creates TE-only DGEList, borrowing library sizes from genes.

### Phase 1: Computation (Cached)
- `compute_te_enrichment(fit, dge, te_level, ...)`: Competitive gene set testing (camera/geneSetTest) for TE families/subfamilies.
- `compute_te_aggregates(dge, ...)`: Aggregates expression by TE family/subfamily.
- `compute_te_de(...)`: Limma-voom analysis specific to TEs.

### Phase 2: Master Tables (Export)
- `format_te_enrichment(...)`: Creates `master_te_enrichment.csv`.
- `format_te_de(...)`: Creates `master_te_de.csv`.

### Phase 3: Visualization (Plotting)
- `plot_te_enrichment_dotplot(...)`: Visualizes enrichment results across contrasts.
- `plot_te_heatmap(...)`: Heatmaps of aggregated TE expression.
- `plot_te_volcano(...)`: Volcano plots highlighting specific TE groups.

---

## Input Requirements

1.  **Gene Count Matrix**: Standard genes x samples matrix.
2.  **TE Count Matrix**: Rows must be TE identifiers. Recommended format: `Subfamily:Family:Class` (e.g., `L1Md_A:L1:LINE`).
3.  **Sample Metadata**: Must match columns of count matrices.

---

## Configuration

Analysis parameters should be defined in your project's `analysis_config.yaml`:

```yaml
te_analysis:
  mode: "combined"
  levels: ["subfamily", "family"]
  exclude_groups: ["Unspecified", "Unknown"]
```

See `inst/config/te_config_template.yaml` for a full example.