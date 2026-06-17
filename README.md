# TE-RNAseq-toolkit

**Version:** 2.0.3
**Date:** 2026-06-17

A toolkit for Transposable Element (TE) analysis in bulk RNA-seq data. 
This module supports both **combined** (genes + TEs) and **separate** (TE-only) analysis modes.

---

## Why this toolkit

It is deliberately minimal: STAR Random-One + featureCounts produce a clean **integer** count matrix of TE subfamilies, 
and the toolkit hands you that matrix and stops — **you own the statistical model.** 

Most TE tools bundle their differential expression and lock it to a two-group comparison (TEtranscripts, SQuIRE) 
or an additive condition-only design (TE-Seq); 
emitting a raw integer matrix instead lets the counts drop straight into edgeR / limma / DESeq2 
with *any* design — 2×2 factorial, interactions, blocking, random effects, arbitrary contrasts. 

That freedom is not unique to featureCounts — atena, Telescope, TEcount, and SQuIRE's `--table_only` all emit matrices 
you can model yourself; what this toolkit adds is a **clean integer** matrix (no rounding of fractional or EM weights) paired with transparent, evidence-graded methodology, a strand sense/antisense channel, and a documented QC gate.

The honest cost of that freedom: for a standard locus-level autonomy call, a baked pipeline (TE-Seq / Telescope) may serve you better — it owns the design you would otherwise rebuild by hand, and it gives you locus identity. Random-One trades locus-level resolution for subfamily-level robustness and simplicity. The crossover is named in [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md#design-freedom): a complex design on a subfamily / strand question lands here; a standard locus-level autonomy call lands on a baked pipeline.

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

For a deep dive into the theoretical background, mapping strategies (Random-One vs Fractional), 
and statistical frameworks (edgeR vs DESeq2 vs limma-voom), please read:

👉 **[Detailed Methodology & Best Practices](docs/METHODOLOGY.md)**

For the **biology** the design rests on — TEs, transcription strand, read-through, bidirectional transcription, derepression — see 
  
👉 **[Biology primer](docs/BIOLOGY.md)** (the biology single-source-of-truth).

Start with the result, then with **[Limitations](docs/LIMITATIONS.md)**: it holds the [escalation map](docs/LIMITATIONS.md#escalation-map) — what a result lets you claim (the [claim ladder](docs/LIMITATIONS.md#claim-ladder)) and when to reach for more (the [decision table](docs/LIMITATIONS.md#decision-table)). It names how far this toolkit takes you and when to escalate to a richer SAF (Simplified Annotation Format), an EM tool, or the wet lab.

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
For valid integer counting of TE subfamilies/families, the aligner must be configured to handle multi-mappers 
by retaining them but emitting only one alignment per read ("Random-One" strategy).

**Key STAR Flags:**
- `--outFilterMultimapNmax 100`: Allow reads to map to up to 100 loci (essential for TEs).
- `--outSAMmultNmax 1`: Emit only **one** alignment per read to the BAM.
- `--outMultimapperOrder Random`: Randomly select the reported alignment.
- `--runRNGseed 777`: Ensure reproducibility.

### 2. Annotation Preparation
To use **Combined Mode**, gene and TE annotations must not overlap. 
Overlapping regions lead to double-counting and inflated library sizes.

**Recommended Workflow:**
1.  **Grouped SAF**: Create a SAF file where `GeneID` is the TE group (e.g., `HAL1:L1:LINE`) rather than the unique locus ID.
2.  **Filter Exons**: Remove TE loci that overlap with annotated gene exons using `bedtools subtract`.

> **Exon subtraction is double-counting hygiene, not autonomy isolation.** It makes genes and TEs mutually exclusive so no read is double-counted — it is *orthogonal* to read-through, which is mostly *intronic*, and intronic TEs are not exonic, so they remain; separating autonomous (intergenic) from read-through-exposed (intronic) TEs is a later genic-context step (see [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md) §2 and the Rung-2 step in §5).

### 3. Feature Counting
This toolkit expects counts generated with specific handling for TEs vs Genes:

| Feature | Strandedness | Multi-mappers | Rationale |
|---------|--------------|---------------|-----------|
| **Genes** | Library-specific (`-s 2` for reverse dUTP) | Excluded | Standard RNA-seq practice. |
| **TEs** | Standalone: `-s 0` (matches the dominant tool's default). Joint stranded matrix: best-practice is to match the genes (`-s 2` here) | Included (`-M`) | See note below. |

**Note on TE strandedness — the field is SPLIT, there is no enchmarked standard:**
The two poles are both legitimate: the **dominant tool TEtranscripts/TEcount defaults to UNSTRANDED** (`--stranded no`), 
while modern ipelines (TE-Seq, Mobile DNA 2025) **recommend stranded** for a directional library so 
TE-derived expression can be resolved from gene expression.

- **For a stranded library + joint gene+TE matrix, the more principled route is stranded TEs matched to the genes** 
(`-s 2` for reverse dUTP). 
This is **best-practice / mechanistically motivated, grade B — not a proven-superior standard** 
[evidence: B — TE-Seq pipeline; the dominant tool TEtranscripts pools/defaults unstranded, so contested].

- **`-s 0` (unstranded) is a defensible standalone choice** — it matches the most-used tool's default. 
Caveat: "unstranded → better TE sensitivity" has **NEVER been directly benchmarked**; the only study simulating both modes 
(Savytska et al. 2022, Front Genet, doi:10.3389/fgene.2022.1026847) found **stranded FDR (54.9%) ≤ unstranded (58.7%)**. 
Unstranded is a **sensitivity-FOR-specificity trade, never quantified as a gain** [evidence: GAP].

- **Bidirectional/antisense TE biology is real but class-specific** — established for the L1 
antisense promoter/ORF0 and LTR/ERV (Criscione 2016; Faulkner 2009), 
weak/passive for SINE (Alu) and intronic passengers (grade C). 
To preserve genuine antisense signal, use a **stranded sense/antisense split** 
(count TEs `-s 1` and `-s 2` as separate feature sets, à la SQuIRE), not a collapse to `-s 0`, 
which discards strand and structurally inflates TE counts. 
The biology is laid out in [`docs/BIOLOGY.md`](docs/BIOLOGY.md).

- (Note: the Teissandier 2019 / Mobile DNA benchmark concerns **multimapper handling**, not strandedness, 
and does not bear on this choice.) 
Sources: 
SQuIRE doi:10.1093/nar/gkz081; 
TE-Seq doi:10.1186/s13100-025-00381-w; 
Savytska 2022 doi:10.3389/fgene.2022.1026847; 
TEtranscripts.

`-M` is **required** under the Random-One strategy (STAR keeps the `NH>1` flag on the single emitted line, 
so featureCounts discards it without `-M`). 
Be aware, though, that full-integer `-M` duplication for TEs against unique-only genes inflates 
the TE contribution to the library by copy number; no canonical TE tool uses full-integer duplication 
(they fractionate 1/N or EM-apportion). 
This compounds the cross-feature non-comparability described below — keep size factors anchored on genes.

**Note on Combined Matrix:**
Removing exonic TEs from the annotation eliminates **double-counting** 
(a multi-mapping read might be excluded from genes due to ambiguity but counted for TEs due to `-M`, 
yet it will never be counted for both). 
This is **necessary but not sufficient** for a valid joint matrix — and it is also not sufficient for *autonomy*: 
read-through TEs survive (see the Annotation-Preparation note above). 
On the measurement side, the `-s`/multimapper asymmetry leaves gene and TE rows on different measurement bases, 
so the combined matrix is valid for **within-feature-type, across-sample differential expression only**, 
and **only if DESeq2/edgeR size factors are estimated from the genes alone** 
(e.g. `estimateSizeFactors(controlGenes = isGene)` / TMM on the gene submatrix). 
It is **not** valid for gene-vs-TE magnitude comparison within a sample, and **"TE % of transcriptome" is not interpretable** 
as biology (treat it as a QC sanity band only).

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
| Exonic TEs removed from SAF | **Combined** | Unified FDR; gene-anchored size factors (estimate from genes only) |
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
