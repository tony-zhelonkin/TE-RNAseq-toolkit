# Theoretical Background & Methodology

This document provides a comprehensive overview of the theoretical considerations, decision points, and alternative strategies for integrating Transposable Element (TE) analysis into RNA-seq workflows. It expands on the "Quick Start" in the README to cover *why* specific choices were made and *how* to interpret the results.

---

## 1. The Challenge of Repetitive Elements

Unlike unique genes, TEs are present in hundreds to thousands of nearly identical copies across the genome. This creates the **Multi-Mapping Problem**: standard RNA-seq aligners (which discard non-unique reads) will systematically blind you to the majority of TE expression.

### The Two Extremes of Bias
1.  **Unique-Only (Undercounting)**: If you discard multi-mappers (default in many gene pipelines), you lose ~80-90% of TE signal. This is generally **unacceptable** for TE analysis [[1]].
2.  **Total-Counts (Overcounting)**: If you count a read for *every* location it maps to (e.g., 20 locations = 20 counts), you inflate library sizes and skew statistical testing.

We need a middle ground: **Every sequenced fragment should contribute exactly 1 count to the final matrix.**

---

## 2. Alignment & Counting Strategies

There are three primary strategies to solve the multi-mapping problem. This toolkit supports **Strategy A** and **Strategy B** most natively.

### Strategy A: "Random-One" (Integer Counts)
***Current Project Implementation** (via STAR parameters)*

*   **Mechanism**: The aligner is configured to identify all valid map locations but report **only one** to the BAM file, chosen at random.
*   **STAR Flags**: `--outFilterMultimapNmax 100 --outMultimapperOrder Random --outSAMmultNmax 1`.
*   **Counting**: Standard `featureCounts` (integers).
*   **Pros**:
    *   Produces **integer counts**, compatible with all DE tools (DESeq2, edgeR, limma).
    *   Keeps BAM files small.
    *   Accurate for **family-level** quantification (random assignment averages out across thousands of reads) [[1]].
*   **Cons**: Locus-level precision is lost (stochastic assignment).

### Strategy B: Fractional Assignment (Weighted Counts)
*   **Mechanism**: The aligner reports *all* map locations. The counter splits the read weight (e.g., if a read maps to 4 loci, each gets 0.25 counts).
*   **Flags**: STAR `--outFilterMultimapNmax 100 --outSAMmultNmax 100` + `featureCounts -M --fraction`.
*   **Pros**:
    *   Deterministic (no random seed needed).
    *   Maximizes information retention.
    *   Theoretically ideal for **family-level** summing.
*   **Cons**:
    *   Produces **non-integer counts**.
    *   Incompatible with integer-only tools (DESeq2 requires rounding/hacks).
    *   Huge BAM files.

### Strategy C: Expectation-Maximization (EM)
*   **Mechanism**: Probabilistic tools (TEtranscripts, Telescope, SalmonTE) iteratively reassign reads based on local coverage abundance.
*   **Pros**: The "Gold Standard" for **locus-specific** resolution (Telescope).
*   **Cons**: Complex to set up; often slower; requires specific reference builds.

---

## 3. Analysis Architecture: Combined vs. Separate

Should genes and TEs be analyzed in the same matrix?

### The Case for Combined Analysis (Recommended)
Literature recommends calling DE on genes and TEs together [[1]].
*   **Unified Normalization**: TMM/Size Factors account for the *total* sequencing depth. If a sample has massive TE derepression (e.g., 20% of reads), separate normalization might artificially inflate gene counts in that sample. Combined normalization corrects for this.
*   **Shared Dispersion**: The ~78k gene features help estimate the variance trend for the ~1.2k TE features, stabilizing p-values.
*   **Single FDR**: Controls the false discovery rate across the entire hypothesis universe.

**Prerequisite**: Annotations must not overlap.
*   *Solution*: This toolkit assumes you have subtracted exon regions from the TE SAF (as done in the `preprocessing` scripts). This ensures no read is double-counted.

### The Case for Separate Analysis
*   **Practicality**: Easier to implement if you have an existing gene count matrix and just want to "add on" TEs.
*   **Independence**: Allows optimizing strandedness differently (e.g., stranded for Genes, unstranded for TEs).
*   **Caveat**: Be mindful that "FDR < 0.05" means something different when testing 1,000 features vs 80,000 features.

---

## 4. Statistical Frameworks for TEs

| Framework | Best For | Fractional Support | Notes |
| :--- | :--- | :--- | :--- |
| **edgeR (QL)** | **Strategy A** (Integer) | Partial | Can accept non-integers, but probability model assumes integers. Robust dispersion is excellent for outliers. |
| **DESeq2** | **Strategy A** (Integer) | No | **Strictly requires integers**. Must round fractional counts (loss of precision) or use tximport-like offsets. |
| **limma-voom** | **Strategy B** (Fractional) | **Yes** | **Gold standard for fractional counts**. Models mean-variance relationship of log-CPM weights, handling non-integers natively. |

**Toolkit Default**: The toolkit's `compute_te_de.R` uses **limma-voom** because it is universally safe: it handles integers (Strategy A) and fractions (Strategy B) equally well without warning.

---

## 5. Biological Interpretation & Hypothesis Generation

Finding differentially expressed TEs is just the start. The next step is connecting them to biological mechanism.

### Interpreting Upregulation
If you see global upregulation of TE subfamilies (e.g., ERVs, L1s), it suggests a breakdown in silencing mechanisms or active immune signaling.
*   **Loss of Repression**: Check for downregulation of known repressors in your gene DE results:
    *   *Epigenetic*: `Setdb1`, `Kap1` (Trim28), `Dnmt1/3a/3b`, `Hdac`s.
    *   *RNAi*: `Mael`, `Piwil` genes (piRNA pathway).
*   **Immune Activation**: TEs often mimic viral infection (dsRNA). Check for upregulation of Interferon Stimulated Genes (ISGs) (e.g., `Oas`, `Ifit`, `Stat1`).

### Linking TEs to Pathways
This toolkit provides functions to correlate TE activity with biological pathways (`compute_te_correlations.R`).
*   **Method**: Correlate the "Eigengene" (PC1) of a TE family with the GSVA/ssGSEA scores of gene pathways.
*   **Example**: "Is the upregulation of *MERVL* correlated with the *WP_ELECTRON_TRANSPORT_CHAIN* pathway across samples?"

### Module Analysis (Advanced)
Since you have combined data, you can perform **WGCNA** (Weighted Gene Co-expression Network Analysis) on the union of genes and TEs.
*   **Hypothesis**: Do TEs form their own co-expression modules, or do they cluster with stress-response gene modules?
*   **Literature**: Studies (e.g., in cotton under stress) have found distinct TE co-expression modules, hinting at specific regulatory programs rather than random noise.

---

## 6. Gaps & Limitations

Be transparent about limitations in your reporting:
*   **Family vs. Locus**: This pipeline focuses on **family/subfamily** enrichment. It cannot definitively claim *which* specific genomic insertion is active (e.g., "L1Md_A on Chr1" vs "L1Md_A on Chr5"). Use **Telescope** if locus resolution is critical.
*   **Ambiguity**: A read mapping to a TE is a probability, not a certainty. "Random-One" assignment works on average but is stochastic for individual reads.
*   **Sensitivity**: Unique-only mapping (if you were to use it) is "conservative" but creates a massive False Negative rate. We accept the uncertainty of multi-mapping to gain sensitivity.

---

## 7. References

1.  **Best Practices Review**: [Tools and best practices for retrotransposon analysis using high-throughput sequencing data](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1) (Mobile DNA, 2019).
2.  **TEtranscripts**: [Including transposable elements in differential expression analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC4757950/) (PMC, 2015).
3.  **Normalization**: [TPM, FPKM, or Normalized Counts?](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w) (J Transl Med, 2021).
4.  **DESeq2 vs edgeR**: [Systematic benchmarking of statistical methods](https://academic.oup.com/bib/article/24/1/bbac612/6966517) (Briefings in Bioinformatics, 2023).