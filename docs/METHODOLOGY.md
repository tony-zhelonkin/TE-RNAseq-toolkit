# Integrating Gene and Transposable Element RNA-seq Analysis

This document is an overview of the theoretical considerations, decision points, and alternative strategies for integrating Transposable Element (TE) analysis into RNA-seq workflows. It expands on the "Quick Start" in the README to cover *why* specific choices were made and *how* to interpret the results.

## Mapping and Counting Strategies for Transposable Elements (TEs)

**Multi-mapping reads require special handling:** Unlike most mRNAs, TE-derived reads often map to many loci due to repetitive sequences. This creates the **Multi-Mapping Problem**: standard RNA-seq aligners (which discard non-unique reads) will systematically blind you to the majority of TE expression. A standard RNA-seq alignment (which keeps only one alignment per read or discards multi-mappers) will **underestimate TE expression** [[1]]. 

### The Two Extremes of Bias
1.  **Unique-Only (Undercounting)**: If you discard multi-mappers (default in many gene pipelines), you lose ~80-90% of TE signal. This is generally **unacceptable** for TE analysis [[1]].
2.  **Total-Counts (Overcounting)**: If you count a read for *every* location it maps to (e.g., 20 locations = 20 counts), you inflate library sizes and skew statistical testing.

## 2. Alignment & Counting Strategies

There are three primary strategies to solve the multi-mapping problem. This toolkit supports **Strategy A** and **Strategy B** most natively.

### Strategy A: "Random-One" (Integer Counts)
Report one random alignment per multi-mapper. This is a middle ground: **Every sequenced fragment should contribute exactly 1 count to the final matrix.**
*   **Mechanism**: The aligner is configured to identify all valid map locations but report **only one** to the BAM file, chosen at random.
In STAR, use `--outMultimapperOrder Random` *and* limit each read to one alignment (e.g., `--outFilterMultimapNmax 100 --outSAMmultNmax 1 --outSAMprimaryFlag OneBestScore`).
*   **STAR Flags**: `--outFilterMultimapNmax 100 --outMultimapperOrder Random --outSAMmultNmax 1`, i.e. limit each read to one alignment
*   **Counting:** Downstream counting (e.g., featureCounts) can treat all reads normally (integers).
*   **Result:** Each read is assigned wholly to a single locus chosen at random.
*   **Accuracy:** This *random one-location* approach recovers TE expression accurately (correlation ~1 with true values benchmarks) [[1]].
*   **Pros**:
    *   Produces **integer counts**, compatible with all DE tools (DESeq2, edgeR, limma).
    *   Keeps BAM files small.
    *   Accurate for **family-level** quantification (random assignment averages out across thousands of reads) [[1]].
*   **Cons**: Locus-level precision is lost (stochastic assignment).

### Strategy B: Fractional Assignment (Weighted Counts)
Allow all alignments and count fractionally. 
*   **Mechanism**: The aligner reports *all* map locations. The counter splits the read weight (e.g., if a read maps to 4 loci, each gets 0.25 counts).
*   **STAR Flags**: STAR `--outFilterMultimapNmax 100 --outSAMmultNmax 100`, i.e. report **all high-scoring alignments**
*   **Counting:** Use a tool that splits counts fractionally (e.g., `featureCounts -M --fraction`). Each read contributes **1/n** to each of its *n* targets.
*   **Result:** Preserves all mapping possibilities.
*   **Accuracy:** Highly accurate for TE quantification [[1]].
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
*   **Cons**: often slower; might have limitations to the design complexity, e.g. TEtranscripts internally allows for only simple pair-wise comparisons of type A vs B.

**FeatureCounts and unsorted BAMs:**
STAR’s `--outSAMtype BAM Unsorted` is sufficient. FeatureCounts accepts unsorted BAMs and pairs reads internally, but can also work with sorted data, unlike TEtranscripts which is better suited for unsorted outputs [[2]].


---

## 3. Analysis Architecture: Combined vs. Separate

Should genes and TEs be analyzed in the same matrix?

### The Case for Combined Analysis (Recommended)
Literature suggests a combined approach can be beneficial [[1]].

### 1. Shared Normalization and Dispersion
*   **Normalization:** Including genes and TEs allows TMM/DESeq2 size factors to account for total library composition. If a sample has high TE content, separate normalization might mis-estimate sequencing depth for genes.
*   **Dispersion:** ~78k genes inform the mean-variance trend, stabilizing dispersion estimates for the smaller set of TEs.
*   **Shared Dispersion**: The ~78k gene features help estimate the variance trend for the ~1.2k TE features, stabilizing p-values.
### 2. Unified FDR Control
*   Testing genes and TEs in one model ensures a single family-wise error rate (FDR) across all genomic features.
*   Separate analysis yields separate FDRs (e.g., top 5% of TEs vs top 5% of genes), which are different hypothesis universes.
### 3. Practical Considerations
*   **Combined:** Requires non-overlapping annotations or fractional assignment logic. Best for rigorous, cohesive analysis.
*   **Separate:** Simpler to implement with off-the-shelf pipelines. Acceptable if interpreted carefully.

**Prerequisite for the Combined Analysis**: Annotations must not overlap.
*   *Solution*: This toolkit assumes you have subtracted exon regions from the TE SAF. This ensures no read is double-counted.

## Gene vs TE overlapping reads
Reads mapping to TEs within genes (introns/UTRs) can be double-counted.
*   **Combined Analysis:** A read overlapping a gene exon and a TE region could be counted for both if not handled.
*   **Solution:** Remove TE SAF regions that overlap gene exons (as done in this toolkit's recommended workflow) or use tools like TEtranscripts that prioritize gene assignment [[3]].
Most TE tools, like TEtranscripts, TEtranscripts** will not count a read toward a TE if it overlaps a gene exon, to avoid attributing host gene expression to TEs. 

### Locus-level vs family-level annotation feature quantification
*   **Subfamily Level (Recommended):** Most analyses focus on subfamily activity (e.g., all *L1Md_A* copies). Aggregating locus-level counts to subfamily improves robustness.
*   **Locus Level:** Testing individual loci is possible but has lower power due to multi-mapping signal dilution.
*   

### The Case for Separate Analysis
*   **Practicality**: Easier to implement if you have an existing gene count matrix and just want to "add on" TEs.
*   **Caveat**: Be mindful that "FDR < 0.05" means something different when testing 1,000 features vs 80,000 features.






---

## Normalization: Library Size and Composition

**Depth normalization (TMM/RLE)** is critical:
*   **TMM (edgeR):** Trims extreme genes to estimate a scaling factor. Robust to composition differences.
*   **RLE (DESeq2):** Uses median ratios to a geometric mean reference.
*   **Avoid FPKM/TPM for DE:** These are descriptive measures. DE tools expect counts (or estimated counts) to model variance correctly [[4]].

**TE Proportion as QC:**
Check the proportion of reads mapping to TEs. Wild variations between replicates can indicate technical issues. Consistent differences between conditions may reflect biology (TE derepression).

---

## Differential Expression Frameworks

### edgeR (QL or Classic)
*   **Model:** Negative Binomial GLM.
*   **Pros:** Fast, highly configurable, robust dispersion estimation. Quasi-Likelihood (QL) F-test controls false positives well.
*   **TE Suitability:** Handles integer counts natively. Can accept fractional counts (will treat as numeric), effective for TE analysis.

### DESeq2
*   **Model:** Negative Binomial GLM with empirical Bayes shrinkage.
*   **Pros:** User-friendly, built-in independent filtering, `lfcShrink` for visualization.
*   **TE Suitability:** Expects **integers**. Fractional counts must be rounded or scaled, making it less ideal for fractional TE pipelines unless adapted.

### limma-voom
*   **Model:** Linear model on log-CPM with precision weights (mean-variance modeling).
*   **Pros:** Flexible, fast, handles fractional counts natively (great for TE "multi-mode" counts). Extensive support for complex designs (random effects via `dream`).
*   **TE Suitability:** **Excellent** for fractional TE counts. The precision weights should handle the variance structure of low/high counts, but I would recommend comparing the influence on the conclusion between the two models. I had an experience when precision weights flipped log2FC signs in batch-affected samples. 

---

## 4. Statistical Frameworks for TEs

| Framework | Best For | Fractional Support | Notes |
| :--- | :--- | :--- | :--- |
| **edgeR (QL)** | **Strategy A** (Integer) | Partial | Can accept non-integers, but probability model assumes integers. Robust dispersion is excellent for outliers. |
| **DESeq2** | **Strategy A** (Integer) | No | **Strictly requires integers**. Must round fractional counts (loss of precision) or use tximport-like offsets. |
| **limma-voom** | **Strategy B** (Fractional) | **Yes** | **Gold standard for fractional counts**. Models mean-variance relationship of log-CPM weights, handling non-integers natively. |

**Toolkit Default**: The toolkit's `compute_te_de.R` uses **limma-voom** because it is universally safe: it handles integers (Strategy A) and fractions (Strategy B) equally well without warning.
---

## Summary Table: Analysis Workflows

| Workflow | Mapping/Counting | Normalization/DE | Pros | Cons |
| :--- | :--- | :--- | :--- | :--- |
| **Separate Analysis** | High multi-map limits. Separate runs (Genes: unique; TEs: fractional). | Separate TMM/DE runs. | Simple implementation. Optimizes strandedness per feature type. | Risk of double-counting. Normalization discordance. Separate FDRs. |
| **Combined Analysis** | Combined annotation (non-overlapping). | Unified TMM. Single DE model. | Unified library size & FDR. No double-counting. Robust dispersion. | Complexity in prep. Requires checking overlaps. |
| **Random-One** | One random alignment per read. | Integer counts (DESeq2/edgeR). | Simple integer workflow. Accurate for families. Smaller BAMs. | Loss of multi-mapping probabilistic info. Stochastic (needs seed). |
| **Fractional** | All alignments, weighted 1/n. | Fractional counts (limma-voom). | Max info retention. Statistically principled. | Requires tools supporting non-integers (voom). Large BAMs. |
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

1. **Mobile DNA (2019)**: [Tools and best practices for retrotransposon analysis](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1)
2. **Biostars**: [Can featureCounts use sortedbycoordinate files?](https://www.biostars.org/p/488290/)
3. **PMC (2015)**: [TEtranscripts: a package for including transposable elements](https://pmc.ncbi.nlm.nih.gov/articles/PMC4757950/)
4. **J Transl Med (2021)**: [TPM, FPKM, or Normalized Counts?](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w)
