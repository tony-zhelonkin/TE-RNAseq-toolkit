# Integrating Gene and Transposable Element RNA-seq Analysis

## Mapping and Counting Strategies for Transposable Elements (TEs)

**Multi-mapping reads require special handling:** Unlike most mRNAs, TE-derived reads often map to many loci due to repetitive sequences. A standard RNA-seq alignment (which keeps only one alignment per read or discards multi-mappers) will **underestimate TE expression** [[1]]. To address this, two broadly equivalent strategies are recommended:

### 1. Report one random alignment per multi-mapper
In STAR, use `--outMultimapperOrder Random` *and* limit each read to one alignment (e.g., `--outFilterMultimapNmax 100 --outSAMmultNmax 1 --outSAMprimaryFlag OneBestScore`).
*   **Result:** Each read is assigned wholly to a single locus chosen at random.
*   **Counting:** Downstream counting (e.g., featureCounts) can treat all reads normally (integers).
*   **Accuracy:** This *random one-location* approach recovers TE expression accurately (correlation ~1 with true values benchmarks) [[1]].

### 2. Allow all alignments and count fractionally
Configure STAR to report **all high-scoring alignments** (e.g., `--outFilterMultimapNmax=100`, `--winAnchorMultimapNmax=200`).
*   **Result:** Preserves all mapping possibilities.
*   **Counting:** Use a tool that splits counts fractionally (e.g., featureCounts `-M --fraction`). Each read contributes **1/n** to each of its *n* targets.
*   **Accuracy:** Highly accurate for TE quantification [[1]].

**Important STAR settings:**
*   `--outMultimapperOrder Random --runRNGseed 777`: Randomize alignment order (reproducible).
*   `--outFilterMultimapNmax 100`: Allow many mappings.

**FeatureCounts and unsorted BAMs:**
STAR’s `--outSAMtype BAM Unsorted` is sufficient. FeatureCounts accepts unsorted BAMs and pairs reads internally [[2]].

### Gene vs TE overlapping reads
Reads mapping to TEs within genes (introns/UTRs) can be double-counted.
*   **Combined Analysis:** A read overlapping a gene exon and a TE region could be counted for both if not handled.
*   **Solution:** Remove TE SAF regions that overlap gene exons (as done in this toolkit's recommended workflow) or use tools like TEtranscripts that prioritize gene assignment [[3]].

### Locus-level vs family-level counts
*   **Subfamily Level (Recommended):** Most analyses focus on subfamily activity (e.g., all *L1Md_A* copies). Aggregating locus-level counts to subfamily improves robustness.
*   **Locus Level:** Testing individual loci is possible but has lower power due to multi-mapping signal dilution.

---

## Combining Gene and TE Count Data: Separate vs Unified Analysis

Literature suggests a combined approach can be beneficial [[1]]:

### 1. Shared Normalization and Dispersion
*   **Normalization:** Including genes and TEs allows TMM/DESeq2 size factors to account for total library composition. If a sample has high TE content, separate normalization might mis-estimate sequencing depth for genes.
*   **Dispersion:** ~78k genes inform the mean-variance trend, stabilizing dispersion estimates for the smaller set of TEs.

### 2. Unified FDR Control
*   Testing genes and TEs in one model ensures a single family-wise error rate (FDR) across all genomic features.
*   Separate analysis yields separate FDRs (e.g., top 5% of TEs vs top 5% of genes), which are different hypothesis universes.

### 3. Practical Considerations
*   **Combined:** Requires non-overlapping annotations or fractional assignment logic. Best for rigorous, cohesive analysis.
*   **Separate:** Simpler to implement with off-the-shelf pipelines. Acceptable if interpreted carefully.

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
*   **TE Suitability:** **Excellent** for fractional TE counts. The precision weights handle the variance structure of low/high counts effectively.

---

## Summary Table: Analysis Workflows

| Workflow | Mapping/Counting | Normalization/DE | Pros | Cons |
| :--- | :--- | :--- | :--- | :--- |
| **Separate Analysis** | High multi-map limits. Separate runs (Genes: unique; TEs: fractional). | Separate TMM/DE runs. | Simple implementation. Optimizes strandedness per feature type. | Risk of double-counting. Normalization discordance. Separate FDRs. |
| **Combined Analysis** | Combined annotation (non-overlapping). | Unified TMM. Single DE model. | Unified library size & FDR. No double-counting. Robust dispersion. | Complexity in prep. Requires checking overlaps. |
| **Random-One** | One random alignment per read. | Integer counts (DESeq2/edgeR). | Simple integer workflow. Accurate for families. Smaller BAMs. | Loss of multi-mapping probabilistic info. Stochastic (needs seed). |
| **Fractional** | All alignments, weighted 1/n. | Fractional counts (limma-voom). | Max info retention. Statistically principled. | Requires tools supporting non-integers (voom). Large BAMs. |

---

## References

1. **Mobile DNA (2019)**: [Tools and best practices for retrotransposon analysis](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1)
2. **Biostars**: [Can featureCounts use sortedbycoordinate files?](https://www.biostars.org/p/488290/)
3. **PMC (2015)**: [TEtranscripts: a package for including transposable elements](https://pmc.ncbi.nlm.nih.gov/articles/PMC4757950/)
4. **J Transl Med (2021)**: [TPM, FPKM, or Normalized Counts?](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w)
