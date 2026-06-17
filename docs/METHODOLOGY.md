# Integrating Gene and Transposable Element RNA-seq Analysis

This document is an overview of the theoretical considerations, decision points, and alternative strategies for 
integrating Transposable Element (TE) analysis into RNA-seq workflows. 

It expands on the "Quick Start" in the README to cover *why* specific choices were made and *how* to interpret the results.

> **Last updated:** 2026-06-11

## Mapping and Counting Strategies for Transposable Elements (TEs)

**Multi-mapping reads require special handling:** Unlike most mRNAs, TE-derived reads often map to many loci due to repetitive sequences. 
This creates the **Multi-Mapping Problem**: standard RNA-seq aligner parameters (which discard non-unique reads) will by default
systematically blind you to the majority of TE expression. 

A standard RNA-seq alignment (which keeps only one alignment per read or discards multi-mappers) will **underestimate TE expression** [[1]]. 

### The Two Extremes of Bias
1.  **Unique-Only (Undercounting)**: If you discard multi-mappers (default in many gene pipelines), 
you lose a substantial, **young-family-biased** fraction of TE signal — documented multimapper fractions are ~10–21% overall (ENCODE analysis) 
and much higher for young families (e.g. L1HS/AluY/SVA), while old, diverged families are affected modestly. 
The exact loss is family/read-length/library dependent, so a single blanket percentage is not warranted; 
the direction (unique-only systematically undercounts young families) is robust. 
This is generally **unacceptable** for TE analysis [[1]].
2.  **Total-Counts (Overcounting)**: If you count a read for *every* location it maps to (e.g., 20 locations = 20 counts), 
you inflate library sizes and skew statistical testing.

## 2. Alignment & Counting Strategies

There are three primary strategies to solve the multi-mapping problem. 
This toolkit supports **Strategy A** and **Strategy B** most natively.

### Strategy A: "Random-One" (Integer Counts)
Report one random alignment per multi-mapper. 
This is a middle ground: **Every sequenced fragment should contribute exactly 1 count to the final matrix.**
*   **Mechanism**: The aligner is configured to identify all valid map locations but report **only one** to the BAM file, chosen at random.
In STAR, use `--outMultimapperOrder Random` *and* limit each read to one alignment (e.g., `--outFilterMultimapNmax 100 --outSAMmultNmax 1 --outSAMprimaryFlag OneBestScore`).
*   **STAR Flags**: `--outFilterMultimapNmax 100 --outMultimapperOrder Random --outSAMmultNmax 1`, i.e. limit each read to one alignment
*   **Counting:** Downstream counting (e.g., featureCounts) can treat all reads normally (integers).
*   **Result:** Each read is assigned wholly to a single locus chosen at random.
*   **Why Random-One works:** Teissandier 2019 benchmarked multimapper handling in featureCounts and found that random-one integer assignment performs on par with fractional (1/n) for family- and subfamily-level quantification, while unique-only undercounts young families (grade A — see the grading table). The dominant TE tools instead reach for EM or fractional assignment (TEtranscripts/Telescope/SQuIRE use EM; TE-Seq drives Telescope fractionally); Random-One is our integer-compatible choice, chosen so the counts drop straight into DESeq2/edgeR without rounding.
*   **The gene and TE passes use different kernels.** The gene pass is unique-only (no `-M`), since a multimapper is genuinely unassignable to one gene; the TE passes are `-M` integer Random-One. These are different counting kernels, so the gene size factors are not the TE library scale — the basis of the gene-vs-TE caveat in §Combined and the `FLAG-KERNEL-MISMATCH` flag in [`QC.md`](./QC.md).
*   **Each fragment contributes exactly one count.** By construction under `--outSAMmultNmax 1`, each fragment emits exactly one locus while `NH` records its true multiplicity; direct inspection of the alignments confirms no fragment emits more than one locus. The delivered sense (`-s 2`) integer matrix is therefore fragment-weighted: young high-copy families are not alignment-inflated. An alignment-weighted run could instead have inflated the youngest families several-fold (in proportion to their mean multiplicity) — a risk averted by construction, not a realized effect. See [`QC.md` §9](./QC.md) (kernel tracer, fragment-weighting argument, counterfactual-risk pattern) and `FLAG-ALIGNMENT-WEIGHTED`.
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
*   **Cons**: often slower; might have limitations to the design complexity, 
e.g. TEtranscripts internally allows for only simple pair-wise comparisons of type A vs B.
(At least, I was only able to support tricker design questions by having separate runs of A vs B, A vs C; and then
joining resulting matrices by the A key, but this beats the whole idea of EM in a single system, to the best of my understanding).

**FeatureCounts and unsorted BAMs:**
STAR’s `--outSAMtype BAM Unsorted` is sufficient. FeatureCounts accepts unsorted BAMs and pairs reads internally, 
but can also work with sorted data, unlike TEtranscripts which is better suited for unsorted outputs [[2]].

> **Runnable, env-locked counting workflow:** the post-nf-core
> TE+gene featureCounts counting workflow — the two-pass `runFeatureCounts_TE_and_genes.sh`
> driver, the BAM symlink / identical-path Docker staging recipe, and the QC gate — is
> packaged as the **`te-gene-featurecounts`** skill, which ships its own pinned
> container **`te-fc:2.0.2`** (featureCounts **v2.0.2**). For reproducibility, run featureCounts
> from a container that ships subread: the locked `te-fc:2.0.2` (and the legacy
> `scdock-r-dev:v0.2`) carry it, while `scdock-r-dev:v0.5.x` ships only MultiQC's parser. The
> toolkit's own `scripts/runFeatureCounts*.sh` remain the upstream source of those drivers.


### Strandedness: when `-s 0` vs stranded — the field is split

Strand choice is a separate axis from the multimapper choice above. 
There is no benchmarked "field standard" for TE strandedness — the tools genuinely split on their defaults, 
and this section presents the two poles with their evidence grades rather than mandating one.

- **The split:** TEtranscripts/TEcount defaults to **unstranded** (`--stranded "no"`), as do SQuIRE (`-s 0`) and Telescope, 
while atena defaults to **stranded** (`strandMode = 1L`) and TE-Seq **auto-detects** per dataset (RSeQC `infer_experiment`, then propagated to genes and Telescope plus a per-TE antisense feature). 
Modern pipelines (TE-Seq, Mobile DNA 2025) recommend stranded for a directional library so 
TE-derived expression can be resolved from gene expression. 

The defaults stay genuinely split, with no settled standard.
- **For a stranded library + joint gene+TE matrix, the more principled route seems to be stranded TEs matched to the genes** 
(`-s 2` for a reverse dUTP library) — matching tool, strandedness, and multimapper policy across genes and TEs keeps them 
on a comparable basis. 
**(grade B — field split; see grading table).** 
The one study simulating both modes reports stranded FDR (54.9%) ≤ unstranded (58.7%) (Savytska 2022), so the limited direct evidence is consistent with the best-practice choice — though it is FDR-only, single-study, simulation.

- **`-s 0` (unstranded) is a defensible standalone choice** for:
(a) **standalone TE-family quantification** (it matches the dominant tool's default) 
or (b) a **genuinely non-directional library**.

**The "unstranded → better TE sensitivity" claim remains unbenchmarked (grade GAP — see grading table).**

The only both-mode study (Savytska 2022) found stranded FDR ≤ unstranded, 
so unstranded reads as a **sensitivity-for-specificity trade, quantified only on the specificity side**. 
Treat `-s 0` as a defensible standalone option for the two cases above, rather than a universal requirement.

- **Bidirectional / antisense transcription TE biology is class-specific.** 
Established for **L1 antisense promoter (ASP) / ORF0 and LTR/ERV** 
(Criscione 2016, BMC Genomics doi:10.1186/s12864-016-2800-5; Faulkner 2009, Nat Genet doi:10.1038/ng.368); 
weak/passive for SINE (Alu passengers) and intronic passengers 
(grade C — class-specific mechanism; see grading table). 
The commonly cited "70%+ of TE loci are bidirectional" figure has no primary source — it conflates two unrelated genome-scale facts (≈66–69% of the genome is *TE-derived*, de Koning 2011; ≈80% of the genome is *transcribed*, ENCODE 2012) and is not a statement about TE loci 
(grade D; full trace in [`BIOLOGY.md`](./BIOLOGY.md)); leave it out. 

The information-preserving way to keep genuine antisense signal is a **stranded sense/antisense split** — count TEs `-s 1` and `-s 2` as two separate feature sets, as TE-Seq does with its explicit per-TE antisense feature (`__AS`), with SQuIRE as conceptual precedent for per-locus strand resolution **(grade B — field split; see grading table)**. Keeping the split preserves the strand that separates autonomous TE transcription from host read-through; collapsing to `-s 0` discards it and structurally inflates TE counts.
- **The Teissandier 2019 / Mobile DNA benchmark (ref [[1]]) does not bear on this choice:** it benchmarks **multimapper handling** (unique vs random-one vs fractional), and stops there. 
Cite it for Random-One/fractional accuracy only (grade A — see grading table).
- **Empirically, strand choice is nearly count-neutral for genes.** 
On an internal reverse-stranded dataset the gene-level ratio `s0/s2 ≈ 0.95` (an internal measurement, not a literature value): 
unstranded assigns ~5% *fewer* gene reads, because it loses antisense-overlapping gene pairs to ambiguity that `-s 2` disambiguates, while `s1/s2 ≈ 0.10` (forward = noise). 
The "unstranded counts much more" intuition therefore applies to TEs but not genes — the both-strand magnitude question lives on the **TE** side (genuine both-strand transcription), which is exactly why a mismatched gene `-s 2` / TE `-s 0` pairing breaks cross-feature magnitude comparability.

---

## Why featureCounts, and why the toolkit stops at the count matrix (design freedom)

The reason to use STAR Random-One + featureCounts here is less about counting TEs *better* than the EM tools and more about what it lets you do *after* counting. featureCounts emits a clean **integer** matrix and nothing else — the toolkit hands you that matrix and stops, so **you own the statistical model.**

That matters because several TE tools **bundle** their differential expression and lock it to a simple two-group comparison. TEtranscripts and SQuIRE run a fixed treatment-vs-control `~ condition` DESeq2 step; even TE-Seq, which is more flexible, fits an additive `~ [batch] + condition` only — no interactions, no factorial cells, no random effects, no custom contrasts. To ask a richer question through those tools you end up running several pairwise jobs and stitching their count tables together by row id — a workaround, not a design.

Handing back a raw integer matrix removes that ceiling. The counts drop straight into edgeR / limma / DESeq2 with **any** design you can write — 2×2 factorial, interactions, blocking/batch, continuous covariates, random effects (`duplicateCorrelation` / `dream`) — and any `makeContrasts` you need. This design freedom is the toolkit's real contribution. In fairness it is not *unique* to featureCounts: atena's `SummarizedExperiment`, Telescope's counts, TEtranscripts' own `TEcount`, and SQuIRE's `--table_only` all emit matrices you can model yourself. What this toolkit adds is a **clean integer** matrix — no rounding of fractional or EM weights — that feeds your model natively, paired with the strand sense/antisense channel and the QC gate.

The deliberate trade-off: Random-One collapses to a robust **subfamily** read-out and discards locus-level information, so locus questions need reprocessing (Rung 3 — see [`LIMITATIONS.md`](./LIMITATIONS.md)). Design freedom and the locus discard are two sides of the same minimal-by-choice instrument.

---

## 3. Analysis Architecture: Combined vs. Separate

Should genes and TEs be analyzed in the same matrix?

### The Case for Combined Analysis (Recommended)
Literature suggests a combined approach can be beneficial [[1]].

### 1. Shared Normalization and Dispersion
*   **Normalization (anchor size factors on genes) (grade B — contested; see grading table):** Size factors are best estimated from the **genes alone** — DESeq2 `estimateSizeFactors(dds, controlGenes = isGene)` (median-of-ratios on gene rows), or **equivalently in the edgeR/limma idiom** `calcNormFactors(method = "TMM")` computed on a genes-only `DGEList` with its `norm.factors` transplanted onto the combined object (`filterByExpr` applied *per feature type*, then `voom`/`voomLmFit` → `lmFit` → `eBayes`) — rather than pooled across genes+TEs. Same gene anchor, different estimator (median-of-ratios vs TMM). TE-Seq makes exactly this choice, with its verbatim `estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_cts$gene_id)` and the comment "estimate the size factors using genes, and not RTEs" (`deseq.R:163-164`); the anchor is contested because TEtranscripts and SQuIRE instead pool genes and TEs into one DESeq estimate. The long-tailed, multimapper-inflated TE minority violates the "most features are unchanged" assumption that TMM/RLE rely on, and a pooled estimate can drag the gene log-fold-changes. Pooling genes+TEs into one size-factor estimate is safe only when the TE signal is modest. The benefit of combined mode is having the TEs *present in the same object* (so the design and dispersion are shared), which is separate from putting them into the *size-factor estimate*. **Kernel caveat:** the gene matrix is unique-only (no `-M`) while the TE matrices are `-M` integer Random-One. These are different counting kernels. Gene-derived size factors are correct **for gene DE** (multimappers are genuinely unassignable to one gene), but the TE matrix's library scale includes `-M` TE signal the gene factors never saw ⇒ **gene size factors are not the TE library scale**; keep them off the TE matrix unless the mismatch is documented. Normalize TEs on a TE-appropriate (or spike/total) basis, or at minimum record the mismatch. See the QC flag `FLAG-KERNEL-MISMATCH` in [`QC.md`](./QC.md).
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
*   *This is necessary but not sufficient.* Removing double-counting does **not** by itself make the matrix fully comparable across feature types. 
If genes and TEs were counted on different bases (e.g. gene `-s 2` unique-only vs TE `-s 0` `-M`), the gene and TE rows sit on different measurement bases. 
The combined matrix is then valid for **within-feature-type, across-sample differential expression only**, and **only with size factors estimated from genes alone** (see Normalization above). 
It supports within-feature-type, across-sample differential expression; gene-vs-TE magnitude comparison within a sample falls outside it (the no-cross-feature-magnitude rule is **grade C — inference**, 
not a cited TE-primary standard), and **"TE % of transcriptome" reads as a QC sanity band, not as biology** — treat it accordingly. 
(The more principled fix, best-practice grade B and not proven-superior, is to recount TEs at the gene strandedness, `-s 2`, so both feature types share one basis.)
Beyond strandedness, the row-bound joint matrix also **mixes two counting kernels**: the TE rows are `-M` (a multimapper contributes to *all* its loci) and the gene rows are unique-only (a multimapper contributes to *none*). A TE count and a gene count in one column therefore sit on different per-fragment bases ⇒ **any gene-vs-TE ratio or "fraction of library in TEs vs genes" is kernel-inconsistent — flag, do not report (grade A — kernel provenance; see grading table).** (See [`QC.md`](./QC.md) flag `FLAG-KERNEL-MISMATCH`.)

## Gene vs TE overlapping reads
Reads mapping to TEs within genes (introns/UTRs) can be double-counted.
*   **Combined Analysis:** A read overlapping a gene exon and a TE region could be counted for both if not handled.
*   **Solution:** Remove TE SAF regions that overlap gene exons (as done in this toolkit's recommended workflow) or use tools like TEtranscripts that prioritize gene assignment [[3]].
Most TE tools, like TEtranscripts, will not count a read toward a TE if it overlaps a gene exon, to avoid attributing host gene expression to TEs.

**Exon subtraction is double-counting hygiene, not autonomy isolation.** Removing exon-overlapping TEs keeps host *exon* (mature-mRNA) expression from being miscounted as TE; it says nothing about read-through, which is mostly *intronic*. Intronic TEs are not exonic, so they remain in the SAF and still catch host read-through. Isolating the autonomous-candidate signal — intergenic vs intronic — is genic-context stratification (Rung 2 in [`LIMITATIONS.md`](./LIMITATIONS.md)), which the field handles by *labelling* TEs by context (e.g. TE-Seq's Exonic/Intronic/Intergenic), not by carving the annotation. (grade A — TE-Seq stratifies on genic context.) 

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
| **Separate Analysis** | High multi-map limits. Separate runs (Genes: unique; TEs: fractional). | Separate TMM/DE runs. | Simple implementation. | Risk of double-counting. Normalization discordance. Separate FDRs. Per-feature-type strandedness is a *tolerated interim trade-off* that forfeits cross-feature comparability. |
| **Combined Analysis** | Combined annotation (non-overlapping). | Gene-anchored size factors (DESeq2 `controlGenes` / edgeR genes-only TMM). Single DE model. | Unified FDR. No double-counting. Robust dispersion. *(Valid for within-feature-type DE; not gene-vs-TE magnitude.)* | Complexity in prep. Requires checking overlaps; size factors must be gene-anchored, not pooled. |
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

## Field trajectory / long reads

Short-read **subfamily-level** TE quantification (what this toolkit does) remains current and routine in peer review, and new tools still baseline against TEtranscripts/featureCounts (grade A at subfamily level — see grading table). Long-read sequencing (ONT/PacBio) **complements** short reads: long reads own **locus identity, isoform structure, and chimera resolution**, while short reads carry the depth for sensitive differential expression — even the flagship hybrid **LocusMasterTE (Genome Biol 2025, doi:10.1186/s13059-025-03522-9) injects long-read TPM into a short-read EM** precisely because long-read depth is too low for DE. The genuinely **contested** part is short-read **locus-level** quantification (high FDR — best tool TElocal still ~26% FDR in Savytska 2022); subfamily-level is settled. Peer-reviewed work treats long-read as complementary to short-read for TE DE, not a replacement.

---

## Evidence grading & open questions

Every load-bearing claim in this toolkit's documentation is graded so that convention, inference, and benchmarked fact are not confused with one another.

**The A/B/C/D/GAP scale**
- **A** — peer-reviewed standard (multiple reviewed papers/tools, presented as a standard).
- **B** — established tool default OR a single strong pipeline's recommendation (real but narrower basis).
- **C** — plausible mechanistic inference; sound reasoning, not stated as such in a TE primary source.
- **D** — folklore / asserted without primary support / mis-imported.
- **GAP** — no adequate primary source found; the question is open / unmeasured.

**Load-bearing claims with grades + primary citation**

| Claim | Grade | Primary source / basis |
|---|---|---|
| Multimapper kernel: random-one integer ≈ `-M --fraction`; unique-only undercounts (young families severely) | **A** | Teissandier 2019, Mobile DNA, doi:10.1186/s13100-019-0192-1 (cited by TE-Seq `README:79`); dominant tools use EM/fractional instead (TEtranscripts SQUAREM, Telescope/SQuIRE EM, TE-Seq Telescope `average`) |
| `-M` required under STAR Random-One (featureCounts reads `NH`, discards multimappers without it) | **A** | tool semantics; confirmed by direct inspection of the alignments |
| Joint gene+TE matrix is a reviewed construct | **A** | TEtranscripts; SQuIRE (joint, split post-hoc); TE-Seq (separate-then-combine: featureCounts genes + Telescope TEs) |
| DESeq2 needs integers; limma-voom handles fractional counts natively | **A** | tool semantics |
| No TPM/FPKM for TE meta-features; no gene-vs-TE within-sample magnitude comparison | **C — inference** | mechanistic (summed loci have no single length); not stated in any TE primary source |
| Stranded TEs matched to genes for a joint stranded matrix | **B — best-practice, field split** | atena default `strandMode=1L` (stranded); TE-Seq auto-detects per dataset (RSeQC `infer_experiment`); TEtranscripts `--stranded` defaults `"no"`, SQuIRE `-s 0` and Telescope unstranded |
| Sense/antisense split to preserve bidirectional biology | **B — field split** | TE-Seq explicit per-TE antisense feature (`__AS`); SQuIRE resolves per-locus strand (conceptual precedent, ships no labeled sense/antisense matrix) |
| Genes-only (`controlGenes`) size factors | **B — contested** | TE-Seq *for*: `estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_cts$gene_id)`, "estimate the size factors using genes, and not RTEs" (`deseq.R:163-164`); TEtranscripts and SQuIRE pool genes+TEs (against) |
| Genic-context stratification (Exonic/Intronic/Intergenic, gene-adjacent, UTR) | **A** | TE-Seq labels every TE locus by genic context and uses it as a DE stratifier (`annotate_rtes.R:882-1001`, `aggregate_count_statistics.R:96`) |
| Matched strandedness for genes+TEs in single-pass tools | **B — mechanical default** | one global `--stranded` flag (TEtranscripts/atena) |
| "Unstranded → better TE sensitivity" | **GAP — unbenchmarked** | never measured; Savytska 2022 found stranded FDR ≤ unstranded |
| TE bidirectionality is class-specific (L1-ASP/ORF0, LTR/ERV real; SINE/intronic weak) | **C** | Criscione 2016; Faulkner 2009 |
| "70%+ of TE loci bidirectional" | **D — mis-imported** | conflates genome TE-composition (de Koning 2011) with genome transcription (ENCODE 2012); not a TE-loci statistic — see BIOLOGY.md |
| "featureCounts/TEtranscripts is still the primary method" | **B — inferred** | continued routine use; explicit "most widely used" wording only in a 2026 preprint |

### Open gaps / convention-not-evidence

These are where the field runs on convention, simulation, or inference rather than settled evidence. Do not paper over them.

1. **No ground-truth gold standard.** Every TE benchmark rests on **simulation** (Polyester etc.); no experimental dataset has known per-locus TE transcript counts. The 2025 comprehensive benchmark says best practice is "currently lacking owing to the difficulties of evaluation."
2. **TE strandedness is under-benchmarked.** Exactly **one** study (Savytska 2022) simulated both modes, and it reports **FDR only**. The "unstranded → better sensitivity" claim has **never** been measured (**GAP**).
3. **Gene–TE disambiguation is unsolved.** Intron retention, exonized fragments, and read-through inflate TE counts; tools diverge (SQuIRE/TElocal prioritize genic reads, TEspeX masks exonized fragments, TE-Seq only diagnoses). No consensus fix — "a major unresolved challenge."
4. **Young / identical-copy loci are intrinsically unresolvable at locus level by short reads** — an information-theoretic limit, not a tooling gap.
5. **The no-cross-feature-magnitude / no-TPM-for-meta-features rule is inference (C)**, not a cited TE-primary standard.
6. **"featureCounts/TEtranscripts is still primary" is INFERRED** from continued routine use; the explicit "most widely used" wording appears only in a 2026 preprint (MAJEC), not a peer-reviewed verbatim claim.

---

## 7. References

1. **Mobile DNA (2019)**: [Tools and best practices for retrotransposon analysis](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1)
2. **Biostars**: [Can featureCounts use sortedbycoordinate files?](https://www.biostars.org/p/488290/)
3. **PMC (2015)**: [TEtranscripts: a package for including transposable elements](https://pmc.ncbi.nlm.nih.gov/articles/PMC4757950/)
4. **J Transl Med (2021)**: [TPM, FPKM, or Normalized Counts?](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w)

---

## Changelog

### 2.0.2 — 2026-06-15 (kernel/normalization caveat + QC cross-ref)
Documentation-only pass; no changes to `R/` or `scripts/` code. Added the gene-vs-TE kernel caveat and grounded the integer Random-One choice after the strand-split overlap-ambiguity reconciliation:
- Normalization: added the gene (unique-only, no `-M`) vs TE (`-M` integer Random-One) kernel caveat — gene size factors are not the TE library scale; keep them off the TE matrix unless the mismatch is documented.
- Combined matrix: noted that the row-bound joint matrix mixes two kernels (TE `-M` vs gene unique-only), so any gene-vs-TE ratio / "fraction of library" is kernel-inconsistent (flag, do not report).
- Strategy A: stated the counting is integer Random-One (no `--fraction`) and grounded the gene-vs-TE kernel difference (pointing to the existing Teissandier line).
- Cross-referenced the new [`docs/QC.md`](./QC.md) (strand-split invariant, directional meter, alignment witness, working principles, warning-flag taxonomy, GREEN/RED gate) and its `FLAG-KERNEL-MISMATCH` flag.

### 2.0.1 — 2026-06-11 (documentation correction)
Documentation-only pass; no changes to `R/` or `scripts/` code. Removed two myths and added missing caveats after a claim-by-claim reconciliation against the assembled evidence:
- Strandedness: corrected the "TEs must be counted unstranded (`-s 0`) because they are bidirectional" claim; added a dedicated "Strandedness: when `-s 0` vs stranded" subsection. For a stranded library + joint matrix, count TEs at the gene strandedness (`-s 2` here); preserve bidirectional biology via a sense/antisense split.
- Removed the implied Teissandier 2019 / Mobile DNA endorsement of `-s 0` — that benchmark concerns multimapper handling only (still cited for Random-One/fractional accuracy).
- Softened the "lose ~80–90% of TE signal" figure to a young-family-biased fraction (~10–21% overall multimapper fractions; higher for young families); direction unchanged.
- Combined matrix: clarified exon-subtraction is necessary but not sufficient; valid for within-feature-type DE only, with gene-anchored size factors; gene-vs-TE magnitude and "TE %" not interpretable.
- Normalization: size factors must be anchored on genes, not pooled across genes+TEs.
- Reframed "optimizes strandedness per feature type" as a tolerated interim trade-off, not an advantage.

**Transparency / evidence-grading refinement (same 2.0.1 pass, no version bump):**
- Corrected the over-claim that "stranded is the field standard." The field is **SPLIT**: the dominant tool TEtranscripts **defaults to UNSTRANDED** (`--stranded no`); TE-Seq/Mobile DNA 2025 recommend stranded. Stranded-for-joint is now **best-practice / mechanistically motivated (grade B), not a benchmarked standard**.
- Flagged "unstranded → better TE sensitivity" as **unbenchmarked (GAP)** — the only both-mode study (Savytska 2022) found stranded FDR (54.9%) ≤ unstranded (58.7%); it is a sensitivity-for-specificity trade, never a quantified gain.
- Removed the "70%+ TE loci bidirectional" framing (mis-imported genome-wide antisense figure, grade D); TE bidirectionality is class-specific (L1-ASP/ORF0, LTR/ERV real; SINE/intronic weak).
- Added inline **evidence grades (A/B/C/D/GAP)** to load-bearing claims (genes-only size factors = B/contested; sense/antisense split = B/SQuIRE-specific; no-cross-feature-magnitude / no-TPM = C/inference; matched strandedness = B; multimapper kernel/Teissandier = A).
- Added a new **"Evidence grading & open questions"** section (scale, graded claim table, "Open gaps / convention-not-evidence" subsection) and a **"Field trajectory / long reads"** note (short-read subfamily quant current; long-read complements not replaces; locus-level is the contested high-FDR part).
- Stated the in-house `s0/s2 ≈ 0.95` gene result is an **empirical in-house measurement (not literature)**, and that the stranded + sense/antisense split + genes-only size-factor route is the best-practice (not proven-superior) joint option, not a mandate.
