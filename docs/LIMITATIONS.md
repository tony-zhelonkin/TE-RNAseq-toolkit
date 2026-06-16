# Limitations: What This Instrument Can and Cannot Claim

This is the **boundary document**. 
It names the inferential limits of the TE quantification this toolkit performs and where a result stops being supported by the current design. 
For the *why* behind the design choices (strandedness, multimapper kernel, joint matrix, normalization), 
and for the evidence grades and citations those choices rest on, see [`METHODOLOGY.md`](./METHODOLOGY.md). 

Grades below use that document's **A/B/C/D/GAP** scale and its citations are not restated here.

## The instrument, in one line

featureCounts on a grouped, exon-subtracted SAF where 
* `GeneID = Subfamily:Family:Class` (1,243 mouse subfamilies), 
* STAR Random-One + integer `-M` (no `--fraction`), 
* counted in three strand channels — sense `-s 2` (gene-matched), antisense `-s 1` (separate channel), unstranded `-s 0` (comparison). 

The joint gene+TE matrix row-binds genes (`-s 2`) and TE-sense (`-s 2`) as mutually exclusive features, 
with **no normalization applied at this stage — the delivered matrix is raw integer counts.**
The *recommended* downstream normalization anchors per-sample scaling on the **gene rows only** —
DESeq2 `estimateSizeFactors(dds, controlGenes = isGene)`, or equivalently edgeR `calcNormFactors(method = "TMM")`
on a genes-only `DGEList` with its factors transplanted onto the combined object (`filterByExpr` per feature type) —
so a genuine en-masse TE shift cannot be absorbed into the size factor.

**What this buys, and only this:** SUBFAMILY resolution. 
The grouped SAF carries **no genic context** (intergenic / intronic / adjacent) and **no locus identity** (which copy at which address). 
Those axes were deliberately thrown away at SAF-build time. Every limitation below is downstream of the genic context throwaway.

---

## 1. Matching the strand basis is NECESSARY, not SUFFICIENT

Counting TEs at the gene strandedness (`-s 2` sense, matched to the `-s 2` gene denominator) buys two **distinct** things 
that happened to ride on the same `-s 2` flag. It's important not to conflate them.

### (a) MEASUREMENT axis — a valid ruler for normalization. **SURVIVES.**

A per-sample size factor `s_j` says "sample *j*'s gene counts are inflated/deflated by this much, in the units genes are measured in." 
Dividing a TE row by `s_j` **asserts** the TE row is measured in those same units — that genes and TEs accrue counts from true library depth 
at the same per-sample rate. Matching the strand basis seems, at least to me, internally defensible. 

A stranded derepression study should be careful to **not** normalize an unstranded TE numerator against a stranded gene denominator. 
The logic is: are you even using the right *measurement ruler*? The way the size factor was computed determines if you can apply the ruler to the TE rows.
**[evidence: B — matched-basis best-practice; see METHODOLOGY §Strandedness and §Combined.]**

### (b) BIOLOGY axis — the autonomous-vs-read-through verdict. **DOWNGRADED.**

At subfamily resolution the antisense split gives strand **asymmetry**, not an autonomy verdict. 
The `-s 1` channel is a strand-of-read aggregate over a whole subfamily.
It may become a verdict on the autonomous vs read-through TE element only when paired with: 
- the genic-context 
- and locus axes...
... the grouped SAF method I usually employ for subfamily aggregate resolution deliberately discards (§5). 
**[evidence: C — inference; the verdict-grade form is TE-Seq's three-axis design, not this repo.]**

The whole point: matching the strand basis earns you (a) outright and only a *signal* (asymmetry) on (b) — never the *verdict*.

---

## 2. What matching REMOVES vs what it leaves (the residue)

**Removes:**
- antisense-oriented read-through (host transcription on the opposite strand) — it lands in the `-s 1` channel, out of the `-s 2` numerator;
- the normalization-rate mismatch — the ruler problem of §1(a).

**Does NOT remove — the residue:**
- **co-oriented (sense) read-through.** A TE sitting sense within a sense intron dumps host transcription straight into the `-s 2` channel. **No strand flag separates that at subfamily level** — the read is on the "right" strand. This is the locus/context problem, and the current grouped `Subfamily:Family:Class` SAF does not touch it. **[evidence: C — inference; consistent with the unsolved gene–TE disambiguation gap, METHODOLOGY §Open gaps #3.]**

Matching the strand basis is a partial, not total, read-through defense. The part it cannot reach is exactly the part that needs Rung 2 (§5).

---

## 2a. Same-strand overlap silent loss (strand-split overlap-ambiguity)

This is a **NEW** limitation, **distinct from and alongside** the co-oriented-readthrough residue of §2 — both unsolved at subfamily resolution, neither subsuming the other. **Silent loss is an annotation-geometry tax** (same-strand nested overlap dropped to Ambiguity); **co-oriented read-through is a transcription-context confound** (host transcription on the right strand). They are different failure modes; closing one does not close the other.

**What it is.** A read overlapping two **same-strand** features is dropped as Ambiguity by `-s 0` *and* by both stranded passes (the strand-matched pass still sees ≥2 candidates → Ambiguity; the opposite pass sees none → NoFeatures). It contributes 0 to the unstranded, sense, and antisense counts alike. The toolkit's three-channel counting never sees it.

**The witnessed ledger.** Per-read `-R CORE` audit of the s0-ambiguous pile (sample 0019, reproduced on 0006): **AP 25.7% : M 6.1% : P(silent) 68.2%**; library-scale silent loss ≈ **9–10% of assigned TE signal**; antiparallel double-presence ≈ 4%; genuine reclaim ≈ 1%; bilateral negligible at **0.38%** (so the single-strand-drop assumption holds). Silent loss is the **dominant** fate of the ambiguous pile. **[evidence: A — per-read witness, reproduced on a 2nd sample; matches the synthetic truth table.]**

**Conditional losslessness (A1).** A prior reading framed the strand split as "lossless / recovery-dominated" — that the `s0 − sense − anti` residual being one-sided proves no loss. **That framing is refuted.** The residual is lossless **only relative to s0's *assigned* set**; same-strand silent loss is **residual-invisible** by construction (it is 0 in all three channels) and carries **no GeneID in CORE** (target = NA). A zero or one-sided residual certifies *nothing* about loss. The visible recovery the residual *can* see is the antiparallel and asymmetric slice only; it sits on top of an unseen silent-loss pile roughly twice its size.

**The operational axioms live in [`QC.md`](./QC.md).** Axioms A1–A6, the warning-flag taxonomy, the directional meter, and the gate thresholds are owned there; they are not duplicated here. This section names the limitation and its boundary. The four cases this becomes **DANGEROUS** in — each by flag name in [`QC.md`](./QC.md) — are:
- **summing channels** — `(sense + anti)` double-counts the antiparallel + A/A/A piles (`FLAG-SUM-CHANNELS`, A2); not bidirectional signal;
- **s0-denominating the ERV/satellite/DNA tail** — s0 drops their antiparallel reads to Ambiguity, a pathological denominator (`FLAG-S0-DENOM-ERV`, A3); denominate on the sense channel;
- **applying gene size factors to the `-M` TE matrix** / comparing gene-vs-TE magnitude — different kernels (`FLAG-KERNEL-MISMATCH`, A5; see §6, METHODOLOGY §Combined);
- **silent loss aligning with the design axis** — sample-stable loss cancels in cross-sample DE, but P_frac co-varies with net strand-offset (r = −0.50), so it must be checked against the design matrix before declaring condition-independence (`FLAG-SILENTLOSS-DESIGN`, A6).

**Scope — not alarmist.** The young autonomous candidates (L1Md_T/Gf/A, IAPEz) are **🟢 GREEN: 96.65% conserved** (non-conserved 3.35% < 5% gate; IAPEz-int 99.8%). The burden sits on the **old SINE / MaLR-LTR** subfamilies (silent loss) and on **ERVK-LTR / satellite** (antiparallel double-presence), not on the young autonomous tail. See [`QC.md` §7](./QC.md) for the GREEN/RED gate. The young-silent attribution is **now witnessed for the gate set** via the `-O` target-revealer ([`QC.md` §4a/§9](./QC.md)): the young silent-loss share is **1.27–1.37%**, concordant with the geometry prediction, so the GREEN gate's **second leg is closed at read level** (no longer geometry-inference alone) — the genic-context Rung-2 GAP (intron/intergenic stratification on the richer SAF, §5) stays open. **[evidence: A — `-O` per-read young-silent witness, concordant with geometry; the genic-context attribution remains the Rung-2 GAP, §5.]**

**The YELLOW DE-confound (A6).** On the SD criterion silent loss is *not* a confound (P_frac-of-signal SD 0.0028 < net-offset SD 0.0081), so it cancels in cross-sample contrasts under the usual conditions. **BUT** P_frac correlates with per-sample net strand-offset (**r = −0.50**): silent loss is not a flat tax, it co-varies with sample-level strand behaviour. This **must** be checked against the experimental design matrix before treating silent loss as condition-independent. The check is **pick-up-ready** in the dataset-side `featurecounts_results/DE_precheck/` (the design matrix is not currently loaded). **[evidence: B — sample-stable on SD; YELLOW pending the design-matrix check.]**

---

## 3. The worked case: how `-s 0` fabricates derepression

One TE subfamily. No real change in autonomous transcription. Ctrl vs Senescent, equal depth, so gene size factors `s_j ≈ 1.0` in both.

| Reads accrued | Ctrl | Sen |
|---|---|---|
| own-strand (autonomous, **unchanged**) | 100 | 100 |
| antisense read-through (**↑ in senescence**) | 20 | 60 |

- **`-s 2`** (matched, sense only): `100 → 100`. Normalized `log2FC = log2(100/100) = 0`. **Truth preserved.**
- **`-s 0`** (mismatched, both strands): `120 → 160`. Normalized `log2FC = log2(160/120) = `**`+0.42`** (≈ 1.3× **fabricated** "derepression").

The gene-anchored size factor **cannot rescue** the `-s 0` case: genes never saw the read-through bump (`-t exon`, on-strand), so `s_j` stays ≈ 1 and leaves the inflation in place. Worse, the fake signal points the **same direction** as the hypothesis (senescence → TE up), so it reads as confirmation.

This only bites when **strand-capture varies across samples** — which is exactly the derepression setting: intron retention and read-through rise with senescence / stress / immune activation, correlated with the tested biology. A *constant* rate difference would cancel in a log-fold-change; it is the *covariance with condition* that fabricates signal. **[evidence: C — inference, arithmetic; the underlying disambiguation problem is METHODOLOGY §Open gaps #3, the strand-basis rationale §Strandedness.]**

---

## 4. The antisense channel must DO WORK

Holding stranded counts on disk is necessary but not sufficient. Subfamily-level featureCounts + the stranded sense/antisense design is the **right instrument** — but the discrimination that separates real derepression from read-through is a **condition × channel interaction test** (does the sense/antisense ratio shift with condition?), **not** the mere existence of the counts.

A `-s 2` sense-channel DE result that is simply *called* "derepression" without using the `-s 1` channel to test it is still exposed to every residue in §2–§3. Correct design — but the antisense channel has to do work **before any DE is called derepression**. **[evidence: B — sense/antisense split is SQuIRE-specific design, METHODOLOGY §Strandedness; the interaction test as the discriminator is grade C inference.]**

---

## 5. The 3-rung ladder (resolution vs robustness)

Each rung adds an axis. Higher rung = more specific claim, more fragile inference.

- **Rung 1 — stranded sense/antisense, grouped SAF.** *(← what this toolkit does.)*
  Derepression vs read-through at **subfamily level, statistically**, via the condition × channel interaction (§4). featureCounts + grouped SAF. **ROBUST; current-standard genre.**

- **Rung 2 — + genic-context-stratified SAF (intergenic / intronic / adjacent).**
  Read-through attribution **structurally**: split each subfamily into e.g. `L1HS:intergenic`, `L1HS:intronic`. Intergenic loci have no host to read through, so intergenic-TE derepression is the clean **autonomous-candidate** signal; intronic "derepression" *not* mirrored intergenically is a **read-through flag**. Still featureCounts, richer SAF. Genic context lives in the SAF **annotation**, not a different downstream tool.

- **Rung 3 — EM locus assignment.**
  Names the firing copy (which L1HS at which locus); per-locus context + neighbor-gene correlation. Telescope (TE-Seq wraps it) / SQuIRE. **SPECIFIC but FRAGILE** (per-locus FDR — METHODOLOGY §Field trajectory: TElocal ≈ 26% FDR, Savytska 2022).

**Judgment.** For an **aggregate-derepression** claim, finish Rung 1 properly (the interaction) and reach for **Rung 2 before Rung 3** — a richer SAF buys most of TE-Seq's read-through defense at a fraction of the cost and none of the locus-FDR fragility. Climb to Rung 3 only when the question is genuinely **"WHICH insertion,"** because that is where you trade a robust claim for a specific one.

Compressed:
- aggregate derepression is the **ROBUST** claim, not the compromise;
- genic context is a **richer SAF**, not a different tool;
- only naming the firing copy forces **Telescope / SQuIRE**.

**The three axes, and what each earns:**

| Axes used | Verdict it supports |
|---|---|
| strand alone | strand **ASYMMETRY** — a *signal* |
| strand × genic-context | autonomous vs passive — a *verdict* (TE-Seq) |
| strand × context × locus | which copy is firing — an *address* (Telescope / SQuIRE) |

TE-Seq pairs all three (strand: sense/antisense; genic context: intergenic/intronic/adjacent; locus: Telescope EM) to call autonomous-vs-read-through. **This toolkit has only the strand axis.** So it earns the signal, not the verdict and not the address. **[evidence: B — TE-Seq / SQuIRE three-axis design; the rung judgment is grade C inference.]**

---

## 6. Magnitude caveat

Do **not** promote TE% or gene-vs-TE magnitude to biology. Making the joint matrix normalization-**coherent** (mutual exclusivity via exon subtraction, genes-only size factors) is the most seductive moment to over-read it — **coherence is not magnitude.**

- No TPM / FPKM for TE meta-features (summed loci have no single length — METHODOLOGY §Open gaps #5). **[evidence: C — inference.]**
- "TE % of transcriptome" is a QC sanity band, not a biological fraction (METHODOLOGY §Normalization, §Combined).
- The joint matrix is valid for **within-feature-type, across-sample DE only**, never gene-vs-TE within-sample magnitude.

The shared strand basis earned strand **asymmetry** and a valid normalization **ruler** — **NOT** an autonomy verdict and **NOT** a magnitude statement. **[evidence: C — the no-cross-feature-magnitude / no-TPM rule, METHODOLOGY §Combined and §Open gaps #5.]**

---

## One-line boundary

This instrument measures **subfamily-level, sense-channel, depth-normalized differential abundance**, and — when the antisense channel is made to work (§4) — **strand asymmetry**. It does not, on its own, certify **autonomy**, **read-through-free** signal, **locus identity**, or **magnitude**, and it does not recover **same-strand-overlap silent-loss reads** (residual-invisible, §2a). Each of those requires climbing the ladder (§5).
