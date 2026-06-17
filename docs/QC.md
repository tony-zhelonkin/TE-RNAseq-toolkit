# Strand-split TE counting — QC, invariants, and foot-guns

This document is the **QC gate** and the **warning taxonomy** for the three-channel (unstranded / sense / antisense) TE counting this toolkit performs. 
The document covers: 
* the strand-split invariant, 
* the directional meter, 
* the regime-witness procedure, 
* the working principles, 
* the warning-flag taxonomy, 
* and the GREEN/RED gate thresholds. 

It is downstream of the subfamily-resolution fact stated in [`LIMITATIONS.md`](./LIMITATIONS.md): 
the grouped `Subfamily:Family:Class` SAF carries no genic context and no locus identity, and every boundary below sits on top of that.

For the biology behind the strand channels — sense/antisense, bidirectional transcription, read-through, derepression — see [`BIOLOGY.md`](./BIOLOGY.md).

Grades below use the **A/B/C/D/GAP** scale of [`METHODOLOGY.md`](./METHODOLOGY.md); its citations are not restated here. 
This document **defers** to [`METHODOLOGY.md`](./METHODOLOGY.md) for the multimapper-kernel rationale (Strategy A/B, Teissandier) 
and to [`LIMITATIONS.md`](./LIMITATIONS.md) for the 3-rung ladder and the co-oriented-readthrough residue — links, not duplication.

**convention (reverse-stranded dUTP library, mm39).** For this library type, the featureCounts strand flags map to:
* sense = `-s 2` assigned,
* antisense = `-s 1` assigned,
* unstranded = `-s 0` assigned.

After this block the prose uses the words **sense / antisense / unstranded**.

featureCounts v2.0.2, all TE passes `-M -F SAF -p --countReadPairs -B -C`, no `-O`, no `--fraction`, default `--minOverlap 1`, integer.

---

## 1. The strand-split invariant (conditional losslessness)

Define, per subfamily and library-wide:

```
unstranded ≈ sense + antisense + residual
residual_frac = residual / unstranded
```

**Conditional losslessness (grade A — per-read audit, reproduced on a second sample).** The residual `unstranded − sense − antisense` is lossless **only relative to the unstranded pass's *assigned* set**. A read overlapping two **same-strand** features is dropped as Ambiguity by the unstranded pass **and** by both stranded passes (the strand-matched pass still sees ≥2 candidates → Ambiguity; the opposite pass sees none → NoFeatures). It contributes 0 to all three passes, so `residual = 0`. This **silent loss** is by construction **invisible** to the residual: a zero or one-sided residual certifies **nothing** about loss. A *positive* residual would be a bug — a stranded pass dropping a read that the unstranded pass assigned — and only a vanishing proper-pair `-B` edge (the A/M/M case below) does this.

**On the residual and recovery.** It is tempting to read a balancing residual as proof that the split is "lossless" and that the net strand offset is entirely overlap *recovery*. It is not. The residual only ever sees the antiparallel and asymmetric slices it *can* see; it is structurally blind to same-strand silent loss. The visible recovery is genuine, but it sits on top of a larger unseen silent-loss pile (roughly twice its size — see below). So a one-sided residual is *consistent with* large silent loss, and the right reading is: the split is conditionally lossless relative to the unstranded assigned set, and silent loss is the dominant fate of the ambiguous pile.

---

## 2. The witnessed ledger (one table)

These figures come from an internal per-read audit of a reverse-stranded mouse dataset, on a high-depth sample and reproduced on a low-depth one. (Exact per-sample counts live in `docs/_internal/`; only the rounded conclusions are kept here.)

Each fragment is labelled by how many features the **sense** and **antisense** passes each assign it — written `(sense, antisense)`. A strand-ambiguous fragment then meets one of three fates: **AP — antiparallel** `(1,1)`, kept by *both* passes under *different* subfamilies (double-presence — what summing the channels would double-count); **M — asymmetric** `(2,1)/(1,2)`, where splitting by strand breaks the tie and the fragment is genuinely *reclaimed*; **P — parallel / same-strand** `(n,0)`, separable by no strand flag and so *silently lost* in all three passes — the dominant fate. The two numeric columns below use **different denominators** (share of the ambiguous pile vs share of total assigned signal) — do not mix them.

| Regime | % of ambiguous pile | Library-scale (of assigned signal) | Grade | Meaning |
|---|---:|---:|:---:|---|
| **AP — antiparallel (1,1)** | **~26%** | ≈ **4%** | **A** | double-presence; +1 to EACH stranded pass, **different** GeneIDs (100% cross-family, 0% same-GeneID) |
| **M — asymmetric (2,1)/(1,2)** | **~6%** | ≈ **1%** | **A** | single reclaim = **genuine recovery**; +1 to excess |
| **P — silent loss** (parallel n,0 + bilateral ≥2,≥2) | **~two-thirds** | ≈ **9–10%** | **A** | residual-invisible loss; bilateral negligible → single-strand-drop assumption **HOLDS** |

**Locked one-liner:** of the ambiguous pile, same-strand silent loss is the dominant fate (~two-thirds), antiparallel double-presence ~a quarter, asymmetric reclaim a small remainder. At library scale that is silent loss ≈ **9–10% of assigned TE signal**, double-presence ≈ 4%, genuine reclaim ≈ 1%. Silent loss is the **dominant** fate and is larger than any scalar solve implied. Bilateral overlap is negligible, so the single-strand-drop assumption holds **(grade A)**.

---

## 3. The `excess/ambiguous` directional meter

**Definition.** Over the ambiguous pile, `excess = (sense_Assigned + antisense_Assigned) − unstranded_Assigned`; the meter is `excess / ambiguous ≈ 0.8`, and it is tight and sample-stable across the dataset.

**Graded directional-only (grade C; grade A on the negative — it is not an estimator).** The meter is real and sample-stable, and it tracks the antiparallel-vs-silent balance **directionally**, but it does **not** yield the split. It is contaminated by (i) the **A/A/A leak** (unstranded-Assigned fragments that BOTH stranded passes also assign to a *different* feature — about a quarter of excess), (ii) asymmetric multi-feature events, and (iii) the sub-25bp `--minOverlap` annotation edge (see below, §8). **Use it as a meter, never as an estimator.**

**The methodological lesson (stated once).** An elegant scalar 3×3 solve over excess / ambiguous / ambsum reconciled cleanly to the integer on the identities it could see, yet was **wrong**: it placed the silent and asymmetric shares far from what the per-read witness later showed. The reason is the **A/A/A leak**, which adds +1 each to excess and is folded invisibly into the asymmetric term — the scalar solve cannot see it. Grade the `excess = 2·AP + M` identity **refuted as an estimator**; the mechanism `excess ∈ {0,+1,+2}` is upheld (see below, §4). **Lesson: for overlap accounting, witness beats inference** — a residual that balances does not certify a mechanism; it can balance because two unseen errors cancel.

**The counterfactual-risk pattern — the second cross-cutting lesson (grade A — method).** When a kernel choice could silently inflate or deflate a headline, **quantify the averted-error magnitude explicitly, then confirm the realized number is the safe one.** A pin you cannot put a downside number on is not a QC gate. The load-bearing instance is alignment- vs fragment-weighting of young high-copy families (see §9): the matrix is fragment-weighted, but a mis-configured `--outSAMmultNmax all` + `-M` run *would have* inflated young high-copy families by their mean alignment multiplicity — a large fold for the youngest, highest-copy families. Recording that counterfactual fold is what makes the kernel pin load-bearing: it states what was at stake and that the witness confirms it was averted. This pairs with "witness beats inference" as the second cross-cutting doctrine — prove both that the number is right and that the wrong number you avoided was large.

**Denominator discipline — caution (grade A).** Always print numerator and denominator at the **same scale** and **label the scale** on every ratio (per-sample vs library-summed). A per-sample numerator must never be divided by a library-summed denominator, or vice versa — the two share-of-excess and share-of-ambiguous ratios have different numerators *and* different denominators and cannot be mixed. A/A/A is **unstranded-ASSIGNED**, so `A/A/A ÷ ambiguous` is **not** a real fraction of the ambiguous pile — do not report it. → `FLAG-DENOM-SCALE` (see §10).

---

## 4. The regime-witness procedure (reproducible)

Run featureCounts **3× on the SAME BAM** (unstranded / sense / antisense), **identical kernel** (`-M -F SAF -p --countReadPairs -B -C`, no `-O`, no `--fraction`), each emitting the per-fragment assignment table. A runnable, parameterised version of this whole witness is packaged in the `te-gene-featurecounts` skill (see Resources) — prefer it over re-deriving the procedure by hand.

> **WARNING (load-bearing foot-gun).** Join the three per-fragment tables by **hash-join on fragment ID**, NOT positional `paste`. Read-ID order differs across threaded featureCounts runs; a positional paste silently misaligns fragments and fabricates regimes. This is the single easiest way to corrupt the witness.

Classify each fragment's (unstranded, sense, antisense) Assigned-count status triple into a regime:

| Regime | Triple | Excess | Meaning |
|---|---|:--:|---|
| conserved | (1,0) / (0,1) | 0 | one strand, one feature |
| antiparallel (AP) | (1,1) | +2 | double-presence, different GeneIDs |
| asymmetric (M) | (2,1) / (1,2) | +1 | single reclaim = **genuine recovery** |
| parallel-silent (P) | (n,0) | 0 | **SILENT LOSS** (residual-invisible) |
| bilateral | (≥2,≥2) | 0 | silent loss (negligible) |
| **A/A/A leak** | unstranded Assigned, both stranded passes Assign a *different* feature | +1 each | the bucket the scalar solve folds invisibly into M — **witnessed benign**, see below |

Reconstruct `excess` from the regimes and apply the gate.

**What the A/A/A bucket actually is — the per-read audit resolves it (grade A).** It is tempting to call A/A/A "the arithmetic cross-channel leak," but the per-read witness shows it is **benign strand-resolved antiparallel-feature selection on single-locus reads**. Of the fragments whose sense and antisense GeneIDs *diverge*, the overwhelming majority are uniquely mapped (NH==1) — per-alignment double counting is categorically impossible for a uniquely-mapped fragment — and the unstranded GeneID equals **exactly one** of {sense, antisense} for **all** of them (never both, never neither). Such a fragment overlaps **two annotated TE features on opposite genomic strands** (nested/abutting antiparallel TEs); featureCounts derives fragment strand from R1, so the sense pass legitimately picks the +strand feature and the antisense pass the −strand feature — same locus, one alignment, two strand-consistent features. So a single unstranded GeneID is fully consistent with the two stranded passes legitimately picking the two opposite-strand features; this is exactly what Random-One produces, not evidence against it. The divergence inflates only `sense + antisense` — a sum you must never form (`FLAG-SUM-CHANNELS`) — and never the SENSE matrix itself (each A/A/A fragment is +1 to exactly one sense GeneID). See §9 for the kernel-tracer and the fragment-weighting theorem this rests on.

**Mechanism gate: per-fragment `excess ∈ {0,+1,+2}`.** Violation threshold `< 0.1%`; observed well under that on both audited samples **(grade A)**. Agrees with the synthetic truth table (the AP/M/P cases give +2/+1/0). The lone violation cell (A/M/M, excess −1) is a `-B`/proper-pair overlap edge and is negligible.

This procedure and its truth table are regression-tested in the `te-gene-featurecounts` skill, `tests/strand_qc/` (see Resources).

### 4a. Attributing the silent pile with `-O` (the target-revealer procedure)

**The problem.** Under the default kernel (no `-O`) featureCounts writes the silent pile as `Unassigned_Ambiguity  -1  NA` — the candidate target list is **not emitted**, so silent-loss fragments carry **no GeneID** (the `UNATTRIB_silent` limitation flagged in §7). To attribute the pile by subfamily you must recover the features each silent fragment was ambiguous **among**.

**The fix — a kernel-identical re-run plus one flag (grade A; every silent fragment matched a multi-target record).** Re-run the unstranded pass **kernel-identical except add `-O`** (`-M -O -F SAF -s 0 -p --countReadPairs -B -C -R CORE`). The alignment total reproduces the audited run exactly; every fragment that was `Unassigned_Ambiguity` in the default run becomes `Assigned n_targets≥2` with the comma-separated GeneID list = the candidate features it overlapped. Then look up each silent fragment's target list by fragment ID.

**`-O` is a target-revealer ONLY — it does NOT redefine the silent set.** The silent set stays fixed by the **default-kernel** unstranded/sense/antisense statuses (the P bucket of §2); `-O` only attaches a candidate-target list to fragments already classified silent. Attribution policy: **contains-young** (≥1 young target) = the upper-bound witness; **exclusively-young** (all targets young) = the floor.

> **Foot-gun — never collide the two `-O` regimes.** The production count matrices correctly use **no `-O`** (with `-O` a fragment would double-assign across all its targets and inflate counts). The §4a use of `-O` is the **diagnostic** one, on a **throwaway re-run**, never on the count matrix. Keep the two uses strictly apart.

This converts the young-silent attribution from read-independent geometry inference to a per-read witness (§9 geometry-vs-witness concordance; the young gate-set share is small and concordant). It closes the §8(ii) limitation at subfamily resolution for the witnessed sample; the genic-context stratification (intron/intergenic) remains a Rung-2 GAP.

---

## 5. The warning-flag taxonomy (THE canonical table)

This is the **single source of truth** for the named flags. Sidecars and skill docs cite them **by name only** — never re-defining them. Severities: **🔴 hard** (violating it produces a wrong number or fabricated biology) · **🟡 watch** (conditional — must be checked). Each flag is a measurement **assertion**, never a biological conclusion.

| Flag name | Sev | One-line assertion |
|---|:--:|---|
| `FLAG-SUM-CHANNELS` | 🔴 | Never use `(sense + antisense)` as a quantity — it double-counts the antiparallel pile (+1 to each channel, 100% cross-family) **and** the A/A/A pile; the sum-over-unstranded ratio is a double-counting artifact, not bidirectional signal. |
| `FLAG-S0-DENOM-ERV` | 🔴 | Never use the unstranded channel as denominator for the LTR/ERV/Satellite/DNA subfamilies — it drops their antiparallel reads to Ambiguity (pathological denominator); denominate on the **SENSE** channel. |
| `FLAG-KERNEL-MISMATCH` | 🔴 | Gene matrix is unique-only (no `-M`, multimapper pairs dropped); every TE pass is `-M`. Do not apply gene size factors to the TE matrix unqualified; no gene-vs-TE magnitude / "fraction of library" — not like-for-like. |
| `FLAG-RESIDUAL-NOT-LOSS` | 🟡 | A zero/positive `unstranded − sense − antisense` residual does NOT certify no loss; same-strand silent loss (≈9–10% of assigned TE signal) is residual-invisible by construction (see conditional losslessness). |
| `FLAG-METER-NOT-ESTIMATOR` | 🟡 | `excess/ambiguous ≈ 0.8` is a directional meter, never an estimator of the AP/M/P split (A/A/A + asymmetric + minOverlap-edge contamination). |
| `FLAG-SILENTLOSS-DESIGN` | 🟡 | Silent loss is sample-stable (cancels under usual conditions) BUT its share co-varies with net strand-offset; verify against the design matrix before declaring condition-independence (`DE_precheck/`). |
| `FLAG-ALIGNMENT-WEIGHTED` | 🔴 | The SENSE integer matrix is fragment-weighted under Random-One (`--outSAMmultNmax 1`, matrix==nfrag); if BAMs ever emit >1 locus/fragment (or `--outSAMmultNmax` > 1 with `-M`), young high-NH families inflate up to ~meanNH-fold — re-verify with the kernel-tracer (§9) before trusting young fold-changes. |
| `FLAG-DENOM-SCALE` | 🟡 | Per-sample vs library-summed denominators must never be mixed; label the scale on every ratio. A/A/A is unstranded-ASSIGNED, so `A/A/A ÷ ambiguous` is not a real fraction of the ambiguous pile — do not report it. |
| `FLAG-READTHROUGH-SENTINEL` | 🟡 | A `-s 2` sense-up / antisense-flat "derepression" call is NOT validated by the antisense channel: co-oriented (same-strand) read-through is count-identical to autonomous and unseparable by strand. A broad antisense rise with condition is the global read-through sentinel — it raises read-through suspicion across all sense calls but is necessary-not-sufficient (null on the co-oriented residue); the separation requires the genic-context Rung-2 hand-off (intergenic vs intronic) before the word derepression is used (see [`LIMITATIONS.md` §4–§5](./LIMITATIONS.md)). |

Names are stable identifiers — do not rename downstream.

**The ratio caveat (folds into the table above).** The sum-over-unstranded ratio is **not** bidirectional signal — it double-counts the antiparallel pile **and** the A/A/A pile (`FLAG-SUM-CHANNELS`). For the LTR/ERV/Satellite/DNA subfamilies, denominate on the **SENSE channel**, not the unstranded one, because the unstranded pass drops their antiparallel reads to Ambiguity and gives a pathologically low denominator (`FLAG-S0-DENOM-ERV`). The `sense/unstranded` and `antisense/unstranded` columns in `qc_ratios.tsv` are **QC meters, not normalization denominators**.

---

## 6. The working principles (graded)

These are the principles I lean on when reading TE strand-split counts. They are **working principles, not received truth** — grounded variously in the public literature (cited in [`METHODOLOGY.md`](./METHODOLOGY.md)) and in our own internal per-read data investigations, and I expect to revise them as the evidence improves. Each carries an evidence grade on the same A/B/C/D/GAP scale, and each maps to a named warning-flag in §5. Where a principle rests mostly on our own audit rather than the literature, the grade says so.

- **Conditional losslessness (grade A).** `unstranded − sense − antisense` is lossless **only relative to the unstranded pass's *assigned* set**. Same-strand-overlap reads are dropped identically by all three modes and are **invisible** to the residual. A zero/positive residual does **not** certify no loss. A *positive* residual would be a bug (only a vanishing A/M/M edge does this). → `FLAG-RESIDUAL-NOT-LOSS`.

- **Never use `(sense + antisense)` as a quantity (grade A).** It double-counts the antiparallel pile (+1 to each channel, witnessed 100% cross-family) **and** the A/A/A pile. The sum-over-unstranded ratio is an artifact of double-counting, not "real bidirectional signal." → `FLAG-SUM-CHANNELS`.

- **Never use the unstranded channel as denominator for the LTR/ERV/Satellite/DNA subfamilies (grade A).** The unstranded pass drops their antiparallel reads to Ambiguity, giving a pathologically low (often negative-residual) denominator for these antiparallel-rich families. **Denominate on the SENSE channel** for them. → `FLAG-S0-DENOM-ERV`.

- **`excess/ambiguous` is directional, not an estimator (grade C → A on the negative).** Decomposition requires the per-read witness. The scalar identity is contaminated by the A/A/A leak + asymmetric multi-feature events + the sub-25bp minOverlap edge. → `FLAG-METER-NOT-ESTIMATOR`.

- **Kernel-consistency for normalization (grade A).** The gene matrix is **unique-only** (NO `-M`); every TE pass is **`-M`** (MultiMapping bucket = 0). Different kernels. Do **not** apply gene size factors to the TE matrix unqualified, and make no gene-vs-TE magnitude or "fraction of library" comparison — not like-for-like. → `FLAG-KERNEL-MISMATCH`. See [`METHODOLOGY.md` §Normalization and §Combined](./METHODOLOGY.md).

- **Silent loss is a *conditional* sensitivity tax (grade B, with a YELLOW open item).** Sample-stable silent loss cancels in cross-sample DE. **But** its share co-varies with per-sample net strand-offset (and weakly with depth) → **verify against the design matrix before declaring condition-independence.** The design matrix is not currently loaded; the check is pick-up-ready in the dataset-side `featurecounts_results/DE_precheck/` (README + per-sample strand QC + the check script + a design-matrix template). → `FLAG-SILENTLOSS-DESIGN`. This is the way to clear the YELLOW (see §7, §8).

- **`-O` target-revealer for the silent pile (grade A).** The default kernel writes the silent pile as `Unassigned_Ambiguity -1 NA` (no GeneID). Re-run the unstranded pass **kernel-identical except add `-O`** to reveal each fragment's candidate GeneID list, then attribute the silent set by subfamily (contains-young = upper bound; exclusively-young = floor). **`-O` is a target-revealer ONLY — it does NOT redefine the silent set** (fixed by the default-kernel statuses), and is a **diagnostic throwaway only, never the count matrix** (with `-O` the production matrix would double-assign). Witnessed: every silent fragment carries ≥2 candidate targets, and the young gate-set share is small and concordant with geometry. → procedure §4a. Closes the §8(ii) limitation at subfamily resolution for the witnessed sample (genic-context stratification stays Rung-2).

- **Random-One ⟹ the SENSE matrix is fragment-weighted (grade A).** Under STAR `--outSAMmultNmax 1` each fragment emits exactly **one** locus while `NH` records true multiplicity; confirmed by direct inspection across samples (no multimapping fragment emits more than one locus). Under `-M -R CORE` every Assigned row has `n_targets≡1`, so `matrix == distinct-Assigned-fragment count` for **every** subfamily — structural, generalizes across the dataset. **Counterfactual:** an alignment-weighted run would have inflated young high-copy families by their mean alignment multiplicity — averted, not realized. The A/A/A bucket is the benign antiparallel-feature strand-resolution this rests on (see §4, §9), inflating only `sense + antisense`, never the SENSE matrix. → `FLAG-ALIGNMENT-WEIGHTED`. See §9.

- **The ledger closes to the read (grade A).** The full (unstranded, sense, antisense) status-triple contingency sums **exactly** to N_fragments, and the excess and ambiguous-pile identities close to **0 reads** on the witnessed sample. A nonzero residual (>0.5%) ⟹ a hidden bucket — investigate; the audit's value is that the small cells (A/A/M, A/M/A, M/M/M, A/M/M, S/S/S) are each **named**, not lumped. The earlier first-order `excess = 2·AP + M` disagreement was **not** a hidden bucket — it was the A/A/A leak (+1) and violations (−1), now folded in. **Denominator discipline (sub-clause):** print numerator and denominator at the same scale and label it; per-sample and library-summed denominators must never mix; A/A/A is unstranded-ASSIGNED so `A/A/A ÷ ambiguous` is not a real fraction (see §3 caution). → `FLAG-DENOM-SCALE`. See §10.

---

## 7. GREEN/RED gate thresholds

**Young-autonomous candidates (L1Md_T/Gf/A, IAPEz) — conserved-fraction gate.**
- **GREEN ≥ 90% conserved** (non-conserved < 5%); **hard-RED at non-conserved ≥ 10%**.
- On the witnessed dataset the conserved fraction sits comfortably above the GREEN floor → **🟢 GREEN**; IAPEz-int is essentially immune.
- **The silent-loss leg is now WITNESSED at read level (grade A).** The `-O` target-revealer (§4a) attributes the silent pile by subfamily: the young gate-set share of silent loss is small and **concordant** with the read-independent annotation geometry (see §9 concordance). Folding even the upper-bound contains-young silent loss into the young denominator keeps conserved at the ≥90% floor; the defensible exclusively-young floor leaves it **comfortably GREEN**. The gate's second leg no longer rests on geometry alone — it is a per-read measurement. **(grade A — `-O` attribution, concordant with geometry.)**
- **Remaining open (Rung-2 GAP):** the *genic-context* stratification (intron/intergenic) of the young-silent reads still needs the richer SAF (see §8(ii)); the subfamily-level `-O` attribution is done.

**Mechanism gate.** Per-fragment `excess ∉ {0,+1,+2}` at `≥ 0.1%` → **RED** (instrument bug). Observed well under threshold → **pass** (see §4).

---

## 8. minOverlap sensitivity + the two open GAPs

**minOverlap sensitivity (grade A).** Raising `--minOverlap` to 25 drops roughly a quarter of the ambiguity pile: that fraction is a **sub-25bp annotation-edge artifact**, not biology. The regime split is robust to it. Report the **raw split as primary**; the minOverlap-25 split as the **edge-filtered lower bound** on the ambiguity denominator. This edge contaminates the directional meter (see §3).

**GAP (i) — silent-loss share vs the experimental design matrix.** Resolves the `FLAG-SILENTLOSS-DESIGN` YELLOW. The check is pick-up-ready in the dataset-side `featurecounts_results/DE_precheck/` (README + per-sample strand QC + the check script + a design-matrix template). The design matrix is not currently loaded. **GAP until done.**

**GAP (ii) — genic-context stratification of young-silent loss (re-scoped).** The **subfamily-level** attribution is now **DONE**: the `-O` target-revealer (§4a) witnesses the young gate-set silent-loss share, concordant with geometry (see §9), closing the limitation at subfamily resolution for the witnessed sample. **What remains GAP** is the **genic-context** axis — splitting those young-silent reads into intron / intergenic / adjacent via a **genic-context-stratified SAF** (`-R SAM`/BAM or coordinate intersection on the richer annotation). This is a **Rung-2** question; defer to [`LIMITATIONS.md` §5](./LIMITATIONS.md) for the ladder. Logged as Open Item 2 in the dataset-side `featurecounts_results/DE_precheck/README.md`. **GAP until done.**

---

## 9. Kernel-tracer: is the matrix fragment- or alignment-weighted?

This is the highest-stakes weighting question for young high-copy families: does the delivered SENSE (sense-pass) integer matrix count **one per fragment** (headline-safe) or **one per alignment** (young, highly-multimapping L1Md/IAPEz inflated by copy number)? The answer is **fragment-weighted**, confirmed by direct inspection of the alignments, in three parts.

**(1) Random-One is directly checkable in the alignments, not just a config pin (grade A).** Under STAR `--outMultimapperOrder Random --outSAMmultNmax 1` each fragment emits **exactly one** locus (one R1 + one R2, both `HI:i:1`) while `NH` still records the *true* multiplicity. Verify it directly: restricted to NH>1 fragments, count distinct emitted HI / secondary(0x100) lines / distinct (pos,mate) loci per readname. Across the witnessed samples, **no NH>1 fragment emits more than one locus**. This **upgrades** the "integer Random-One" kernel pin ([`METHODOLOGY.md` §Strategy A](./METHODOLOGY.md)) from a config claim to a witnessed fact.

**(2) Random-One ⟹ the SENSE integer matrix is fragment-weighted (grade A — matrix equals fragment count for essentially every subfamily).** Under `-M -R CORE` every Assigned row has `n_targets≡1`, so `matrix == distinct-Assigned-fragment count` for **every** subfamily (the matrix-over-ΣNH ratio sits far below any alignment-weighting, even for the highest-copy young families). This is **structural** (one emitted locus ⇒ counted once regardless of NH), so it generalizes across the dataset; spot-checking more only re-confirms.

**(3) Counterfactual risk — the headline-safety number (grade A).** Were the matrix alignment-weighted, young families would be **mean-alignment-multiplicity-fold larger** — a large fold for the youngest, highest-copy families, against near-unity for the old/rest. **None of that is realized** — counts equal distinct-fragment counts exactly. Record young-set mean alignment multiplicity as "the factor a *mis*-configured `--outSAMmultNmax all` + `-M` run WOULD have inflated by" — the explicit counterfactual that makes the kernel pin load-bearing (see the counterfactual-risk pattern in §3). → `FLAG-ALIGNMENT-WEIGHTED`: if BAMs ever emit >1 locus/fragment, re-verify before trusting young fold-changes.

**The A/A/A bucket is the benign mechanism this rests on.** A/A/A (unstranded-Assigned, and both stranded passes Assigned) is benign strand-resolved antiparallel-feature selection on single-locus reads, almost entirely uniquely-mapped (full account in §4) — consistent with Random-One, not evidence against it. Each A/A/A fragment is +1 to exactly one sense GeneID; the divergence inflates only `sense + antisense` (`FLAG-SUM-CHANNELS`), never the SENSE matrix.

### 9a. Geometry-vs-witness concordance (validating the read-level gate)

**The method (grade A — per-read `-O` attribution vs a read-independent bedtools self-intersect).** Validate a read-level gate (here the young-silent attribution of §4a/§7) against **read-independent annotation geometry** — the bedtools self-intersect of the SAF (antiparallel:parallel pair/bp ratios, per-class/age-bin localization). **Concordance corroborates** the gate (the witness is not an artifact of read placement); **divergence flags** and must be resolved — as either a real read-density effect (expression × geometry) or a bin-definition artifact — *before* the number is trusted.

**Worked example.** The young-silent witness (the gate-set silent-loss share) is **concordant** with the geometry prediction across the young classes, close to exact for the youngest ERVK/IAP and gate-L1 bins. The one apparent divergence — a broad-L1Md bin running high against geometry — resolves **entirely** as the broad annotation bin pulling in older L1Md subfamilies that are not gate-young; restricting to the gate set brings the share back below geometry. The residual direction (reads slightly over-represented in older expressed copies vs flat geometry) is the same **expression × geometry ≠ geometry-alone** effect seen in the antiparallel-to-parallel ratio. **Rule: concordance is corroboration; divergence is either a real read-density effect or a bin-definition artifact — diagnose which before trusting it.** This is the read-level evidence behind the §7 GREEN verdict.

---

## 10. Closure-table audit (the ledger closes to the read)

After the per-read witness (§4), **prove the ledger closes to the read.** Build the full (unstranded, sense, antisense) status-triple contingency and check **(grade A — zero-read residual on the audited sample):**

1. it **sums exactly to N_fragments** (every fragment in exactly one cell, no remainder);
2. `excess = 2·AP + M + AAA − violations` closes to **0 reads**;
3. `ambiguous = AP + M + P(parallel) + bilateral` closes to **0 reads** — confirming A/A/A is **NOT** part of the ambiguous pile (it is unstranded-ASSIGNED);
4. `sense_Amb` / `antisense_Amb` reconstruct from the same cells (residual 0 each).

**A nonzero residual (>0.5% threshold) ⟹ a hidden bucket — investigate it.** The audit's value is that the small cells are each **NAMED**, not lumped: `A/A/M` and `A/M/A` are unstranded-Assigned with one stranded arm Ambiguous (they enter sense_Amb/antisense_Amb but **not** the ambiguous pile or excess); `M/M/M` bilateral (in the ambiguous pile, 0 to excess); `A/M/M` the truth-table violation (the −1 debit); `S/S/S` the singleton floor (off-ledger). Every fragment lands in a named cell — there is no extra bucket. The earlier first-order `excess = 2·AP + M` disagreement was simply the A/A/A leak (+1) and violations (−1), both folded into identity (2) above, which then closes exactly.

**Denominator discipline (sub-clause; see §3 caution).** Always print numerator and denominator at the **same scale** and **label the scale**. A per-sample numerator and a library-summed denominator are two different scales and must never be combined; nor may the two share-of-excess and share-of-ambiguous ratios be mixed. A/A/A is **unstranded-ASSIGNED**, so `A/A/A ÷ ambiguous` is **not** a real fraction of the ambiguous pile — do not report it. → `FLAG-DENOM-SCALE`.

> **Note (read-level closure of record).** The high-depth read-level closure is the closure of record; a second sample's marginals are consistent with the same regime structure, but its per-read assignment files were not retained, so its full status-triple cannot be rebuilt from witnessed data — re-run the per-read witness on its BAM to reproduce.

---

## Resources

- [`METHODOLOGY.md`](./METHODOLOGY.md) — multimapper kernel rationale (Strategy A/B, Teissandier), normalization, combined-matrix validity.
- [`LIMITATIONS.md`](./LIMITATIONS.md) — the 3-rung ladder (§5), the co-oriented-readthrough residue (§2), the necessary-not-sufficient antisense channel and global read-through sentinel (§4), and the same-strand silent loss limitation (§2a).
- [`BIOLOGY.md`](./BIOLOGY.md) — the biology behind the strand channels: sense/antisense, class-specific bidirectional transcription, read-through, derepression.
- `te-gene-featurecounts` skill — strand-split QC regression: `tests/strand_qc/run_regression.sh` (synthetic truth table in `te-fc:2.0.2` + the `-R CORE` regime-classifier fixture); reference `references/strand-split-qc.md`.
- `te-gene-featurecounts` skill `qc/` — the runnable, parameterized "QC a new TE dataset end-to-end" suite (`qc/run_qc.sh`: BAM dir + SAF + strand + outdir → witness tables + GREEN/RED verdicts). Operationalizes the `-R CORE` regime witness (§4), the kernel-tracer / fragment-weighting check (§9), the closure-table audit (§10), and the `-O` silent-loss attribution (§4a). See `qc/README.md`.
- `featurecounts_results/DE_precheck/` (dataset-side, not in this repo) — the pick-up-ready DE-confound check for GAP (i) / `FLAG-SILENTLOSS-DESIGN`.
- `featurecounts_results/MATRICES.md` (dataset-side, not in this repo) — the five-matrix manifest with per-matrix kernel provenance and do-not flags.
