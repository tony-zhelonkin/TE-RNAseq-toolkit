# Strand-split TE counting — QC, invariants, and foot-guns

This document is the **QC gate** and the **warning taxonomy** for the three-channel (s0 / sense / anti) TE counting this toolkit performs. It owns the strand-split invariant, the directional meter, the `-R CORE` regime-witness procedure, the numbered axioms, the warning-flag taxonomy, and the GREEN/RED gate thresholds. It is downstream of the subfamily-resolution fact stated in [`LIMITATIONS.md`](./LIMITATIONS.md): the grouped `Subfamily:Family:Class` SAF carries no genic context and no locus identity, and every boundary below sits on top of that.

Grades below use the **A/B/C/D/GAP** scale of [`METHODOLOGY.md`](./METHODOLOGY.md); its citations are not restated here. This document **defers** to [`METHODOLOGY.md`](./METHODOLOGY.md) for the multimapper-kernel rationale (Strategy A/B, Teissandier) and to [`LIMITATIONS.md`](./LIMITATIONS.md) for the 3-rung ladder and the co-oriented-readthrough residue — links, not duplication.

> **One thing to fix up front.** Same-strand-overlap silent loss is a **NEW** limitation that lives **alongside** — not replacing — the co-oriented-readthrough residue of [`LIMITATIONS.md` §2](./LIMITATIONS.md). The two are distinct (silent loss is an annotation-geometry tax; read-through is a transcription-context confound) and neither subsumes the other. See [`LIMITATIONS.md` §2a](./LIMITATIONS.md).

**Convention (reverse-stranded dUTP library, 45 samples, mm39).** `sense = -s2 Assigned`, `anti = -s1 Assigned`, `s0 = -s0` unstranded. featureCounts v2.0.2, all TE passes `-M -F SAF -p --countReadPairs -B -C`, no `-O`, no `--fraction`, default `--minOverlap 1`, integer.

---

## 1. The strand-split invariant (conditional losslessness)

Define, per subfamily and library-wide:

```
s0 ≈ sense + anti + residual
residual_frac = residual / s0
```

with `sense = -s2 Assigned`, `anti = -s1 Assigned`, `s0 = -s0 Assigned`.

**Axiom A1 — conditional losslessness. [evidence: A — per-read `-R CORE` witness, reproduced on a 2nd sample.]** The residual `s0 − sense − anti` is lossless **only relative to s0's *assigned* set**. A read overlapping two **same-strand** features is dropped as Ambiguity by `-s 0` **and** by both stranded passes (the strand-matched pass still sees ≥2 candidates → Ambiguity; the opposite pass sees none → NoFeatures). It contributes 0 to s0, sense, and anti alike, so `residual = 0`. This **silent loss** is by construction **invisible** to the residual. A zero or one-sided (never positive) residual therefore certifies **nothing** about loss. A *positive* residual would be a bug — a stranded pass dropping a read that s0 assigned — and only the 0.0004% A/M/M `-B` proper-pair edge does this.

**The refuted framing, stated plainly.** The original QC report (`strand_split_invariant_QC.md`, a dataset-side artifact, not in this repo) concluded the split is **"lossless"** and that the −13% net offset is **"entirely overlap RECOVERY."** **That framing is REFUTED.** The residual sees only the antiparallel (+2) and asymmetric (+1) slices it *can* see; it is structurally blind to same-strand silent loss. The visible −13% recovery is real but sits on top of an **unseen silent-loss pile roughly twice its size** (§2). The residual being one-sided is *consistent with*, not evidence against, large silent loss. The authoritative reconciliation is the dataset-side `_strand_invariant/SYNTHESIS_verdict.md` (not in this repo).

---

## 2. The witnessed ledger (one table)

Per-read `-R CORE` witness, sample 0019 (high-depth); reproduced on 0006 (low-depth). The s0-ambiguous pile = 657,884 fragments in 0019, ≈17.5M library-wide.

Each fragment is labelled by how many features the **sense** and **antisense** passes each assign it — written `(sense, antisense)`. A strand-ambiguous fragment then meets one of three fates: **AP — antiparallel** `(1,1)`, kept by *both* passes under *different* subfamilies (double-presence — what summing the channels would double-count); **M — asymmetric** `(2,1)/(1,2)`, where splitting by strand breaks the tie and the fragment is genuinely *reclaimed*; **P — parallel / same-strand** `(n,0)`, separable by no strand flag and so *silently lost* in all three passes — the dominant fate. The two numeric columns below use **different denominators** (share of the ambiguous pile vs share of total assigned signal) — do not mix them.

| Regime | % of s0_Amb | Library-scale (of s0_Assigned) | Grade | Meaning |
|---|---:|---:|:---:|---|
| **AP — antiparallel (1,1)** | **25.7%** (0006: 24.7%) | ≈ **4%** | **A** | double-presence; +1 to EACH stranded pass, **different** GeneIDs (100% cross-family, 0% same-GeneID) |
| **M — asymmetric (2,1)/(1,2)** | **6.1%** (0006: 5.9%) | ≈ **1%** | **A** | single reclaim = **genuine recovery**; +1 to excess |
| **P — silent loss** (parallel n,0 + bilateral ≥2,≥2) | **68.2%** (0006: 69.3%) | ≈ **9–10%** (per-sample 8.8–10.9%) | **A** | residual-invisible loss; bilateral only **0.38%** → single-strand-drop assumption **HOLDS** |

**Locked one-liner:** s0-ambiguous pile = **AP 25.7% : M 6.1% : P(silent) 68.2%**; library-scale silent loss ≈ **9–10% of assigned TE signal**, double-presence ≈ 4%, genuine reclaim ≈ 1%. Silent loss is the **dominant** fate and is **larger** than any scalar solve implied. Bilateral at 0.38% means the single-strand-drop assumption holds **[evidence: A]**.

---

## 3. The `excess/s0_Amb` directional meter

**Definition.** Over the s0-ambiguous pile, `excess = (sense_Assigned + anti_Assigned) − s0_Assigned`; the meter is `excess / s0_Amb ≈ 0.78` (mean 0.779, SD 0.015, range [0.739, 0.809], n = 45).

**Graded DIRECTIONAL-ONLY (A4). [evidence: C — scalar/directional; A on the negative claim that it is not an estimator.]** The meter is real and sample-stable, and it tracks the antiparallel-vs-silent balance **directionally**, but it does **not** yield the split. It is contaminated by (i) the **A/A/A leak** (s0-Assigned fragments that BOTH stranded passes also assign to a *different* feature — 25.4% of excess), (ii) asymmetric multi-feature events, and (iii) the sub-25bp `--minOverlap` annotation edge (24% of the pile, §8). **Use it as a meter, never as an estimator.**

**The methodological lesson (stated once).** An elegant scalar 3×3 solve — `excess = 2·AP + M`, `s0_Amb = AP + M + P`, `ambsum = M + P` ⇒ AP 24% / M 30% / P 46% — **reconciled to the integer** on the s0_Amb identity (+0.000%) and to ~2.5% on ambsum, yet was **wrong**: witnessed P is 68%, not 46%, and witnessed M is ~6%, not ~30%. The scalar solve cannot see the A/A/A leak, which adds +1 each to excess and is folded invisibly into M; that leak is the sole reason `excess = 2·AP + M` fails by **−25.31%** (`2·AP 338,064 + M 40,362 + A/A/A 128,422 − violations 192 = 506,656` = measured excess exactly). **Grade the `excess = 2·AP + M` identity REFUTED as an estimator; the mechanism `excess ∈ {0,+1,+2}` is upheld (§4).** **Lesson: for overlap accounting, witness > inference** — a residual that balances does not certify a mechanism; it can balance because two unseen errors cancel.

**The counterfactual-risk pattern (the second cross-cutting lesson). [evidence: A — method.]** When a kernel choice could silently inflate or deflate a headline, **quantify the averted-error magnitude explicitly, then witness that the realized number is the safe one.** A pin you cannot put a downside number on is not a QC gate. The load-bearing instance is alignment- vs fragment-weighting of young high-copy families (§9): the matrix is fragment-weighted, but a mis-configured `--outSAMmultNmax all` + `-M` run *would have* inflated young families by their mean alignment multiplicity (young-set `meanNH = 5.08`; IAPEz **16.6×**, MMERVK10C 9.4×). Recording that counterfactual fold is what makes the kernel pin load-bearing — it states what was at stake and that the witness confirms it was averted. This pairs with "witness > inference" as the second cross-cutting doctrine: do not just prove the number is right; prove the wrong number you avoided was large.

**Denominator discipline (caution). [evidence: A.]** Always print numerator and denominator at the **same scale** and **label the scale** on every ratio (per-sample vs library-summed). The per-sample-0019 "25.4%" (A/A/A ÷ per-sample excess) and "25.7%" (AP ÷ per-sample s0_Amb) are two different per-sample numerators over two different per-sample denominators; the library-summed **17.5M** s0_Amb (all 45 samples) must **never** denominate a single-sample (0019) numerator. A/A/A is **s0-ASSIGNED**, so `A/A/A ÷ s0_Amb` is **not** a real fraction of the ambiguous pile — do not report it. → `FLAG-DENOM-SCALE` (A9, §10).

---

## 4. The `-R CORE` regime-witness procedure (reproducible)

Run featureCounts **3× on the SAME BAM** (`-s 0` / `-s 2` / `-s 1`), **identical kernel** (`-M -F SAF -p --countReadPairs -B -C`, no `-O`, no `--fraction`), each with `-R CORE` to emit the per-fragment assignment table.

> **WARNING (load-bearing foot-gun).** Join the three per-fragment tables by **hash-join on fragment ID**, NOT positional `paste`. Read-ID order differs across threaded featureCounts runs; a positional paste silently misaligns fragments and fabricates regimes. This is the single easiest way to corrupt the witness.

Classify each fragment's `(s0, s2, s1)` Assigned-count status triple into a regime:

| Regime | Triple | Excess | Meaning |
|---|---|:--:|---|
| conserved | (1,0) / (0,1) | 0 | one strand, one feature |
| antiparallel (AP) | (1,1) | +2 | double-presence, different GeneIDs |
| asymmetric (M) | (2,1) / (1,2) | +1 | single reclaim = **genuine recovery** |
| parallel-silent (P) | (n,0) | 0 | **SILENT LOSS** (residual-invisible) |
| bilateral | (≥2,≥2) | 0 | silent loss (negligible, 0.38%) |
| **A/A/A leak** | s0 Assigned, both stranded passes Assign a *different* feature | +1 each | the bucket the scalar solve folds invisibly into M — **witnessed benign**, see below |

Reconstruct `excess` from the regimes and apply the gate.

**What the A/A/A bucket actually is (witnessed, not arithmetic-only). [evidence: A — whole-BAM `-R CORE` hash-join, sample 0019.]** Earlier wording called A/A/A only "the arithmetic cross-channel leak." The per-read witness resolves it: A/A/A is **benign strand-resolved antiparallel-feature selection on single-locus reads**. Of the 114,648 fragments whose sense and antisense GeneIDs *diverge*, **99.1% are NH==1 (uniquely mapped)** — per-alignment double counting is categorically impossible for a uniquely-mapped fragment — and the s0 (unstranded) GeneID equals **exactly one** of {sense, anti} for **100%** of them (≈47% == sense, ≈53% == anti, never both, never neither). The single fragment overlaps **two annotated TE features on opposite genomic strands** (nested/abutting antiparallel TEs); featureCounts derives fragment strand from R1, so the sense pass legitimately picks the +strand feature and the antisense pass the −strand feature — same locus, one alignment, two strand-consistent features. This **REFUTES** the reviewer premise that "s0 assigns ⟹ a single GeneID ⟹ sense/anti cannot diverge ⟹ the BAMs are not Random-One": unstranded *can* assign one feature while the two stranded passes legitimately pick the two opposite-strand features. The divergence inflates only `sense + anti` — a sum you must never form (`FLAG-SUM-CHANNELS`, A2) — and **never the SENSE matrix itself** (each A/A/A fragment is +1 to exactly one sense GeneID). See §9 for the kernel-tracer and the fragment-weighting theorem this rests on.

**Mechanism gate: per-fragment `excess ∈ {0,+1,+2}`.** Violation threshold `< 0.1%`; observed **0.00042%** (0019) / 0.00026% (0006) **[evidence: A]**. Agrees with the synthetic truth table (cases 2/4/3,6 → +2/+1/0). The single violation cell (A/M/M, excess −1, 192 fragments) is a `-B`/proper-pair overlap edge and is negligible.

This procedure and its truth table are regression-tested in the `te-gene-featurecounts` skill, `tests/strand_qc/` (see §-Resources).

### 4a. Attributing the silent pile with `-O` (the target-revealer procedure)

**The problem.** Under the default kernel (`-R CORE`, no `-O`) featureCounts writes the silent pile as `Unassigned_Ambiguity  -1  NA` — the candidate target list is **not emitted**, so silent-loss fragments carry **no GeneID** (the `UNATTRIB_silent` limitation flagged in §7). To attribute the pile by subfamily you must recover the features each silent fragment was ambiguous **among**.

**The fix (kernel-identical re-run + one flag). [evidence: A — per-read `-O` CORE on the same BAM; 448,490/448,490 silent fragments matched a multi-target record (sample 0019).]** Re-run the s0 pass **kernel-identical except add `-O`** (`-M -O -F SAF -s 0 -p --countReadPairs -B -C -R CORE`). The alignment total reproduces the audited run exactly; every fragment that was `Unassigned_Ambiguity` in the default run becomes `Assigned n_targets≥2` with the comma-separated GeneID list = the candidate features it overlapped. Then look up each silent fragment's target list in the `-O` CORE by fragment ID.

**`-O` is a target-revealer ONLY — it does NOT redefine the silent set.** The silent set stays fixed by the **default-kernel** s0/s2/s1 statuses (the P bucket of §2); `-O` only attaches a candidate-target list to fragments already classified silent. Attribution policy: **contains-young** (≥1 young target) = the upper-bound witness; **exclusively-young** (all targets young) = the floor.

> **Foot-gun — never collide the two `-O` regimes.** The production count matrices correctly use **no `-O`** (with `-O` a fragment would double-assign across all its targets and inflate counts). A7/§4a is the **diagnostic** use of `-O` on a **throwaway CORE re-run**, never on the count matrix. Keep the two uses strictly apart.

This converts the young-silent attribution from read-independent geometry inference to a per-read witness (§9 geometry-vs-witness concordance; the gate-set share is **1.27–1.37%**, concordant). It closes the §8(ii) CORE-NA limitation at subfamily resolution for the witnessed sample; the genic-context stratification (intron/intergenic) remains a Rung-2 GAP.

---

## 5. The warning-flag taxonomy (THE canonical table)

This is the **single source of truth** for the named flags. Sidecars and skill docs cite them **by name only** — never re-defining them. Severities: **🔴 hard** (violating it produces a wrong number or fabricated biology) · **🟡 watch** (conditional — must be checked). Each flag is a measurement **assertion**, never a biological conclusion.

| Flag name | Sev | Axiom | One-line assertion |
|---|:--:|:--:|---|
| `FLAG-SUM-CHANNELS` | 🔴 | A2 | Never use `(sense + anti)` as a quantity — it double-counts the antiparallel pile (+1 to each channel, 100% cross-family) **and** the A/A/A pile; `(s+a)/s0 = 1.13` is a double-counting artifact, not bidirectional signal. |
| `FLAG-S0-DENOM-ERV` | 🔴 | A3 | Never s0-denominate the LTR/ERV/Satellite/DNA subfamilies — s0 drops their antiparallel reads to Ambiguity (pathological denominator); denominate on the **SENSE** channel. |
| `FLAG-KERNEL-MISMATCH` | 🔴 | A5 | Gene matrix is unique-only (no `-M`, 111,888,517 MM pairs dropped); every TE pass is `-M`. Do not apply gene size factors to the TE matrix unqualified; no gene-vs-TE magnitude / "fraction of library" — not like-for-like. |
| `FLAG-RESIDUAL-NOT-LOSS` | 🟡 | A1 | A zero/positive `s0 − sense − anti` residual does NOT certify no loss; same-strand silent loss (≈9–10% of assigned TE signal) is residual-invisible by construction. |
| `FLAG-METER-NOT-ESTIMATOR` | 🟡 | A4 | `excess/s0_Amb ≈ 0.78` is a directional meter, never an estimator of the AP/M/P split (A/A/A + asymmetric + minOverlap-edge contamination). |
| `FLAG-SILENTLOSS-DESIGN` | 🟡 | A6 | Silent loss is sample-stable (cancels under usual conditions) BUT P_frac co-varies with net strand-offset (r = −0.50); verify against the design matrix before declaring condition-independence (`DE_precheck/`). |
| `FLAG-ALIGNMENT-WEIGHTED` | 🔴 | A8 | The SENSE integer matrix is fragment-weighted under Random-One (`--outSAMmultNmax 1`, matrix==nfrag); if BAMs ever emit >1 locus/fragment (or `--outSAMmultNmax` > 1 with `-M`), young high-NH families inflate up to ~meanNH-fold (IAPEz 16.6×) — re-verify with the kernel-tracer (§9) before trusting young fold-changes. |
| `FLAG-DENOM-SCALE` | 🟡 | A9 | Per-sample vs library-summed denominators must never be mixed; label the scale on every ratio. A/A/A is s0-ASSIGNED, so `A/A/A ÷ s0_Amb` is not a real fraction of the ambiguous pile — do not report it. |

Names are stable identifiers — do not rename downstream.

**The ratio caveat (folds into the table above).** `(s+a)/s0 = 1.1287` is **NOT** bidirectional signal — it double-counts the antiparallel pile **and** the A/A/A pile (`FLAG-SUM-CHANNELS` / A2). For the LTR/ERV/Satellite/DNA subfamilies, denominate on the **SENSE channel**, not s0, because s0 drops their antiparallel reads to Ambiguity and gives a pathologically low denominator (`FLAG-S0-DENOM-ERV` / A3). The `sense/unstr` and `anti/unstr` columns in `qc_ratios.tsv` are **QC meters, not normalization denominators**.

---

## 6. The axioms A1–A6 (numbered, graded)

- **A1 — Conditional losslessness (grade A).** `s0 − sense − anti` is lossless **only relative to s0's *assigned* set**. Same-strand-overlap reads are dropped identically by all three modes and are **invisible** to the residual. A zero/positive residual does **NOT** certify no loss. A *positive* residual would be a bug (only the 0.0004% A/M/M edge does this). → `FLAG-RESIDUAL-NOT-LOSS`.

- **A2 — Never use `(sense + anti)` as a quantity (grade A).** It double-counts the antiparallel pile (+1 to each channel, witnessed 100% cross-family) **and** the A/A/A pile. `(s+a)/s0 = 1.13` is an artifact of double-counting, not "real bidirectional signal." → `FLAG-SUM-CHANNELS`.

- **A3 — Never s0-denominate the LTR/ERV/Satellite/DNA subfamilies (grade A).** s0 drops their antiparallel reads to Ambiguity, giving a pathologically low denominator (CYRA11 `residual_frac −4.185`; RMER6BA −1.245; RLTR/RMER ERVK −0.44 to −0.59). **Denominate on the SENSE channel** for these families. → `FLAG-S0-DENOM-ERV`.

- **A4 — `excess/s0_Amb` is directional, not an estimator (grade C → A on the negative).** Decomposition requires the per-read `-R CORE` witness. The scalar identity is contaminated by the A/A/A leak + asymmetric multi-feature events + the sub-25bp minOverlap edge. → `FLAG-METER-NOT-ESTIMATOR`.

- **A5 — Kernel-consistency for normalization (grade A).** Gene matrix is **unique-only** (NO `-M`, 111,888,517 multimapper pairs dropped); every TE pass is **`-M`** (MultiMapping bucket = 0). Different kernels. Do **not** apply gene size factors to the TE matrix unqualified; no gene-vs-TE magnitude or "fraction of library" comparison — not like-for-like. → `FLAG-KERNEL-MISMATCH`. See [`METHODOLOGY.md` §Normalization and §Combined](./METHODOLOGY.md).

- **A6 — Silent loss is a *conditional* sensitivity tax (grade B, with a YELLOW open item).** Sample-stable silent loss (P_frac-of-signal SD 0.0028, below net-offset SD 0.0081) cancels in cross-sample DE. **BUT** P_frac co-varies with per-sample net strand-offset (**r = −0.50**) and weakly with depth (r = −0.11) → **verify against the design matrix before declaring condition-independence.** The design matrix is not currently loaded; the check is pick-up-ready in the dataset-side `featurecounts_results/DE_precheck/` (README + `per_sample_strand_qc.tsv` + `check_silentloss_vs_design.py` + `design_matrix.TEMPLATE.tsv`). → `FLAG-SILENTLOSS-DESIGN`. This is the way to clear the YELLOW (§7, §8).

- **A7 — `-O` target-revealer for the silent pile (grade A).** The default kernel writes the silent pile as `Unassigned_Ambiguity -1 NA` (no GeneID). Re-run the s0 pass **kernel-identical except add `-O`** to reveal each fragment's candidate GeneID list, then attribute the silent set by subfamily (contains-young = upper bound; exclusively-young = floor). **`-O` is a target-revealer ONLY — it does NOT redefine the silent set** (fixed by the default-kernel statuses), and is a **diagnostic throwaway CORE only, never the count matrix** (with `-O` the production matrix would double-assign). Witnessed: 448,490/448,490 silent fragments carry ≥2 candidate targets; young gate-set share 1.27–1.37%, concordant with geometry. → procedure §4a. Closes the §8(ii) CORE-NA limitation at subfamily resolution for the witnessed sample (genic-context stratification stays Rung-2).

- **A8 — Random-One ⟹ the SENSE matrix is fragment-weighted (grade A).** Under STAR `--outSAMmultNmax 1` each fragment emits exactly **one** locus while `NH` records true multiplicity; BAM-witnessed (0/3,775,221 NH>1 frags emit >1 locus on 0019; 0 on 0001/0033). Under `-M -R CORE` every Assigned row has `n_targets≡1`, so `matrix == distinct-Assigned-fragment count` for **every** subfamily (matrix/ΣNH = 0.06 for IAPEz) — structural, generalizes to all 45 samples. **Counterfactual:** an alignment-weighted run would have inflated young families by their meanNH (young set 5.08; IAPEz 16.6×, MMERVK10C 9.4×) — averted, not realized. The A/A/A bucket is the benign antiparallel-feature strand-resolution this rests on (§4, §9), inflating only `sense + anti`, never the SENSE matrix. → `FLAG-ALIGNMENT-WEIGHTED`. See §9.

- **A9 — The ledger closes to the read (grade A).** The full `(s0,s2,s1)` 14-cell status-triple contingency sums **exactly** to N_fragments, and the identities `excess = 2·AP + M + AAA − violations` and `s0_Amb = AP + M + P + bilateral` close to **0 reads** (sample 0019). A nonzero residual (>0.5%) ⟹ a hidden fifth bucket — investigate; the audit's value is that the small cells (A/A/M, A/M/A, M/M/M, A/M/M, S/S/S) are each **named**, not lumped. The earlier `excess = 2·AP + M` disagreement was **not** a hidden bucket — it was the A/A/A leak (+1) and violations (−1), now folded in. **Denominator discipline (sub-clause):** print numerator and denominator at the same scale and label it; per-sample-0019 vs library-summed (17.5M) must never mix; A/A/A is s0-ASSIGNED so `A/A/A ÷ s0_Amb` is not a real fraction (§3 caution). → `FLAG-DENOM-SCALE`. See §10.

---

## 7. GREEN/RED gate thresholds

**Young-autonomous candidates (L1Md_T/Gf/A, IAPEz) — conserved-fraction gate.**
- **GREEN ≥ 90% conserved** (non-conserved < 5%); **hard-RED at non-conserved ≥ 10%**.
- Observed **96.65% conserved / 3.35% non-conserved** → **🟢 GREEN**. IAPEz-int 99.8% conserved (essentially immune).
- **The silent-loss leg is now WITNESSED at read level (grade A).** The `-O` target-revealer (§4a, A7) attributes the silent pile by subfamily: the young gate-set share is **1.27–1.37%** of silent loss (exclusively-young floor 0.285%), **concordant** with the read-independent ev-annotation geometry (~1.0–1.6%; IAPEz 0.62% vs 0.6%, gate-L1 0.93% vs ~1.0% — see §9 concordance). Folding even the upper-bound contains-young silent loss into the young denominator keeps conserved at **90.4% (≥90% floor)**; the defensible exclusively-young floor leaves it **comfortably GREEN (95.3% conserved)**. The gate's second leg no longer rests on geometry alone — it is a per-read witness. **[evidence: A — `-O` per-read witness, concordant with geometry.]**
- **Remaining open (Rung-2 GAP):** the *genic-context* stratification (intron/intergenic) of the young-silent reads still needs the richer SAF (§8(ii)); the subfamily-level `-O` attribution is done.

**Mechanism gate.** Per-fragment `excess ∉ {0,+1,+2}` at `≥ 0.1%` → **RED** (instrument bug). Observed 0.00042% → **pass** (§4).

---

## 8. minOverlap sensitivity + the two open GAPs

**minOverlap sensitivity (grade A).** At `--minOverlap 25`, s0_Amb falls **657,884 → 500,141**: **24.0% of the ambiguity pile is a sub-25bp annotation-edge artifact**, not biology. The regime split is robust to it. Report the **raw split as primary**; the mO25 split as the **edge-filtered lower bound** on the ambiguity denominator. This edge contaminates the directional meter (§3).

**GAP (i) — P_frac vs the experimental design matrix.** Resolves the A6 / `FLAG-SILENTLOSS-DESIGN` YELLOW. The check is pick-up-ready in the dataset-side `featurecounts_results/DE_precheck/` (README + `per_sample_strand_qc.tsv` + `check_silentloss_vs_design.py` + `design_matrix.TEMPLATE.tsv`). The design matrix is not currently loaded. **GAP until done.**

**GAP (ii) — genic-context stratification of young-silent loss (re-scoped).** The **subfamily-level** attribution is now **DONE**: the `-O` target-revealer (§4a, A7) witnesses the young gate-set silent-loss share at **1.27–1.37%**, concordant with geometry (§9), closing the CORE-NA limitation at subfamily resolution for the witnessed sample. **What remains GAP** is the **genic-context** axis — splitting those young-silent reads into intron / intergenic / adjacent via a **genic-context-stratified SAF** (`-R SAM`/BAM or coordinate intersection on the richer annotation). This is a **Rung-2** question; defer to [`LIMITATIONS.md` §5](./LIMITATIONS.md) for the ladder. Logged as Open Item 2 in the dataset-side `featurecounts_results/DE_precheck/README.md`. **GAP until done.**

---

## 9. Kernel-tracer: is the matrix fragment- or alignment-weighted?

This is the highest-stakes weighting question for young high-copy families: does the delivered SENSE (`-s 2`) integer matrix count **one per fragment** (headline-safe) or **one per alignment** (young, highly-multimapping L1Md/IAPEz inflated by copy number)? The answer is **fragment-weighted**, BAM-witnessed in three parts.

**(1) Random-One is BAM-witnessable, not just pinned. [evidence: A — whole-BAM, 3 samples.]** Under STAR `--outMultimapperOrder Random --outSAMmultNmax 1` each fragment emits **exactly one** locus (one R1 + one R2, both `HI:i:1`) while `NH` still records the *true* multiplicity. Verify it directly: restricted to NH>1 fragments, count distinct emitted HI / secondary(0x100) lines / distinct (pos,mate) loci per readname. Observed **0/3,775,221 fragments emit >1 locus on 0019; 0 on 0001 and 0033**. This **upgrades** the "integer Random-One" kernel pin ([`METHODOLOGY.md` §Strategy A](./METHODOLOGY.md)) from a config claim to a witnessed fact.

**(2) Theorem — Random-One ⟹ the SENSE integer matrix is fragment-weighted. [evidence: A — matrix==nfrag for all 1,140/1,146 subfamilies (0019/0033).]** Under `-M -R CORE` every Assigned row has `n_targets≡1`, so `matrix == distinct-Assigned-fragment count` for **every** subfamily (matrix/ΣNH = **0.06** for IAPEz — far below any alignment-weighting). This is **structural** (one emitted locus ⇒ counted once regardless of NH), so it generalizes to all 45 samples; spot-checking more only re-confirms.

**(3) Counterfactual risk (the headline-safety number). [evidence: A.]** Were the matrix alignment-weighted, young families would be **meanNH-fold larger**: IAPEz **16.6×**, MMERVK10C **9.4×**, young set ~5× (meanNH/assigned-frag = **5.08**; old/rest meanNH 1.18). **None of that is realized** — counts equal distinct-fragment counts exactly. Record young-set `meanNH` as "the factor a *mis*-configured `--outSAMmultNmax all` + `-M` run WOULD have inflated by" — the explicit counterfactual that makes the kernel pin load-bearing (§3 counterfactual-risk pattern). → `FLAG-ALIGNMENT-WEIGHTED` (A8): if BAMs ever emit >1 locus/fragment, re-verify before trusting young fold-changes.

**The A/A/A bucket is the benign mechanism this rests on.** A/A/A (s0-Assigned and both stranded passes Assigned) is **not** divergent double-counting and **not** evidence against Random-One — it is benign strand-resolved antiparallel-feature selection on single-locus reads (99.1% NH==1), refuting the not-Random-One premise (full witness in §4). Each A/A/A fragment is +1 to exactly one sense GeneID; the divergence inflates only `sense + anti` (`FLAG-SUM-CHANNELS`, A2), never the SENSE matrix.

### 9a. Geometry-vs-witness concordance (validating the read-level gate)

**The method. [evidence: A — per-read `-O` witness vs read-independent bedtools self-intersect.]** Validate a read-level gate (here the young-silent attribution of §4a/§7) against **read-independent annotation geometry** — the bedtools self-intersect of the SAF (antiparallel:parallel pair/bp ratios, per-class/age-bin localization). **Concordance corroborates** the gate (the witness is not an artifact of read placement); **divergence flags** and must be resolved — as either a real read-density effect (expression × geometry) or a bin-definition artifact — *before* the number is trusted.

**Worked example.** The young-silent witness (gate set **1.27–1.37%**) is **concordant** with the geometry prediction (~1.0–1.6%; YOUNG_ERVK_IAP 0.62% witnessed vs 0.6% geometry ≈ exact; gate-L1 Tf/Gf/A 0.93% vs ~1.0%). The lone apparent divergence — broad-L1Md **2.52% vs geometry 1.0%** (~2.5×) — resolves **entirely** as ev-annotation's broad `YOUNG_L1Md` bin pulling in the OLDER L1MdF/V/Mus/Fanc (top drivers L1MdF_IV, L1MdV_III, L1MdMus_I, L1MdFanc_II — all older, not gate-young); restricting to the gate set gives 0.93%, *below* geometry. The residual direction (reads slightly over-represented in older expressed copies vs flat geometry) is the same **expression × geometry ≠ geometry-alone** effect seen in the antiparallel:parallel ratio (0.27 obs vs 0.39 geom). **Rule: concordance is corroboration; divergence is either a real read-density effect or a bin-definition artifact — diagnose which before trusting it.** This is the read-level evidence behind the §7 GREEN verdict.

---

## 10. Closure-table audit (the ledger closes to the read)

After the `-R CORE` witness (§4), **prove the ledger closes to the read.** Build the full `(s0, s2, s1)` status-triple contingency (14 cells) and check: **[evidence: A — full 14-cell contingency, 0-read residual on 0019.]**

1. it **sums exactly to N_fragments** (every fragment in exactly one cell — 45,913,779 on 0019, no remainder);
2. `excess = 2·AP + M + AAA − violations` closes to **0 reads** (`2·169,032 + 40,362 + 128,422 − 192 = 506,656` = measured excess exactly);
3. `s0_Amb = AP + M + P(parallel) + bilateral` closes to **0 reads** (`169,032 + 40,362 + 446,012 + 2,478 = 657,884`) — confirming A/A/A is **NOT** part of s0_Amb (it is s0-ASSIGNED);
4. `sense_Amb` / `anti_Amb` reconstruct from the same cells (185,315 / 315,806, residual 0 each).

**A nonzero residual (>0.5% threshold) ⟹ a hidden fifth bucket — investigate it.** The audit's value is that the small cells are each **NAMED**, not lumped: `A/A/M` (4,742) and `A/M/A` (4,665) are s0-Assigned with one stranded arm Ambiguous (they enter sense_Amb/anti_Amb but **not** s0_Amb or excess); `M/M/M` bilateral (2,478, in s0_Amb, 0 to excess); `A/M/M` the truth-table violation (192, the −1 debit); `S/S/S` the singleton floor (128,065, off-ledger). **There is no fifth bucket** on 0019. The earlier first-order `excess = 2·AP + M` disagreement was **not** a hidden bucket — it was the A/A/A leak (+1) and violations (−1), both folded into identity (2) above, which then closes exactly.

**Denominator discipline (A9 sub-clause; see §3 caution).** Always print numerator and denominator at the **same scale** and **label the scale**. The per-sample-0019 "25.4%" (A/A/A ÷ per-sample excess) and "25.7%" (AP ÷ per-sample s0_Amb) are two different per-sample numerators over two different per-sample denominators; the library-summed **17.5M** s0_Amb (all 45 samples) must **never** denominate a single-sample (0019) numerator (against it the solve-level AP is 23.8%, not 25.7%). A/A/A is **s0-ASSIGNED**, so `A/A/A ÷ s0_Amb` is **not** a real fraction of the ambiguous pile — do not report it. → `FLAG-DENOM-SCALE` (A9).

> **Note (read-level closure of record).** The 0019 read-level closure is the closure of record; a second sample's marginals (0006) are consistent with the same regime structure (`excess ≈ 0.74 · s0_Amb`), but its per-read assignment files were not retained, so its full status-triple cannot be rebuilt from witnessed data — re-run `-R CORE` on its BAM to reproduce.

---

## Resources

- [`METHODOLOGY.md`](./METHODOLOGY.md) — multimapper kernel rationale (Strategy A/B, Teissandier), normalization, combined-matrix validity.
- [`LIMITATIONS.md`](./LIMITATIONS.md) — the 3-rung ladder (§5), the co-oriented-readthrough residue (§2), and the same-strand silent loss limitation (§2a).
- `te-gene-featurecounts` skill — strand-split QC regression: `tests/strand_qc/run_regression.sh` (synthetic truth table in `te-fc:2.0.2` + the `-R CORE` regime-classifier fixture); reference `references/strand-split-qc.md`.
- `te-gene-featurecounts` skill `qc/` — the runnable, parameterized "QC a new TE dataset end-to-end" suite (`qc/run_qc.sh`: BAM dir + SAF + strand + outdir → witness tables + GREEN/RED verdicts). Operationalizes the `-R CORE` regime witness (§4), the kernel-tracer / fragment-weighting check (§9), the closure-table audit (§10), and the `-O` silent-loss attribution (§4a). See `qc/README.md`.
- `featurecounts_results/DE_precheck/` (dataset-side, not in this repo) — the pick-up-ready DE-confound check for GAP (i) / `FLAG-SILENTLOSS-DESIGN`.
- `featurecounts_results/MATRICES.md` (dataset-side, not in this repo) — the five-matrix manifest with per-matrix kernel provenance and do-not flags.
