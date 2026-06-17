# Limitations: What This Instrument Can and Cannot Claim

This is the **boundary document**. 
It names the inferential limits of the TE quantification this toolkit performs 
and where a result stops being supported by the current design. 
For the *biology* these limits rest on — strand, read-through, bidirectional transcription, derepression — see [`BIOLOGY.md`](./BIOLOGY.md). 
For the *why* behind the design choices (strandedness, multimapper kernel, joint matrix, normalization), 
and for the evidence grades and citations those choices rest on, see [`METHODOLOGY.md`](./METHODOLOGY.md). 

> **Last updated:** 2026-06-17

Grades below use that document's **A/B/C/D/GAP** scale and its citations are not restated here.

## The instrument, in one line

featureCounts on a grouped, exon-subtracted SAF (Simplified Annotation Format) where 
* `GeneID = Subfamily:Family:Class` (the full set of mouse subfamilies), 
* STAR Random-One + integer `-M` (no `--fraction`), 
* counted in three strand channels — sense `-s 2` (gene-matched, the DE numerator), antisense `-s 1` (separate channel), unstranded `-s 0` (QC reconciliation only, never a DE numerator). 

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

## The escalation map: what you can claim, and how to reach further

<a id="escalation-map"></a>

This is the spine of the document, and the part I reach for first when someone asks "so what does this number *mean*?"
Read this before the technical sections below — they justify what this map asserts.
The rule throughout: **lead with the action, then the one-line reason.**

The map answers three questions in order.
What does a result let me say (the claim ladder)?
How far do my current tools take me, and where do they stop (the kernel, the SAF)?
When I need more, which of the three exits do I take, and what does each actually resolve?

### A. The claim ladder — what a result lets you say
<a id="claim-ladder"></a>

A stranded sense hit on a subfamily lets you make a *stack* of claims, each stronger than the last.
The discipline is simple: start at the top, climb only as far as your checks support, and **stop at the first claim you cannot back.**
Each step states the claim in plain words, then names what you would need to go one step further.

This is the **claim ladder** — what a result entitles you to *say*.
It is distinct from the **method ladder** (the three rungs in §5), which is about the *tools and resolution* you bring to bear.
Claims here; methods there. They climb in parallel, but do not confuse one for the other.

Safest first:

1. **"A subfamily's sense-strand read count rose, reproducibly, across my samples."**
   This is the literal output — the same feature, counted the same way, with size factors anchored on genes, so the per-feature biases cancel in the ratio.
   *To go past this you need:* nothing — this is what the instrument hands you.

2. **"That subfamily's sense-strand RNA actually went up."**
   The count rise reflects a real rise in steady-state [sense-strand](./BIOLOGY.md) RNA, not a counting artifact, because each fragment contributes exactly one count.
   *To go past this you need:* the antisense channel checked (next step).

3. **"The rise is strand-specific."**
   It shows up on the sense strand and not as a matching antisense rise.
   *To go past this you need:* the [antisense](./BIOLOGY.md) channel actually tested — a condition-by-channel interaction. You cannot claim this if you never looked at antisense.

4. **"This subfamily is coupled to the pathway I'm studying."**
   Its pattern across conditions tracks the gene markers of that pathway (for an immune contrast, the interferon-stimulated genes in the same matrix).
   This is **coupling, not causation** — the subfamily moves *with* the pathway, which does not mean it drives it.
   *To go past this you need:* the relevant gene markers to have moved with it, as a positive control — and then you are at the line.

*— the line falls here. Everything above is licensed by this toolkit. Everything below needs more. —*

5. **"This subfamily is autonomously derepressed — its own promoter fired."**
   *You cannot claim this from a stranded count.* Same-strand host read-through looks identical to autonomous transcription on the sense strand.
   *To reach it you need:* **Rung 2** — genic context, which tests whether the signal survives in copies that have no host gene to read through. See §5 and exit D.1 below.

6. **"This specific copy, at this address, is the one firing."**
   *You cannot claim this.* The kernel discards locus identity by design — it tells you the subfamily moved, never which copy.
   *To reach it you need:* **Rung 3** — an expectation-maximization tool (Telescope, SQuIRE). See §5 and exit D.2.

7. **"This TE drives the phenotype."**
   *You cannot claim this from any steady-state count.* Causation is a perturbation question, and "did the promoter fire" is a nascent-transcription question.
   *To reach it you need:* the **wet lab** — a perturbation or a nascent-transcription assay. See exit D.3.

**The young-vs-old caveat — read this before you quote a fold-change.**
For **old, uniquely-mapping families**, the fold-change is trustworthy in **both direction and size** — quote it as a number.
For **young, high-copy families** (L1Md, IAP, active ERVK), only the **direction** is safe — say "up" or "down," not "up 3-fold."
The reason: a young-family read maps to many near-identical copies, and the kernel places it at one at random. That scatter cancels in the subfamily total *on average*, so the direction is robust, but it adds counting noise to the *size* at any single run. Direction survives; magnitude carries multimapper noise.
**(grade A for old families; grade B for the young-family direction-only call — see the METHODOLOGY grading table.)**

### B. How far the current kernel takes you
<a id="kernel-reach"></a>

The kernel is STAR Random-One: each fragment is placed at one of its valid locations, chosen at random, and counted as one integer.

**Good for:** subfamily-level direction. The random placement averages out across a subfamily, so the *direction* of a subfamily shift is trustworthy down into the young families.

**Stops at:**
- **which locus** — it discards locus identity by design, so it can never tell you which copy fired;
- **magnitude for young, high-copy families** — the random placement adds counting noise to the size (not the direction) for families with many near-identical copies.

**You change the kernel** (to an expectation-maximization tool) when the question is **"which locus,"** or when you need an **autonomy call on a young family** — neither of which Random-One can answer. That is exit D.2.

### C. How far the current SAF takes you
<a id="saf-reach"></a>

The SAF is grouped, exon-subtracted, and keyed by subfamily.

**Good for:** subfamily abundance — how much of each subfamily's RNA is present, across samples.

**Stops at:** the autonomy question. It cannot separate a TE transcribed from its **own** promoter from a TE that is just caught in a host gene's **read-through**. That separation needs **genic context** — where each copy sits relative to genes — and genic context lives in the **annotation** (Rung 2, exit D.1), not in a different counting tool. You expand the SAF; you do not swap the kernel.

### D. The three ways to reach further
<a id="three-exits"></a>

When the claim you need sits below the line, there are exactly three exits. Each does some things and not others — pick by what you actually need to resolve.

#### D.1 — Expand the SAF (Rung 2)
<a id="exit-saf"></a>

**When:** you need an **autonomy call**, or you suspect read-through.

**How (plainly):** tag each TE copy by its genic context — intergenic, intronic, or gene-adjacent (treat a copy within ~10 kb of a gene as adjacent, and require ≥20 kb from any gene for a *clean* intergenic set). Then **re-count the same alignments** through the tagged annotation. No realignment — the BAMs do not change, only how you bin them.

**Resolves:** for **old, uniquely-mapping families**, whether the signal survives in the host-free intergenic copies (an **autonomous candidate**) or shows up only in intronic copies that correlate with their host gene (**read-through**).

**Does NOT resolve:**
- **young families** — the multimapper scatter blurs the intergenic/intronic split, so hand these off to Rung 3 (exit D.2);
- **interferon-inducible enhancer-RNA LTRs** — these survive the intergenic test yet are not autonomous retrotransposition (the enhancer is firing, not the element's own promoter);
- **distal read-through** — host transcription can run tens of kilobases past a gene under stress, reaching past the adjacency window.

This is **Rung 2** of the method ladder (§5). The full intergenic-survival logic and its conditions live there.

#### D.2 — Change the kernel / use a baked pipeline
<a id="exit-kernel"></a>

**When:** **"which copy"** questions, or a **young-family autonomy** call.

**Tools:** an expectation-maximization kernel (Telescope, SQuIRE) or a bundled pipeline (TE-Seq).

**What you get:** **locus-level calls** — claims about individual copies, which Random-One cannot make.

**The cost:**
- locus-level false-discovery is **high for everyone**, and worst for young families (see the METHODOLOGY field-trajectory note);
- the bundled tools **fix you to a two-group design** — you give up the flexible statistical model this toolkit lets you write.

This is **Rung 3** (§5).

#### D.3 — Go to the wet lab
<a id="exit-wetlab"></a>

**When:** the question is whether a promoter **actually fired**, whether the TE is **causal**, or **which copy is active in vivo** — none of which a steady-state count can answer, at any rung.

**Which assay:**
- **"did the promoter fire"** — nascent transcription (PRO-cap / GRO-seq / TT-seq);
- **"is this specific element expressed"** — targeted RT-qPCR, 5′-RACE, or ORF1p protein;
- **"is it causal"** — a perturbation (knockdown, inhibitor, or an earlier timepoint).

### E. Design freedom vs a baked pipeline — when freedom becomes a burden
<a id="design-freedom"></a>

**The premise.** This toolkit hands you a clean integer matrix and stops, so **you own the statistical design.** You can fit a factorial model, an interaction, a blocked or repeated-measures design — whatever your experiment needs. The baked pipelines cannot: they bake in a fixed two-group test.
**(grade A — design freedom; see the METHODOLOGY grading table and §design-freedom.)**

**The honest turn.** That freedom becomes a *burden* at one specific point: when answering your question would force you to rebuild a baked pipeline's machinery by hand — locus-level expectation-maximization *plus* genic-context stratification — just to keep your flexible design alive. At that point the freedom is costing you more than it buys.

**What to do there.** Use the baked pipeline and accept its design, or run both and reconcile. Do not invent one-off tricks to preserve flexibility you no longer need.

**The crossover, named:**
- a **complex design** on a **subfamily / strand** question → **this toolkit** (its design freedom is the whole point);
- a **standard locus-level autonomy** call → a **baked pipeline** (let it own the design; you were going to rebuild it anyway).

### F. When to skip computation and go straight to the wet lab
<a id="skip-to-wetlab"></a>

If your question is fundamentally **did the promoter fire**, **is the TE causal**, or **which copy is active in vivo**, then climbing rungs will not close it — no steady-state count can.
The assay is **faster and more definitive** than the computation here.
Say so to yourself early, so you do not spend weeks on convoluted analysis for a question that a single experiment answers cleanly. The rungs are for the claims *above* the line; these questions live below it.

### G. The decision table
<a id="decision-table"></a>

One screen, to find your question and read off where it lands.

| The question you're asking | Furthest tool that answers it | What you still cannot claim |
|---|---|---|
| Did a subfamily's abundance shift? | This toolkit (Rung 1) | That it's strand-specific, coupled, or autonomous |
| Is the shift strand-specific? | This toolkit, antisense channel tested (Rung 1) | That it's autonomous vs co-oriented read-through |
| Is the subfamily coupled to my pathway? | This toolkit + the gene markers in the same matrix | Causation — coupling is not driving |
| Is an **old** family autonomously derepressed? | Rung 2 (genic-context SAF) | That a *specific copy* fired; that the promoter fired |
| Is a **young** family autonomously derepressed? | Rung 3 (EM tool) — Rung 2 blurs here | A clean locus call (locus FDR is high for young families) |
| Which copy / locus is firing? | Rung 3 (Telescope / SQuIRE) | A low-FDR answer; in-vivo activity |
| Does the TE drive the phenotype? | Wet lab (perturbation / nascent assay) | Anything from steady-state counts alone |

---

## 1. Matching the strand basis is NECESSARY, not SUFFICIENT

Counting TEs at the gene strandedness (`-s 2` sense, matched to the `-s 2` gene denominator) buys two **distinct** things 
that happened to ride on the same `-s 2` flag. It's important not to conflate them.

### (a) MEASUREMENT axis — a valid ruler for normalization. **SURVIVES.**

A per-sample size factor `s_j` says "sample *j*'s gene counts are inflated/deflated by this much, in the units genes are measured in." 
Dividing a TE row by `s_j` **asserts** the TE row is measured in those same units — that genes and TEs accrue counts from true library depth 
at the same per-sample rate. Matching the strand basis seems, at least to me, internally defensible. 

In a stranded derepression study, normalize the TE numerator on the same strand basis as the gene denominator — keep an unstranded TE numerator off a stranded gene denominator. 
The logic is: are you even using the right *measurement ruler*? The way the size factor was computed determines if you can apply the ruler to the TE rows.
**(grade B — matched-basis best-practice; see METHODOLOGY §Strandedness and §Combined.)**

One caveat on this ruler, and it matters most in exactly the strong-perturbation studies I run. 
It assumes the gene compartment is mostly stable. 
Under a transcriptome-wide response — an interferon/ISG program is the case in point — the gene size factors can start to move *with* the contrast, 
and then a TE fold-change is read against a ruler that is itself bending with the biology. 
So before I trust TE fold-changes in such a contrast, I check the gene-derived size factors are not themselves correlated with it 
(and that they agree with a housekeeping subset); if they bend, I anchor on an explicit stable-control set instead. 
The gene rows stay a valid *observed* control even when they are no longer a valid *ruler*. 
**(grade B — see METHODOLOGY §Normalization.)**

### (b) BIOLOGY axis — the autonomous-vs-read-through verdict. **DOWNGRADED.**

At subfamily resolution the antisense split gives strand **asymmetry**, autonomy verdict is a separate stage. 
The `-s 1` channel is a strand-of-read aggregate over a whole subfamily.
It may become a verdict on the autonomous vs read-through TE element only when paired with 
- the genic-context 
- and locus axes 

— and the grouped SAF method I usually employ for subfamily aggregate resolution deliberately discards those (§5). 
**(grade C — inference; the verdict-grade form is TE-Seq's three-axis design.)**

The whole point: matching the strand basis earns you (a) outright and only a *signal* (asymmetry) on (b) — never the *verdict*.

---

## 2. What matching REMOVES vs what it leaves (the residue)

**First, what "read-through" means.** A TE puts reads into your data two ways.
* *Autonomously* — its own promoter (an L1's 5′ promoter, an LTR) fires and Pol II transcribes the TE itself; 
this is genuine TE expression.
* *Passively* — the transcription machinery of a nearby **host gene** runs past its normal limits 
and into a TE sitting inside or beside it, so the reads over the TE are really the host's transcript, 
not the TE's promoter. 
That second case is **read-through**, and §2–§5 are about telling the two apart.

**Exon subtraction is double-counting hygiene, not autonomy isolation.**
The SAF is built by subtracting gene *exons* from the TE intervals.

That keeps host *exon* (mature-mRNA) expression from being miscounted as TE and makes genes and TEs mutually exclusive features —
but it does **not** select for autonomous TEs, and it is orthogonal to read-through.

Read-through is mostly *intronic*, and intronic TEs are not exonic, so they survive the subtraction:
a substantial share of the retained TE territory still sits inside gene bodies (read-through-exposed), 
the rest is intergenic (the cleaner autonomous-candidate territory — cleaner for old, uniquely-mapping families, ≥20 kb from any gene and strand-checked, §5).

Separating the two is genic-context stratification — Rung 2 (§5) — which exon subtraction cannot reach.
**(grade A — the field separates autonomy by genic-context stratification, e.g. TE-Seq's Exonic/Intronic/Intergenic labelling, not by exon subtraction.)**

**Matching REMOVES:**
- antisense-oriented read-through (host transcription on the opposite strand) — it lands in the `-s 1` channel,
out of the `-s 2` numerator;
- the normalization-rate mismatch — the ruler problem mentioned above.

**Matching DOES NOT remove — the residue:**
- **co-oriented (sense) read-through.** A TE sitting sense within a sense intron dumps host transcription 
straight into the `-s 2` channel. 

**No strand flag separates that at subfamily level** — the read is on the "right" strand. 
This is the locus/context problem, and the current grouped `Subfamily:Family:Class` SAF does not touch it. 
**(grade C — inference; the unsolved gene–TE disambiguation gap, METHODOLOGY §Open gaps #3.)**

Matching the strand basis may be only a partial read-through defense.
The part it cannot reach is exactly the part separated by **genic context** — does the apparent signal survive 
in the read-through-free **intergenic** copies? 

That is Rung 2 of the ladder, which I come back to in §5.

---

## 2a. Same-strand overlap silent loss (strand-split overlap-ambiguity)

This is another potential limitation, **distinct from and alongside** the one described above.

**Silent loss is an annotation-geometry tax** (same-strand nested overlap dropped to Ambiguity); 
**co-oriented read-through is a transcription-context confound** (host transcription on the right strand). 

They are different failure modes; closing one does not close the other, both unresolved at subfamily resolution. 

**What it is.** A read overlapping two **same-strand** features is:
* dropped as Ambiguity by `-s 0` 
* *and* by both stranded passes (the strand-matched pass still sees ≥2 candidates → Ambiguity; 
the opposite pass sees none → NoFeatures). 
It contributes 0 to the unstranded, sense, and antisense counts alike. 
The toolkit's three-channel counting never sees it.

**The three fates of a strand-ambiguous read.** When the unstranded pass (`-s 0`) cannot 
uniquely assign a read because it overlaps two annotated features at once, re-counting that read *by strand* 
sends it to exactly one of three fates:
- **AP — antiparallel.** The two features are on **opposite strands**; the sense pass keeps one and the antisense pass keeps the other, 
so the read is counted **twice, under two different subfamilies** ("double-presence"). 
This is exactly what would be double-counted if the sense and antisense matrices were ever summed.
- **M — asymmetric (genuine reclaim).** Splitting by strand **breaks the tie** — only one of the two features is 
on the matched strand — so the read becomes **uniquely assignable** and is legitimately **recovered**.
- **P — parallel, same-strand (silent loss).** Both features are on the **same strand**, so no strand flag can separate them; 
the read is dropped by **all three** passes and never appears in any count — the case described just above. 
**Residual-invisible**, and the **dominant** fate.

**The witnessed ledger.** A per-read audit of the strand-ambiguous reads (on internal data) broke down 
two ways — and **the two groupings below use different denominators, so my recommendation is not to mix them**: 
* *Share of the ambiguous pile* (which fate befell an ambiguous read): 
    - silent loss (P) is **≈ two-thirds**, 
    - antiparallel (AP) **≈ a quarter**, 
    - genuine reclaim (M) the small remainder — silent loss dominates; 
* *Share of total assigned TE signal* (what actually reaches DE): 
    - silent loss (P) ≈ **9–10%**; 
* antiparallel double-presence (AP) a few percent; 
* genuine reclaim (M) about a percent; 
* the rare "bilateral" case (a read overlapping ≥2 features on *both* strands) negligible 
(so the single-strand-drop assumption holds). 

Silent loss is the **dominant** fate of the ambiguous pile. 
**(grade A — per-read witness, reproduced on internal data.)**

**The conditional-losslessness principle.** The `s0 − sense − anti` residual is 
lossless **only relative to s0's *assigned* set**. 

Same-strand silent loss is residual-invisible by construction: 
it is 0 in all three channels and carries no feature assignment, so the per-read audit records it as unassigned. 

The residual can therefore see only the antiparallel and asymmetric recovery slice, 
and that slice sits on top of an unseen silent-loss pile roughly twice its size. 

So the residual being one-sided tells you about recovery, and says nothing about silent loss. 
(The full statement of this principle lives in [`QC.md`](./QC.md).)

**The operational principles live in [`QC.md`](./QC.md).** 
The working principles, the directional meter, 
and the gate thresholds are owned there; they are not duplicated here. 
This section names the limitation and its boundary. The four cases this becomes **DANGEROUS** in — each stated as a principle in [`QC.md`](./QC.md) — are:
- **summing channels** — `(sense + anti)` double-counts every read that lands in **both** channels under different subfamilies — the antiparallel (AP) double-presence pile; never use the channel sum as a quantity, and do not read it as bidirectional signal;
- **s0-denominating the ERV/satellite/DNA tail** — s0 (unstranded) drops their antiparallel reads to Ambiguity, a pathological denominator; denominate on the sense channel for these families;
- **applying gene size factors to the `-M` TE matrix** / comparing gene-vs-TE magnitude — genes and TEs are counted on different kernels, so this comparison is not like-for-like (the kernel-mismatch caveat; see [`QC.md` §6](./QC.md), METHODOLOGY §Combined);
- **silent loss aligning with the design axis** — sample-stable loss cancels in cross-sample DE, but the silent-loss share co-varies with the net strand-offset, so check it against the design matrix before declaring it condition-independent.

**Scope.** The young autonomous candidates (L1Md_T/Gf/A, IAPEz) sit **🟢 GREEN: well past the conservation gate (>90%)**. 
The burden sits on the **old SINE / MaLR-LTR** subfamilies (silent loss) and on **ERVK-LTR / satellite** (antiparallel double-presence), 
not on the young autonomous tail. 
See [`QC.md` §7](./QC.md) for the GREEN/RED gate. 
The young-silent attribution is **witnessed for the gate set** via the `-O` target-revealer ([`QC.md` §4a/§9](./QC.md)): 
the young silent-loss share is a small fraction of a percent, concordant with the geometry prediction, 
so the GREEN gate's **second leg is closed at read level** (no longer geometry-inference alone) — the genic-context 
Rung-2 GAP (intron/intergenic stratification on the richer SAF, §5) stays open. 
**(grade A — `-O` per-read young-silent witness, concordant with geometry; the genic-context attribution remains the Rung-2 GAP, §5.)**

**The YELLOW DE-confound.** Sample-to-sample the silent-loss share varies less than the net strand-offset does, 
so it cancels in cross-sample contrasts under the usual conditions. 
It does, however, co-vary with per-sample net strand-offset — silent loss is not a flat tax, it tracks sample-level strand behaviour. 
So it **must** be checked against the experimental design matrix before treating it as condition-independent. 
**(grade B — sample-stable; YELLOW pending internal matrix design-checks.)**

---

## 3. The worked case: how unstranded counting `-s 0` can potentially fabricate derepression

One TE subfamily. No real change in autonomous transcription. 
Ctrl vs Senescent, equal depth, so gene size factors `s_j ≈ 1.0` in both.

| Reads accrued | Ctrl | Sen |
|---|---|---|
| own-strand (autonomous, **unchanged**) | 100 | 100 |
| antisense read-through (**↑ in senescence**) | 20 | 60 |

- **`-s 2`** (matched, sense only): `100 → 100`. Normalized `log2FC = log2(100/100) = 0`. **Truth preserved.**
- **`-s 0`** (mismatched, both strands): `120 → 160`. Normalized `log2FC = log2(160/120) = `**`+0.42`** (≈ 1.3× **fabricated** "derepression").

The gene-anchored size factor **cannot rescue** the `-s 0` case: genes never saw the read-through bump (`-t exon`, on-strand), 
so `s_j` stays ≈ 1 and leaves the inflation in place. 
Worse, the fake signal points the **same direction** as the hypothesis (senescence → TE up), so it reads as confirmation.

This only bites when **strand-capture varies across samples** — which is exactly the derepression setting: 
intron retention and read-through rise with senescence / stress / immune activation, 
correlated with the tested biology. 
A *constant* rate difference would cancel in a log-fold-change; it is the *covariance with condition* that 
fabricates signal. 
**(grade C — inference, arithmetic; the disambiguation problem is METHODOLOGY §Open gaps #3, the strand-basis rationale §Strandedness.)**

---

## 4. The antisense channel must DO WORK

Holding stranded counts on disk is only the first step. 
Subfamily-level featureCounts + the stranded sense/antisense design is the **right instrument** — 
but the counts have to be *used*, and it matters precisely what each channel can and cannot settle.

**What the antisense channel settles.** A condition × channel interaction test (does the sense/antisense 
ratio shift with condition?) catches the **antisense flavour** of read-through — opposite-strand host 
transcription that lifts antisense while sense stays flat — and confirms a **bidirectional-class** signal 
when sense and antisense co-move (real for L1-ASP/ORF0 and LTR/ERV — grade C class-specific). 
Run library-wide, a broad antisense rise with condition is a **global read-through sentinel**: read-through 
pressure is up genome-wide, so the sense calls inherit suspicion.

**What it cannot settle — and this is the load-bearing point.** 
Co-oriented (same-strand) read-through and genuine autonomous derepression are **both sense-up, antisense-flat**. 
They produce the *same* ratio, so the interaction test cannot separate them, and it is 
silent on bidirectional autonomous derepression (sense and antisense rise together, ratio flat). 

Matching the strand basis and testing the channel interaction is therefore **necessary, not sufficient**: 
it removes the antisense flavour, raises the prior, and flags asymmetry — it does **not** certify a 
sense-up call is autonomous rather than co-oriented read-through.

**Where the verdict comes from.** The co-oriented residue is separated by **genic context, not by strand** — 
Rung 2 of the ladder (§5): does the apparent derepression survive in the read-through-free **intergenic** 
copies? That is the hand-off. The antisense channel does its work; the autonomy verdict waits on Rung 2.

A `-s 2` sense-channel result simply *called* "derepression" without the antisense channel doing this work, 
and without the Rung-2 hand-off for any sense-up / antisense-flat call, is still exposed to every residue 
in §2–§3. 
**(grade B — explicit sense/antisense channel is TE-Seq's design, SQuIRE per-locus strand as precedent, 
METHODOLOGY §Strandedness; the interaction test as a necessary-not-sufficient sentinel — separation 
deferred to genic context — is grade C inference.)**

---

## 5. The 3-rung ladder (resolution vs robustness)
<a id="method-ladder"></a>

This is the **method ladder** — the tools and the resolution they buy.
It is distinct from the **claim ladder** at the top of the document (the [escalation map](#claim-ladder)), which is about what a result entitles you to *say*.
Methods here; claims there. The claim ladder tells you when you are forced onto a higher rung; this ladder tells you what each rung costs and delivers.

Each rung adds an axis. Higher rung = more specific claim, more fragile inference.

- **Rung 1 — stranded sense/antisense, grouped SAF.** *(← what this toolkit does.)*
  Strand asymmetry and the antisense-flavour / bidirectional read-through, at **subfamily level, statistically**, via the condition × channel interaction (§4). 
The co-oriented residue is **not** settled here — it hands off to Rung 2. 
featureCounts + grouped SAF. **ROBUST; current-standard genre.**

- **Rung 2 — + genic-context-stratified SAF (intergenic / intronic / adjacent).**
  Read-through attribution **structurally**: split each subfamily into e.g. `L1HS:intergenic`, `L1HS:intronic`. 
**Important:** Genic context lives in the SAF **annotation**, not a different downstream tool. 

Intergenic loci have no host to read through, so for old, uniquely-mapping families — at least ~20 kb from any gene and strand-checked, 
and where the locus is not an active enhancer LTR — intergenic-TE derepression is the **cleaner** autonomous-candidate signal; 
Intronic "derepression" *not* mirrored intergenically is a **read-through flag**. 

Still can be done with featureCounts, but needs richer SAF. 
This is TE-Seq's design — it labels every TE locus Exonic/Intronic/Intergenic (plus gene-adjacent and UTR) 
and stratifies its DE on that axis — so the rung is grounded in a real pipeline, not only inference.

- **Rung 3 — EM locus assignment.**
  Names the firing copy (which L1HS at which locus); per-locus context + neighbor-gene correlation. 
Telescope (TE-Seq wraps it) / SQuIRE. **SPECIFIC but FRAGILE** 
(per-locus FDR — METHODOLOGY §Field trajectory: TElocal ≈ 26% FDR, Savytska 2022).

**Why young families hand off from Rung 2 to Rung 3 — the two things that leak into a stranded TE count.**
Two different things put unwanted reads into a sense channel, and they have different fixes.

The first is **real read-through biology** — a host gene's transcription depositing RNA over intronic copies. 
Genic-context stratification (Rung 2) attacks this one directly, because intergenic copies have no host to read through: if the signal survives in the host-free intergenic stratum, it was not read-through.

The second is the **multimapper placement itself**. A young-family read has many equally-good positions, and the kernel places it at one of them. A context-stratified re-count uses the **same alignments** — it inherits that placement and cannot move the read back to where it really came from. 

So the intergenic test is a **clean** test for **old, uniquely-mapping families** (no placement ambiguity to inherit) and a **blurred** one for **young, high-copy families** (the placement noise leaks across strata). 
That is exactly why young families hand off to EM at Rung 3: only re-placing the reads — which Rung 2's same-alignment re-count cannot do — can clean the split. 
**(grade A for the two-channel mechanism; the exact young-family leakage severity is a GAP — see the METHODOLOGY grading table.)**

**Judgment.** For an **aggregate-derepression** claim, 
finish Rung 1 properly (the interaction) and 
reach for **Rung 2 before Rung 3** — a richer SAF buys most of TE-Seq's read-through defense 
at a fraction of the cost and none of the locus-FDR fragility. 

Climb to Rung 3 when the question is genuinely **"WHICH insertion,"** because that is where 
you trade a robust claim for a specific one — or for a **young-family autonomy** call, where Rung 2's same-alignment re-count cannot un-blur the strata (above). 

For any **load-bearing sense-up / antisense-flat** call — and for every call once 
the global read-through sentinel fires — Rung 2 is **required**, not optional: it is the only step 
that separates co-oriented read-through from autonomous derepression. 

Stop at Rung 1 only for an honest strand-asymmetric signal, never an autonomy verdict.

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

TE-Seq pairs all three (strand: sense/antisense; genic context: intergenic/intronic/adjacent; locus: Telescope EM) 
to call autonomous-vs-read-through. 
**This toolkit has only the strand axis.** So it earns the signal, and the verdict and address stay out of reach. 
**(grade A for the genic-context rung — TE-Seq's Exonic/Intronic/Intergenic stratifier; the rung-ordering judgment is grade C inference.)**

---

## 6. Magnitude caveat

Do **not** promote TE% or gene-vs-TE magnitude to biology. 
Making the joint matrix normalization-**coherent** (mutual exclusivity via exon subtraction, genes-only size factors) 
is the most seductive moment to over-read it — **coherence is not magnitude.**

- No TPM / FPKM for TE meta-features (summed loci have no single length — METHODOLOGY §Open gaps #5).
**(grade C — inference.)**
- "TE % of transcriptome" is a QC sanity band, not a biological fraction (METHODOLOGY §Normalization, §Combined).
- The joint matrix is valid for **within-feature-type, across-sample DE only**, 
never gene-vs-TE within-sample magnitude.

The shared strand basis earned two things and only two: strand **asymmetry** and a valid normalization **ruler**. 
The autonomy verdict and the magnitude statement stay beyond it. 
**(grade C — the no-cross-feature-magnitude / no-TPM rule, METHODOLOGY §Combined and §Open gaps #5.)**

---

## One-line boundary

This instrument measures **subfamily-level, sense-channel, depth-normalized differential abundance**, 
and — when the antisense channel is made to work (§4) — **strand asymmetry**. 

It does not, on its own, certify **autonomy**, **read-through-free** signal, **locus identity**, or **magnitude**, 
and it does not recover **same-strand-overlap silent-loss reads** (residual-invisible, §2a). 
Each of those requires climbing the ladder (§5).
