# Biology primer: transposable elements, strand, and "derepression"

This is the **biology single-source-of-truth** for the toolkit. The other documents reason about
*measurement* — [`METHODOLOGY.md`](./METHODOLOGY.md) (why the counting choices), [`LIMITATIONS.md`](./LIMITATIONS.md)
(what a result can and cannot claim), [`QC.md`](./QC.md) (the QC gate). 

This one explains the *biology* those choices rest on, for a reader who does not already work on transposable elements.

Claims carry the same **A/B/C/D/GAP** evidence grades used across the toolkit: 
* A = peer-reviewed standard, 
* B = single strong source / tool default, 
* C = mechanistic inference, 
* D = folklore / mis-imported, GAP = no adequate source). 

> **Last updated:** 2026-06-11

Citations are listed at the end; textbook fundamentals are left uncited by design.

---

## 1. Strands and transcription direction

DNA is double-stranded and antiparallel; both strands are read 5′→3′. For any given gene, 
RNA polymerase reads one strand (the **template**) and the resulting mRNA matches the other 
(the **coding**, or **sense**, strand). Genes on the `+` and `−` chromosome strands therefore point in
opposite directions — `+`/`−` is just a coordinate label for the chromosome, not where a gene sits.

Two definitions the rest of this document leans on, and they are always **relative to a chosen
feature** (a gene, or a TE):

- **sense** — a read in the *same* orientation as that feature's transcript.
- **antisense** — a read in the *opposite* orientation.

A read is only "sense" or "antisense" *with respect to a particular feature*. The same read can be
sense to one feature and antisense to an overlapping one. That relativity is the root of every
strand subtlety below. **[textbook — grade A.]**

---

## 2. RNA-seq strandedness, and how it reaches the counting flags

Making cDNA normally erases which strand the RNA came from. 
**Stranded** library preps add a step that records it. The dominant protocol is **dUTP second-strand marking** 
(mark the second cDNA strand with uracil, then destroy it before sequencing), 
which yields a **reverse-stranded** library:
* read 2 reports the original transcript's strand (sense), 
* read 1 the opposite. 
**[grade B — Parkhomchuk 2009; Levin 2010.]**

`featureCounts` exposes this as the `-s` flag: 
* `-s 0` unstranded, 
* `-s 1` stranded/forward, 
* `-s 2` reverse-stranded. 

For a reverse dUTP library, 
**`-s 2` counts reads that are SENSE to a feature and `-s 1` counts reads  that are ANTISENSE**. 
**[grade B — Liao 2014; dUTP→`-s 2` convention.]**

That is the whole reason the toolkit counts TEs in three channels: 
* `-s 2` (sense), 
* `-s 1` (antisense), 
* `-s 0` (unstranded, kept only for QC reconciliation). 

What each channel *means biologically* is §4.

---

## 3. What transposable elements are

TEs are mobile, repetitive DNA sequences that make up roughly **45% of the human genome** 
(and a similar share of the mouse genome). 
**[grade A — Lander 2001.]** 
They fall into a few classes:

- **LINEs** (Long Interspersed Nuclear Elements; e.g. **L1**) — **autonomous** retrotransposons: a
  full-length L1 encodes the machinery (ORF1p, ORF2p) to copy itself via an RNA intermediate.
- **SINEs** (Short Interspersed Nuclear Elements; e.g. **Alu** in human, **B1/B2** in mouse) —
  **non-autonomous**; transcribed by Pol III from an internal promoter, and they borrow L1's
  machinery to mobilise. **[grade A — Dewannieux 2003.]**
- **LTR elements / ERVs** (endogenous retroviruses) — retrovirus-derived; their long terminal
  repeats (LTRs) carry promoter and enhancer activity.
- **DNA transposons** — move by a cut-and-paste DNA mechanism; in mammals these are mostly ancient
  and immobile.

Two facts matter for interpretation. 
First, **most TE copies are likely to be old, mutated, and fragmented** — they are genomic fossils, not active elements. 
**Fewer than ~1% are transposition-competent**; in a human, only ~80–100 L1s are active 
and a handful of "hot" L1s account for most activity. 
**[grade A — Brouha 2003.]** 
Second, families differ in **age**: young families (human L1HS, AluY; mouse L1Md, IAP) 
are near-identical across their many copies, which is exactly why short reads cannot tell those
copies apart 
(the locus-resolution limit in [`LIMITATIONS.md`](./LIMITATIONS.md)).

---

## 4. How a TE puts reads into your data: autonomous vs read-through

A TE generates reads two ways, and telling them apart is the central interpretive problem.

- **Autonomous transcription.** The TE's *own* promoter fires (an L1's 5′ promoter, an LTR) and Pol
  II transcribes the TE itself. This is genuine TE expression — what a "derepression" claim is
  really about.
- **Read-through (passive) transcription.** The transcription machinery of a nearby **host gene**
  runs past its normal limits and into a TE inside or beside it, so the reads over the TE are really
  the host's transcript. The TE looks expressed, but it is a passenger.

Where the TE sits decides how exposed it is to read-through:

- **Intergenic** copies have no host gene transcribing through them → their reads can only be
  autonomous → the clean **autonomous-candidate** signal.
- **Intronic** copies sit inside a host gene → host pre-mRNA / intron-retained transcription runs
  through them → their "expression" may be read-through.

Now combine this with strand (§1–2). 
A **sense** read over a TE is 
* either autonomous TE transcription 
* **or** same-strand (**co-oriented**) host read-through
... and at the subfamily level these are **count-identical**, so strand alone cannot separate them. 

An **antisense** read is 
* either opposite-strand host read-through 
* **or** genuine TE antisense transcription (§5). 
This is why matching the gene strandedness (`-s 2`) removes the *antisense* read-through confound but leaves the
*co-oriented* one — the residue that needs genic context 
(Rung 2 in [`LIMITATIONS.md`](./LIMITATIONS.md)), not strand, to resolve. 
**[grade C — mechanistic; the measurement consequences are worked through in LIMITATIONS §2–§4.]**

---

## 5. Antisense and "bidirectional" TE transcription — the real, class-specific story

TEs can be transcribed in both directions, but this is class-specific biology, not a blanket property of all TEs:

- **LINE-1 — a genuine antisense promoter.** The L1 5′UTR carries a sense promoter and a separate
  **antisense promoter (L1-ASP)**, which drives antisense transcription and the **ORF0** product and
  can form chimeric transcripts with neighbouring genes. 
**[grade A existence / B prevalence — Speek 2001; Nigumann 2002; Criscione 2016; Denli 2015.]**

- **LTR / ERV — bidirectional promoter activity.** LTRs can initiate transcription in both
  orientations; documented cases drive divergent host genes and antisense transcripts. 
**[grade A/B — Dunn 2005; Gogvadze 2009; Faulkner 2009.]**

- **SINE / Alu — weak / passive.** Alu have a directional internal Pol III promoter and are
  non-autonomous; only a small minority are detectably Pol III-transcribed, and there is **no
  described autonomous antisense promoter**. Antisense signal over Alu (and over intronic TEs
  generally) is best treated as **host read-through / passenger transcription** unless shown
  otherwise. 
**[grade A non-autonomy; B no antisense promoter; C read-through attribution — Dewannieux 2003; Conti 2015.]**

The one well-quantified bidirectional-initiation finding concerns **TE-associated regulatory elements** 
(DNase-hypersensitive sites) in mouse ES cells, where the majority show balanced
bidirectional initiation, ERV-led and cell-type-specific.
**[grade B — Bakoulis & Andersson 2022.]**

**Takeaway for measurement:** antisense reads over the L1 5′UTR and over LTRs *can* be real
TE-driven biology; antisense reads over Alu and intronic TEs should be read as read-through until
strand-plus-TSS evidence says otherwise. This is why the toolkit keeps the antisense channel as a
*separate* feature set rather than summing it into the sense counts.

---

## 6. Silencing and derepression — why TE expression is interesting

Because active TEs are mutagenic, cells silence them. 
The main layers: 
* **DNA methylation** **[grade A — Walsh & Bestor 1998]**; 
* **H3K9me3** deposited by **SETDB1**, recruited (with **TRIM28/KAP1**) by **KRAB-zinc-finger proteins** 
**[grade A — Rowe 2010; Matsui 2010; Imbeault 2017]**; 
* and, in the germline, the **piRNA pathway** **[grade A — Brennecke 2007]**.

**Derepression** is the loss of that silencing, and TE transcription rising as a result. 
It is observed 
* in **early development**, 
* in **cancer** **[grade A — Burns 2017]**, 
* in **senescence and aging** **[grade B — De Cecco 2019]**, 
* and on treatment with **DNA-methyltransferase inhibitors**.
(A caveat: methylation loss often triggers *re-silencing* by the H3K9 machinery rather than permanent derepression. 
**[grade B — Walter 2016].**)

Why it matters functionally: derepressed TEs — especially ERVs — can form **double-stranded RNA**
that the innate-immune sensors **MDA5/RIG-I** detect as if it were viral, triggering an interferon
response. This **"viral mimicry"** links TE derepression to immune signalling and to the response to
epigenetic therapy. 
**[grade A — Chiappinelli 2015; Roulois 2015.]**

---

## 7. What an aggregate subfamily read-out can and cannot say

This toolkit measures **subfamily-level, strand-resolved differential abundance**. 
Against the biology above, that means:

**It can** 
* detect a net shift in a subfamily's transcription between conditions (a candidate derepression signal), 
* separate it from *antisense* read-through via the strand channel, 
* and serve as a state biomarker that pairs with the silencing-gene and interferon-response read-outs in the same
matrix.

**It cannot**, currently as of June 2026 at least, on its own: 
* prove the signal is **autonomous** rather than co-oriented host read-through (needs genic context — Rung 2); 
* name **which locus** is firing (needs EM — Rung 3);
* prove **bidirectional** TE biology from a sense/antisense ratio alone; 
* demonstrate **dsRNA** or transposition; 
* or say **which silencing lesion** changed. 

Each of those is a separate experiment or a higher rung. 
See [`LIMITATIONS.md`](./LIMITATIONS.md) for where each boundary sits.

---

## References

Grades reflect how the source is used here, not the paper's quality.

- **Parkhomchuk et al. 2009**, *Nucleic Acids Res* — dUTP stranded RNA-seq. doi:10.1093/nar/gkp596. [B]
- **Levin et al. 2010**, *Nat Methods* — comparison of stranded protocols. doi:10.1038/nmeth.1491. [B]
- **Liao, Smyth & Shi 2014**, *Bioinformatics* — featureCounts. doi:10.1093/bioinformatics/btt656. [B]
- **Lander et al. (IHGSC) 2001**, *Nature* — human genome; ~45% TE-derived. doi:10.1038/35057062. [A]
- **Dewannieux et al. 2003**, *Nat Genet* — Alu mobilised by L1, non-autonomous. doi:10.1038/ng1223. [A]
- **Brouha et al. 2003**, *PNAS* — few transposition-competent L1s. doi:10.1073/pnas.0831042100. [A]
- **Speek 2001**, *Mol Cell Biol* — L1 antisense promoter. PMID 11238933. [A]
- **Nigumann et al. 2002**, *Genomics* — L1-ASP chimeric transcripts. doi:10.1006/geno.2002.6758. [B]
- **Criscione et al. 2016**, *BMC Genomics* — genome-wide L1 antisense chimerics. doi:10.1186/s12864-016-2800-5. [B]
- **Denli et al. 2015**, *Cell* — L1 ORF0. doi:10.1016/j.cell.2015.09.025. [A/B]
- **Dunn et al. 2005**, *Gene* — bidirectional LTR (DSCR4/DSCR8). PMID 16288839. [B]
- **Gogvadze et al. 2009**, *J Virol* — intronic antisense LTR effects. doi:10.1128/JVI.00123-09. [B]
- **Faulkner et al. 2009**, *Nat Genet* — FANTOM4 CAGE; TE-driven bidirectional transcription. doi:10.1038/ng.368. [B]
- **Conti et al. 2015**, *Nucleic Acids Res* — few Alu detectably Pol III-transcribed. doi:10.1093/nar/gku1361. [B]
- **Bakoulis & Andersson 2022**, *Nucleic Acids Res* — bidirectional initiation at TE-associated regulatory elements. doi:10.1093/nar/gkac088. [B]
- **Walsh & Bestor 1998** — DNA methylation silences TEs. [A]
- **Rowe et al. 2010**, *Nature* (Trono lab) — TRIM28/KAP1 silencing of ERVs. doi:10.1038/nature08674. [A]
- **Matsui et al. 2010**, *Nature* — SETDB1/ESET silences ERVs. doi:10.1038/nature08858. [A]
- **Imbeault, Helleboid & Trono 2017**, *Nature* — KRAB-ZFP repertoire. doi:10.1038/nature21683. [A]
- **Brennecke et al. 2007**, *Cell* — piRNA pathway. doi:10.1016/j.cell.2007.01.043. [A]
- **Burns 2017**, *Nat Rev Cancer* — TEs in cancer. doi:10.1038/nrc.2017.35. [A]
- **De Cecco et al. 2019**, *Nature* — L1 activation in senescence drives interferon. doi:10.1038/s41586-018-0784-9. [B]
- **Walter et al. 2016**, *eLife* — re-silencing after methylation loss. doi:10.7554/eLife.11418. [B]
- **Chiappinelli et al. 2015**, *Cell* — ERV dsRNA / viral mimicry. doi:10.1016/j.cell.2015.07.011. [A]
- **Roulois et al. 2015**, *Cell* — dsRNA from ERVs after DNMT inhibition. doi:10.1016/j.cell.2015.07.056. [A]
