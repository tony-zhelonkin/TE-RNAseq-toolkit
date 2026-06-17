## Changelog

### 2.0.3 — 2026-06-17 (clarity & scope pass)
Documentation-only pass; no code changes.
- LIMITATIONS: added a legible **escalation map** as the spine — the claim ladder (what a result lets you say, and where the line falls), how far the Random-One kernel and the current SAF reach, the three exits (expand the SAF, change the kernel, go to the wet lab), the design-freedom-vs-baked-pipeline crossover, and a one-screen decision table. Added the gene-anchored-ruler audit caveat to the normalization axis.
- METHODOLOGY: replaced the shared-dispersion claim with "join for size factors, split for dispersion"; added the gene-anchored-ruler audit caveat; scoped Teissandier's "on par" result to the subfamily aggregate; sharpened the design-freedom note.
- QC: moved the operational flag checklist beside the data and rewrote the working principles in plain language.
- README: surfaced the escalation map and the honest design-freedom framing.

### 2.0.2 — 2026-06-15 (kernel / normalization caveat)
Documentation-only pass; no code changes. Added the gene-vs-TE kernel-mismatch caveat, grounded the integer Random-One choice after the strand-split overlap reconciliation, and cross-referenced the strand-split QC. See the METHODOLOGY changelog for detail.

### 2.0.1 — 2026-06-11 (documentation correction)
Documentation-only correction pass; no code changes. Removed two myths and added missing caveats after a claim-by-claim reconciliation against the assembled evidence:
- Corrected the "TEs must be counted unstranded (`-s 0`) because they are bidirectional" claim — for a stranded library + joint gene+TE matrix the more principled best-practice is **stranded TEs matched to the genes** (`-s 2` here; grade B, not a benchmarked standard — see the transparency refinement below); `-s 0` is a defensible standalone choice or for non-directional libraries, and bidirectional biology is preserved by a sense/antisense split, not by collapsing strand.
- Removed the implication that Teissandier 2019 / the Mobile DNA benchmark supports `-s 0`; that benchmark concerns **multimapper handling only** (it remains cited solely for the Random-One/fractional accuracy point).
- Clarified that exon-subtraction (mutual exclusivity) is **necessary but not sufficient**: the combined matrix is valid for within-feature-type across-sample DE only, and only with **gene-anchored size factors**; gene-vs-TE magnitude and "TE %" are not interpretable.
- Added that joint size factors must be **anchored on genes** (`controlGenes` / TMM on the gene submatrix), not pooled across genes+TEs.
- Softened the "lose ~80–90% of TE signal" figure to a young-family-biased fraction (documented multimapper fractions ~10–21%, higher for young L1/SVA/ERV families); the direction (unique-only undercounts) is unchanged.
- Reframed "optimizes strandedness per feature type" as a tolerated interim trade-off (forfeits cross-feature comparability), not a benefit; added a strandedness guidance subsection to the methodology.

**Transparency / evidence-grading refinement (same 2.0.1 pass, no version bump):**
- Corrected the over-claim that "stranded is the field standard." The field is **SPLIT**: the dominant tool TEtranscripts **defaults to UNSTRANDED** (`--stranded no`); best-practice pipelines (TE-Seq, Mobile DNA 2025) recommend stranded. Stranded-for-joint is now presented as **best-practice / mechanistically motivated (grade B), not a benchmarked standard**.
- Flagged the unstranded-sensitivity claim as **unbenchmarked (GAP)**: never measured head-to-head; the only both-mode study (Savytska 2022) found stranded FDR ≤ unstranded — it is a sensitivity-for-specificity trade, not a quantified gain.
- Stated TE bidirectionality as **class-specific** (L1-ASP/ORF0, LTR/ERV established; SINE/intronic weak).
- Added inline **evidence grades** (A/B/C/D/GAP) to load-bearing claims and a new **"Evidence grading & open questions"** section in METHODOLOGY.md (with an "Open gaps / convention-not-evidence" subsection) plus a "Field trajectory / long reads" note.
- Stated plainly that the toolkit's `-s 0` TE default is a defensible standalone choice matching the dominant tool's default, that the in-house `s0/s2 ≈ 0.95` gene result is an empirical in-house measurement (not literature), and that stranded counting + sense/antisense split + genes-only size factors is the best-practice (not proven-superior) joint route — an option with its grade, not a mandate.
