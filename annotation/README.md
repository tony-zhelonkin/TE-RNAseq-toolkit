# annotation

Builds the reference artifacts the TE methodology consumes. Logic lives here and is
git-tagged; the bytes live in a reference cache and carry a manifest pointing back at the
tag that produced them.

## Layers

| Layer | Artifact | Consumed by | Script |
|---|---|---|---|
| **0** sources | genome FASTA, gene GTF, TE GTF | everything | `fetch_mm39.sh` |
| **1** subfamily | grouped, exon-subtracted SAF | gene + TE sense/antisense counting | *pending* |
| **2** genic context | intronic / adjacent / intergenic SAFs + distance table | the intergenic-survival test | *pending* |
| **3** locus | TElocal `.locInd` + strand-reversed twin | locus-level EM | *pending* |

Layer 1 is the production annotation. Layers 2 and 3 re-count or re-align the same data to
attack read-through contamination; see `docs/LIMITATIONS.md`.

## Contig namespace

The genome and gene GTF are `chr`-prefixed. The TE GTF is Ensembl-style unprefixed.
**All derived artifacts normalise to unprefixed.**

This is load-bearing. featureCounts reconciles the prefix silently, so Layer 1 tolerates a
mismatch. `bedtools`, TElocal and Telescope do not — they return zero overlaps and no
warning. Layer 2 and Layer 3 builds therefore assert the namespace rather than assume it.

## Provenance

Every source file is md5-pinned and every derived artifact records the md5s of its inputs
plus the git SHA of the script that produced it. A mismatch is a hard error.

This exists because the previously delivered SAF was exon-subtracted against a filtered
gene GTF from a project output directory that no longer exists, with no record. That build
is not reproducible and must not be reused.
