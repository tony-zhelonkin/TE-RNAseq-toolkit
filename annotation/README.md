# annotation

**The TE annotation build lives in its own repo: https://github.com/Mogilenko-Lab/te-tracks**

`te-tracks` builds reproducible TE annotation tracks from RepeatMasker. One canonical table holds
one row per TE locus, and every delivered track is a projection of it, so bulk RNA-seq, single-cell
RNA and single-cell ATAC share one subfamily universe and one locus identity.

This directory holds the pointer. The science stays here, in
`docs/{METHODOLOGY,BIOLOGY,LIMITATIONS,QC}.md`.

## Consuming the tracks

Reference roots carry a dated snapshot and a `current` symlink:

    /data2/users/shared/refcache/te_mm39/current          local
    /gpfs/data/rathmell-lab/data/refdata/te_mm39/current  CRI

`MANIFEST.json` in each snapshot records the source md5s, the output md5s and the builder git SHA.
`te-tracks/docs/OUTPUTS.md` states every file, its columns and its coordinate convention.
