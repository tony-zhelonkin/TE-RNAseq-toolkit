# TE-RNAseq-toolkit


## Rough notes

Prepare STAR parameters 

**About STAR options for TEs**  
We’ll pass them **from CLI** via `--extra_star_align_args` (this parameter is officially supported in the STAR-Salmon route). Use:

```bash
--extra_star_align_args '--twopassMode Basic \
  --alignEndsType EndToEnd \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped None \
  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMismatchNmax 999 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
  --outFilterScoreMinOverLread 0.4 \
  --outFilterMatchNminOverLread 0.4 \
  --outFilterType BySJout \
  --outFilterMultimapNmax 100 \
  --winAnchorMultimapNmax 200 \
  --outMultimapperOrder Random \
  --outSAMmultNmax 1 \
  --outSAMprimaryFlag OneBestScore \
  --runRNGseed 777'
```


Brief interpretation
```bash
--twopassMode Basic                    # discover novel SJs in pass1; improves splice accuracy
--alignEndsType EndToEnd               # no soft-clipping; reduces opportunistic partial hits in repeats
--outSAMunmapped None                  # don’t emit unmapped reads (smaller BAM; QC modules won’t use them)
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1  # reasonable SJ support (discovery + annotated)
--outFilterMismatchNoverLmax 0.06      # ≤6% mismatches per read; keeps alignments clean
--outFilterMismatchNmax 999            # (paired with the relative cap above)
--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000  # mammalian intron/gap ranges
--outFilterScoreMinOverLread 0.4       # ≥40% of read must be mappable by score
--outFilterMatchNminOverLread 0.4      # ≥40% of read must match; avoids tiny anchors in repeats
--outFilterType BySJout                # prefer alignments consistent with known/learned SJs; fewer spurious SJs
--outFilterMultimapNmax 100            # retain multi-mappers (needed for TE signal)
--winAnchorMultimapNmax 200            # allow many candidate sites in repetitive sequence
--outMultimapperOrder Random           # randomize among equally good hits
--outSAMmultNmax 1                     # **emit ONE alignment per read** (random-one strategy)
--outSAMprimaryFlag OneBestScore       # never write multiple primaries on ties
--runRNGseed 777                       # reproducible random selection across runs
```

Tell STAR to **emit a genome BAM**  
Add back: `--outSAMtype BAM Unsorted`

#### STAR parameter and different annotation tracks
##### 1) **Exons (gene-level)**

- **Good**: `--twopassMode Basic` + `--outFilterType BySJout` improves splice junction calling and reduces mis-spliced alignments.
- **EndToEnd** ensures reads won’t be soft-clipped to fit repeats; if adapters/low-quality tails remain, they simply fail (nf-core trims by default, so you’re fine).
- **Multi-mappers**: harmless for exons because you emit **one** best/random alignment per read. You’ll get consistent, integer counts with featureCounts (no `-M`, no `--fraction`).
    

##### 2) **Introns**

- If you derive introns as **genebody − exons** (recommended) or from an intron SAF, this mapping favors clean separation of spliced (exonic) vs unspliced (intronic) signal.
- Two-pass + SJ filtering reduces the chance that intronic reads get “forced” into dubious splice alignments.
    

##### 3) **Genebody**

- Broad first-to-last exon spans benefit from permissive intron windows, so reads across large introns remain alignable. Your `alignIntronMax 1e6` / `MatesGapMax 1e6` cover that.
    

##### 4) **TE (family/subfamily/class; integer counts)**
- This is exactly the **random-one** policy: keep many candidate loci (`Nmax 100`) but **write one** per read at random (`outSAMmultNmax 1` + `OneBestScore` + `Random` + fixed seed).
- That harmonizes with featureCounts **without** `--fraction` to produce **integer** family/subfamily counts, avoiding multi-mapping inflation.
- Caveat: **don’t** use this BAM for **locus-specific** TE inference (assignments are random). If you later need locus resolution, do a second, unique-only STAR run.
    

---

##### Why these choices

- **Two-pass + BySJout**: higher confidence in splice junctions improves all gene-centric analyses (exon/intron/genebody) and prevents spurious intron signal arising from mis-spliced alignments.
- **EndToEnd**: avoids soft-clipping that can create misleading partial matches in repetitive DNA—useful when you also care about TE signal and clean intron/exon boundaries.
- **Relative mismatch cap (6%)**: keeps alignments specific without being so strict that you lose real signal. If your RNA is lower quality/longer reads, you can relax to **0.08–0.10**.
- **High multimap retention + single emitted alignment**: preserves TE information during mapping and implements the **random-one** TE strategy at the alignment stage, so downstream counting is simple and reproducible.
- **Reproducibility**: fixed seed guarantees identical random assignments across re-runs.
    

##### Small optional tweaks 

 **Relax mismatch proportion if needed**  
    If you see lower than expected alignment rates (e.g., older libraries, long reads), try:
    
```
--outFilterMismatchNoverLmax 0.08
```

**Junction permissiveness**  
    Your `BySJout` is conservative (good default). If you work with unusual splicing (e.g., immune hyper-editing, many non-canonical SJs), you might drop it—but only if you have a reason.
    

##### Summary: green/yellow/red flags

- ✅ **Green (keep):** two-pass + BySJout; EndToEnd; intron windows; multimap 100; random-one trio; RNG seed.
- 🟨 **Yellow (situational):** `outFilterMismatchNoverLmax 0.06` (relax to 0.08–0.10 if needed); `outSAMunmapped None` (flip to `Within` for QC).
- 🟥 **Red (don’t mix):** Don’t use this BAM for **TE locus-specific** counting or for **fractional TE** counting—its policy is random-one. If you want **fractional**, re-run STAR allowing multiple alignments (`--outSAMmultNmax >1`) and then count TEs with `featureCounts -M --fraction`.
    

With this configuration, you can align **once** and then:
- count **exons/introns/genebodies** with featureCounts on disjoint gene tracks, and
- count **TE families/subfamilies** (integer) with your grouped SAF (no `--fraction`).


# Vignette — TE (family/subfamily) + gene quantification with featureCounts

This guide shows how we:
1. convert the TE GTF into a **grouped SAF** so featureCounts outputs rows like `HAL1:L1:LINE`,
2. (optionally) remove TE loci that overlap annotated exons with **bedtools**,
3. run **TE counts** (fractional multi-mappers; fragments; unstranded),
4. run **gene counts** per library strandedness,
5. do quick **QA checks** so we know we did the right thing.

Everything below uses your actual paths.

---

## 0) Inputs

- **Gene GTF**  
```bash
GENE_GTF=/data1/users/antonz/data/DM_summer_2025/13441-DM/results_mm39_TE/genome/gencode.vM37.primary_assembly.annotation.filtered.gtf
```

- **TE GTF** (from TEtranscripts website, mm39 RepeatMasker in GTF form)  
```bash
TE_GTF=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_Ensembl_rmsk_TE.gtf.gz
```

- **BAMs** (per project):
    - 13441-DM (reverse-stranded genes):  
        `/data1/users/antonz/data/DM_summer_2025/13441-DM/results_mm39_TE/star_salmon/bam_fin/*.bam`
    - 13444-DM (unstranded genes):  
        `/data1/users/antonz/data/DM_summer_2025/13444-DM/results_mm39_TE/star_salmon/bam_fin/*.bam`
        

---

## 1) Build a grouped TE SAF (TEtranscripts-like labels)

### 1A. Make a **grouped SAF** (one row per **locus**, but `GeneID` is the **group label**)

> This lets featureCounts **sum all loci** for the same TE group into one **meta-feature** row in the output.  
> 
> We keep class info so labels look like `HAL1:L1:LINE`.

```bash
TE_GTF=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_Ensembl_rmsk_TE.gtf.gz
SAF_ALL=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_all.saf

zcat "$TE_GTF" |
awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     $0 !~ /^#/ {
       match($0,/gene_id "([^"]+)"/,g);      # subfamily / repName
       match($0,/family_id "([^"]+)"/,f);    # family
       match($0,/class_id "([^"]+)"/,c);     # class (LINE/LTR/SINE/DNA/RC/...)
       gid = g[1] ":" f[1] ":" c[1];         # TEtranscripts-like label
       print gid, $1, $4, $5, $7
     }' > "$SAF_ALL"
```

**Checks:**
```bash
head -3 "$SAF_ALL"
# GeneID  Chr  Start  End  Strand
# Lx2B2:L1:LINE    1   8387807  8388657  +
# B3:B2:SINE       1   41942995 41943142 +
wc -l "$SAF_ALL"                                  # rows = loci + 1 header
cut -f1 "$SAF_ALL" | tail -n +2 | sort -u | wc -l # unique TE groups (≈ 1–2k; you saw 1243)
```

#### 1B. (Recommended) Keep **retrotransposons only** (drop Simple_repeat, etc.)

```bash
SAF_RETRO=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_retro.saf

zcat "$TE_GTF" |
awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     $0 !~ /^#/ {
       match($0,/class_id "([^"]+)"/,cl); cls = cl[1]
       if (cls=="LINE" || cls=="SINE" || cls=="LTR" || cls=="RC") {
         match($0,/gene_id "([^"]+)"/,g);
         match($0,/family_id "([^"]+)"/,f);
         gid = g[1] ":" f[1] ":" cls;
         print gid, $1, $4, $5, $7
       }
     }' > "$SAF_RETRO"
```

**Why**: This mirrors common TE-expression analyses and avoids low-complexity/simple repeats and RNA genes that can dominate counts.

#### 1C. (Optional) Remove TE loci overlapping **gene exons** (conservative clean-up)

> We use **bedtools** to subtract exon intervals from TE intervals, then rebuild the SAF.

**Make exon BED (0-based) from the **gene** GTF:**

```bash
GENE_GTF=/data1/users/antonz/data/DM_summer_2025/13441-DM/results_mm39_TE/genome/gencode.vM37.primary_assembly.annotation.filtered.gtf
awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5,".",".",$7}' "$GENE_GTF" > exons.bed
```

#### How do I check which classes are available to make sure what I\`m excluding here?

###### 1) See which `class_id` values exist (with counts)
```bash
TE_GTF=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_Ensembl_rmsk_TE.gtf.gz

# Unique classes, descending by locus count
zcat "$TE_GTF" |
awk '$0!~/^#/{ if (match($0,/class_id "([^"]+)"/,m)) print m[1]; else print "NA" }' |
sort | uniq -c | sort -nr
```

That will show things like `LINE, SINE, LTR, DNA, RC, Simple_repeat, Low_complexity, Satellite, rRNA, tRNA, snRNA, scRNA, srpRNA, Unknown, ...` (exact list depends on the file).

###### 2) (Optional) Drill down by family within a class

```bash
# families per class, with counts
zcat "$TE_GTF" |
awk '$0!~/^#/{ match($0,/family_id "([^"]+)"/,f); match($0,/class_id "([^"]+)"/,c);
      if(f[1]!="" && c[1]!="") print c[1]"\t"f[1]; }' |
sort | uniq -c | sort -nr | head
```

###### 3) Build the whitelist and filter (retrotransposons only)

Once you’ve seen what’s present, set a **regex** for the classes to keep:

```bash
KEEP_CLASSES='^(LINE|SINE|LTR|RC)$'   # edit if you want DNA or others

SAF_RETRO=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_retro.saf
zcat "$TE_GTF" |
awk -v pat="$KEEP_CLASSES" 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     $0!~/^#/{
       match($0,/class_id "([^"]+)"/,cl); cls=cl[1];
       if (cls ~ pat){
         match($0,/gene_id "([^"]+)"/,g);
         match($0,/family_id "([^"]+)"/,f);
         print g[1] ":" f[1] ":" cls, $1, $4, $5, $7
       }
     }' > "$SAF_RETRO"
```

###### 4) Sanity-check the filtered SAF really contains only your classes

```bash
# classes present in the SAF (3rd token in GeneID)
cut -f1 "$SAF_RETRO" | tail -n +2 | awk -F: '{print $3}' | sort -u
```

If that prints only `LINE, LTR, RC, SINE` (or whatever you chose), you’re good.

###### 5) (Optional) See examples for a given class

```bash
zcat "$TE_GTF" | grep 'class_id "LTR"' | head
```

###### 6) If you later want to exclude more (e.g., “DNA”, “Unknown”)

Just extend the regex, e.g.:

```bash
KEEP_CLASSES='^(LINE|SINE|LTR)$'      # drop RC (Helitrons)
```

That’s all you need to **inspect** and **control** class filtering with confidence.


#### **Make TE BED from SAF (0-based) and normalize contig names if needed:**

```bash
IN_SAF="$SAF_ALL"  # or $SAF_RETRO
awk 'BEGIN{OFS="\t"} NR>1{print $2,$3-1,$4,$1,".",$5}' "$IN_SAF" > te_grouped.bed


# Make sure we normalized contig names the same way in both inputs before running bedtools
# (either both with 'chr' or both without)
grep -m1 '^chr' exons.bed       || echo "exons: no chr prefix"
grep -m1 '^chr' te_grouped.bed  || echo "te:    no chr prefix"

# If exons.bed uses chr* but TE does not (or vice versa), normalize:
# Strip 'chr' from exons to match TE (Ensembl mm39 often has no 'chr'):
# (If your exons DO have 'chr' and TE does NOT, use this)
sed -E 's/^chr//' exons.bed > exons.tmp && mv exons.tmp exons.bed
```

**Sort both files (bedtools is happier & faster on sorted inputs):**

```bash
LC_ALL=C sort -k1,1 -k2,2n exons.bed        -o exons.bed
LC_ALL=C sort -k1,1 -k2,2n te_grouped.bed   -o te_grouped.bed
```

**Measure overlaps, subtract, verify:**

```bash
bedtools intersect -a te_grouped.bed -b exons.bed -u | wc -l   # how many TE loci overlap exons?
bedtools subtract  -a te_grouped.bed -b exons.bed > te_grouped_noExon.bed
bedtools intersect -a te_grouped_noExon.bed -b exons.bed -u | wc -l   # should be 0
```

**Rebuild a SAF from the exon-removed BED:**

```bash
#SAF_ALL=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_all.saf
SAF_NOEXON=/data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_GROUPED_all_noExon.saf
awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     {print $4,$1,$2+1,$3,$6}' te_grouped_noExon.bed > "$SAF_NOEXON"

# Check counts again
wc -l "$SAF_NOEXON"
cut -f1 "$SAF_NOEXON" | tail -n +2 | sort -u | wc -l   # groups should remain ~ the same (≈1243)
```

> **How bedtools works (brief):**
> - `intersect -u` reports intervals from **A** that have **any** overlap with **B** (we used it to _count_ overlaps).
> - `subtract` removes the parts of **A** that overlap **B** (we used it to drop exon-overlapping TE loci segments).
> - We used 0-based BED (hence `start-1` when converting from GTF/SAF to BED, and `+1` when converting back).