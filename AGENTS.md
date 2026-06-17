# AGENTS.md

TE-RNAseq-toolkit is a set of scripts and wiki guidances on analyzing transposable elements (TEs) in short-read bulk RNA-seq data.

## Project Structure
- Core R functions in `R/` grouped by phase
- CLI helper scripts and TE-specific workflows are under `scripts/`; example end-to-end pipelines are in `examples/` 
- Configuration in `inst/config/` (see `te_config_template.yaml`)
- narrative docs are in `docs/`
- Tests in `tests/, runnable via `Rscript tests/<file>.R`

## Architecture

Most scripts in repo deal with convenience parsing, analysis and viz of already annotated count matrices

### Three-Phase Design
1. **Phase 1 (Compute)**: `R/compute_*.R` - Expensive calculations cached via `load_or_compute()`
2. **Phase 2 (Tables)**: `R/format_*.R` - Export to standardized "Master Tables" (CSV)
3. **Phase 3 (Viz)**: `R/plot_*.R` - Pure visualization, reads from checkpoints/tables only

For details on preprocessing reqs see below, when necessary — read the in-depth guide of docs/METHODOLOGY.md

A different repo, SciAgent-toolkit, packages convenience post-nf-core env-locked TE+gene featureCounts-based counting
(based on `scripts/runFeatureCounts_TE_and_genes.sh` two-pass driver + BAM recipe + QC gate) in a **`te-gene-featurecounts`** skill, 
with its own pinned container **`te-fc:2.0.2`** (featureCounts v2.0.2). Use that image to run counting on correctly aligned data. 
See alignment pre-reqs below.

### Core Function Categories

| Category | Files | Purpose |
|----------|-------|---------|
| Factory | `create_combined_dge.R`, `create_separate_dge.R` | Build DGEList objects |
| Utils | `te_utils.R`, `validate_te_input.R` | Parsing, annotation, validation |
| Compute | `compute_te_de.R`, `compute_te_enrichment.R`, `compute_te_aggregates.R`, `compute_te_proportions.R`, `compute_te_correlations.R` | Analysis |
| Format | `format_te_de.R`, `format_te_enrichment.R`, `format_te_aggregates.R` | Master table export |
| Plot | `plot_te_*.R`, `plot_utils.R` | Visualization |

All analysis scripts must source the toolkit functions manually (until this is a proper package):
```r
TOOLKIT_PATH <- "path/to/TE-RNAseq-toolkit/R"
source(file.path(TOOLKIT_PATH, "te_utils.R"))
source(file.path(TOOLKIT_PATH, "create_combined_dge.R"))
# ... source other modules as needed
```

### Design Decisions

- **TE ID Format**: Expects TEtranscripts-style `Subfamily:Family:Class` (e.g., `L1Md_A:L1:LINE`), `te_utils.R` functions for parsing (`parse_te_id`)
- **Stats**: Prefer limma-voom (handles both integer and fractional counts)
- **Enrichment**: Uses `limma::geneSetTest` for TE set testing

## Preprocessing requirements
Before using this toolkit, ensure:
1. **STAR alignment** 
For convenience, most personal projects used `--outFilterMultimapNmax 100 --outSAMmultNmax 1 --outMultimapperOrder Random` for Random-One strategy
Can potentially support Fractional (weighted counts) strategy. See docs/METHODOLOGY.md for details
2. **TE annotation** has exonic regions removed (bedtools subtract) for Combined Mode
3. **featureCounts** used `-M` for TEs (retain Random-One multimappers); UNSTRANDED (`--stranded no` follows TETtranscripts defaults) or stranded counting TE at the gene strandedness (`-s 2` for reverse dUTP, follows modern pipelines) 
see docs/METHODOLOGY.md for caveats, controversies and in-depth review

## Coding Style & Naming Conventions
- Use lower_snake_case for files, functions, and objects; keep verbs up front (`compute_te_*`, `format_te_*`).
- Prefer 4-space indentation, explicit `return()` only when clarifying branches.
- Document functions with docstring blocks; 
- keep parameter names aligned with existing factories and `compute_*` helpers.
- Avoid hard-coding paths—pass them as arguments or read from config YAML.

## Common Workflow
Copy `inst/config/te_config_template.yaml` to your project and customize analysis parameters.

## Key Constants (te_utils.R)
```r
TE_CYTOPLASMIC_CLASSES <- c("LINE", "SINE", "LTR", "Retroposon")  # Retrotransposons
TE_NUCLEAR_CLASSES <- c("DNA", "RC")                               # DNA transposons
```

```r
# 1. Create combined DGEList
dge <- create_combined_dge(gene_counts, te_counts, samples)

# 2. Set up design and contrasts
design <- model.matrix(~ 0 + group, data = dge$samples)
contrast_matrix <- makeContrasts(Treatment_vs_Control = Treatment - Control, levels = design)

# 3. Fit model
fit <- compute_te_de(dge, design, contrast_matrix)

# 4. Compute enrichment
enrichment <- compute_te_enrichment(fit, dge, te_level = "family")

# 5. Format and export
write_csv(format_te_enrichment(enrichment), "master_te_enrichment.csv")
```


## Testing Guidelines
- For data-dependent tests, gate with clear checks and log requirements; default to synthetic fixtures where possible.
- Add assertions that validate TE/gene ID collisions, sample alignment, and library-size ratios to mirror existing safety checks.
- Surface human-readable messages (`message()`), and fail loudly on validation errors.

## Commit & Pull Request Guidelines
- Follow existing history: short, present-tense subjects starting with a verb (e.g., `Add`, `Update`, `Fix`); keep scope narrow.
- Reference issues when relevant, and describe data dependencies or expected inputs in the PR body.
- Include what you tested (commands above) and attach key plots/tables for visualization changes.
- Avoid committing generated outputs or large datasets; link to reproducible steps instead.

## Configuration & Data Tips
- Start from `inst/config/te_config_template.yaml`; commit templates, not run-specific copies.
- Keep alignment and counting assumptions consistent with `docs/METHODOLOGY.md` (Random-One STAR params, non-overlapping annotations).
