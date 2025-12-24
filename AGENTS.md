# Repository Guidelines

## Project Structure & Module Organization
- Core R functions live in `R/` grouped by phase: `compute_*` (Phase 1 compute/cache), `format_*` (Phase 2 master tables), `plot_*` (Phase 3 viz), plus factories (`create_*_dge.R`) and utils.
- CLI helper scripts and TE-specific workflows are under `scripts/`; example end-to-end pipelines are in `examples/` (combined vs separate modes).
- Configuration templates sit in `inst/config/` (see `te_config_template.yaml`), and narrative docs are in `docs/`.
- Integration and validation checks are in `tests/`. Keep new assets (figures, CSVs) out of the repo unless essential.

## Build, Test, and Development Commands
- Source toolkit in R: `source_toolkit("TE-RNAseq-toolkit")` or `source(file.path("R", "<file>.R"))` for iterating locally.
- Run example pipelines: `Rscript examples/1.5.te_processing_combined.R` (genes+TEs) or `Rscript examples/1.5.te_processing_separate.R` (TE-only).
- Run tests (fast, synthetic): `Rscript tests/test_validation_synthetic.R`.
- Run integration test (requires `/workspaces/AdaW_eWAT_WL_2025/00_data/processed`): `Rscript tests/test_integration_real_data.R`.

## Coding Style & Naming Conventions
- Use lower_snake_case for files, functions, and objects; keep verbs up front (`compute_te_*`, `format_te_*`).
- Prefer 4-space indentation, explicit `return()` only when clarifying branches, and pipe-friendly tidyverse patterns.
- Document functions with roxygen2 blocks; keep parameter names aligned with existing factories and `compute_*` helpers.
- Avoid hard-coding paths—pass them as arguments or read from config YAML.

## Testing Guidelines
- Name tests `test_*.R` and make them runnable via `Rscript tests/<file>.R` without extra flags.
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
- When sharing pipelines, point to example scripts rather than duplicating logic to avoid drift.
