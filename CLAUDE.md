# multispectral_drone_svh

Code supporting Gemmell et al., "Applying the spectral variability hypothesis to arid shrublands, using multispectral drone imagery."

The analysis tests whether spectral heterogeneity from drone-borne multispectral imagery (5 bands: blue, green, red, red-edge, NIR) predicts taxonomic plant diversity across four AusPlot sites in arid NSW (NSABHC0009–0012), each surveyed as a 5×5 grid of 20 m subplots.

## Repo layout

- `_targets.R` — **canonical pipeline driver** (Phase 4). `Rscript -e 'targets::tar_make()'` runs the whole analysis. See "Phase 4 notes" below for design context.
- `continuous_metrics_analysis.R` — **legacy driver**, retained as a fallback during the Phase 4 transition. Computes CV, SV, 5D CHV using rarefaction (n=999); fits per-taxonomic-metric × per-spectral-metric mixed models. Its model-fitting logic now lives in `funx.R::fit_spectral_biodiversity_model()` and is called from `_targets.R`. Plan to retire after Phase 5.3 lands `report.Rmd`.
- `spectral_species_analysis.R` — **legacy driver**, same status. Per-seed RF + k-means clustering logic is now in `funx.R::spectral_species_one_seed()` and called from `_targets.R` via a `tar_map` over 20 seeds.
- `figures-for-publication.R` — **legacy figure generator**, not yet migrated. Reads `data_out/` outputs and produces Figure 5 (`masked_24_plot.png`) and Figure 6 (`cv_band_combinations_plot.png`). Re-fits the 12 mixed models inline (redundant with the cached fits in the targets DAG). Phase 5.3 (`report.Rmd`) should absorb these as `tar_render()` chunks reusing cached fits via `predict()`.
- `funx.R` — analysis-specific helpers and saltbush-delegating wrappers. The three "main" reusable functions (`calculate_field_diversity`, `extract_pixel_values`, `calculate_spectral_metrics`) are now thin wrappers around `saltbush::` equivalents — the wrappers handle AusPlots-specific subplot binning and the project's `subplot_id` / `CHV_nopca` column conventions. What's still genuinely local: `calculate_coefficient_of_variance` (band-subset CV; no saltbush equivalent), `calculate_cv` (called only by that), `bin_survey_subplots` (transect → 5×5 grid), `download_zenodo_rasters` (Zenodo fetcher), and the preprocessing-only helpers `create_masked_raster` / `find_optimum_thresholds`.
- `data/ausplots_march_24.csv` — field survey hits.
- `data/fishnets/NSABHC00{09..12}_fishnet.shp` — the 5×5 subplot grids per site.
- `data/raster_images/` — gitignored; populated on first run by `download_zenodo_rasters()` from Zenodo 10.5281/zenodo.17089161 (~4.4 GB, 4 files: 0009 1.5 GB, 0010 554 MB, 0011 957 MB, 0012 1.5 GB).
- `_targets/` — gitignored; the canonical pipeline cache (Phase 4). Content-addressed; reproduced exactly from the lockfile on a fresh checkout via `tar_make()`. Inspect via `targets::tar_load(<name>)` / `targets::tar_meta()`.
- `data_out/` — gitignored; **legacy outputs** from the top-level analysis scripts. Will become redundant once Phase 5.3 (`report.Rmd`) replaces `figures-for-publication.R` as the figure source and the legacy scripts are retired.
  - `spectral_taxonomic_diversity.rds`, `cv_values.rds` — extracted pixel-derived metrics joined to taxonomic diversity. Equivalent of `tar_read(spectral_taxonomic_diversity)`.
  - `spectral_biodiversity_model_results.rds`, `cv_biodiversity_model_results.rds` — 12-model and CV-band-combo loop output. Equivalent of the eponymous `tar_read()` targets.
  - `model_checks/` — diagnostic per-model `.rds` (mixed-effect specification only). Use these to inspect random-effect variances before deciding whether a refit is warranted. No targets equivalent yet — produced only by the legacy script.
  - `model_fits/` — final fits one `.rds` per taxonomic × spectral metric, suffixed `_mixed` or `_fixed` depending on whether the singular-fit refit kicked in, plus `model_summaries.csv` and `top_significant_models.csv`. Now superseded by the `model_<tax>_<spec>` and `summary_<tax>_<spec>` targets.
- `reports/` — tracked: `interim_progress.Rmd` + its rendered `.html`, plus `main_results_summary.txt`. The directory itself is **not** in `.gitignore` (only `*.log` is); only transient figure drafts (`workflow_*.png`, `workflow_*.mmd`) are intentionally untracked. Read `interim_progress.Rmd` for the current running state of the analysis. To be revisited in Phase 5.8 once `report.Rmd` (Phase 5.3) decides the fate of the interim report.
- `maps_graphs/` — gitignored; output directory for `figures-for-publication.R` (`masked_24_plot.png`, `cv_band_combinations_plot.png`).

## Required R packages

The canonical pinned dependency set lives in `renv.lock` (committed in Phase 5.1). Use `renv::restore()` to install. The README has a grouped human-readable summary. The headline packages, for orientation: `sf` + `terra` (spatial), `glmmTMB` + `performance` + `vegan` + `geometry` (modelling, with `geometry` doing the 5D CHV), `randomForest` + `cluster` (spectral species), `targets` + `tarchetypes` (pipeline), `ggraph` + `tidygraph` + `ggnewscale` + `ggh4x` (figures), `pROC` (mask thresholds, used in preprocessing helpers only), `saltbush` from GitHub at commit `aece6a1`.

## Conventions

- Per-subplot IDs are `"row_col"` (e.g. `"3_2"`); after joining across sites the analysis prefixes them with a site letter (`E`/`G`/`S`/`C`) to keep them unique.
- Bands must be stacked in wavelength order (blue, green, red, red-edge, NIR) — several functions assume this implicitly.
- Site names are the AusPlot IDs (`NSABHC0009`..`0012`) and are extracted from filenames via `str_extract(basename(x), "^[^_]+")`.
- **Rarefaction seed:** `RAREFACTION_SEED <- 42` is set at the top of **both** `_targets.R` and the legacy `continuous_metrics_analysis.R`, threaded into `calculate_spectral_metrics()` and `calculate_coefficient_of_variance()` so CV/CHV are bit-exact reproducible across runs. The spectral-species analysis loops over seeds 1–20 internally and does not need a single global seed.
- **Singular-fit refit:** the three `pielou_evenness` mixed models (CV, SV, log.CHV) consistently produce site random-intercept std.dev ≈ 1e-6 (`performance::check_singularity() == TRUE`). The convention is to refit those models as fixed-effect linear models (no `(1 | site)`). This logic now lives in `funx.R::fit_spectral_biodiversity_model()` (auto-detects singularity and returns an `lm` instead of a `glmmTMB` fit); the model-results tables tag the row with `model_kind = "fixed"` vs `"mixed"`. In the legacy script's `data_out/` outputs the same models are persisted with suffix `_fixed` next to their `_mixed` diagnostic counterparts.

## Current analysis state

For the running state of the analysis (which scripts have been run end-to-end, current model results, decisions awaiting co-author review) see `reports/interim_progress.Rmd`. That report is the canonical place for transient state. CLAUDE.md only documents durable structure, conventions, and the phase plan.

---

# Phase plan

Six phases total. Phases 0–4 are complete; Phase 5 is in progress. The two function refactors (Phases 1 and 3) are bracketed by reproducibility work — seed wiring before swapping anything, `renv` lock after the structure stabilises.

| Phase | Status | Summary |
|---|---|---|
| 0 — Test scaffolding | ✓ done | `testthat` + Tier 1/2 fixtures + `seed = NULL` on rarefaction functions |
| 1 — `raster` → `terra` | ✓ done | `create_masked_raster()` migrated; no `raster::` calls remain |
| 2 — Reproducibility groundwork | ✓ done | seed wired into analysis, LICENSE (MIT + CC-BY-4.0), CITATION.cff, README setup |
| 3 — `funx.R` → `saltbush` | ✓ done | three wrappers delegate to saltbush; band-subset CV stays local |
| 4 — `targets` pipeline | ✓ done | full DAG; first cold end-to-end run 2026-05-22 → 2026-05-24, ~44 h wall-clock, 84 targets |
| 5 — Publication prep | in progress | 5.1 done; see below |

Each phase ends with `testthat::test_dir("tests")` passing.

---

# Testing strategy

Phase 0 added `testthat` so subsequent refactors are verified mechanically by snapshot diff.

## Tiers

- **Tier 1 — every change (~1 s).** Taxonomic-diversity snapshot. Uses only `data/ausplots_march_24.csv`, no rasters. Exercises `calculate_field_diversity()`, which both analysis scripts call. Fixture: `tests/fixtures/taxonomic_diversity_baseline.rds`.
- **Tier 2 — per phase, manual (~minutes).** One-site `extract_pixel_values` + `calculate_spectral_metrics` snapshots on NSABHC0010 (smallest raster at 554 MB). Uses `n = 99` for rarefaction in tests, `n = 999` in the real analysis. Fixture generated once, re-checked after each function swap. Lives in `tests/test-spectral-metrics.R` and is **gated by an env var** — skipped by default; opt in with `RUN_TIER2=true Rscript -e 'testthat::test_dir("tests")'`. The test also skips automatically if the Zenodo raster isn't downloaded or if `geometry` isn't installed, so adding it doesn't break Tier 1 runs on a fresh checkout.
- **Tier 3 — release / merge to main (~hours).** Full end-to-end on all four sites. Acceptance check, not a per-commit gate.

## Tooling

- `testthat` for structure; `waldo::compare()` for readable diffs.
- Snapshots checked into the repo as `.rds` files in `tests/fixtures/` — small (summary tibbles, not rasters).
- Project is not an R package (no `DESCRIPTION`), so tests live in flat `tests/` (no `tests/testthat/` subdir). Run with `Rscript -e 'testthat::test_dir("tests")'`.

## When a snapshot diff appears

A refactor changing outputs isn't automatically wrong — `saltbush::calculate_spectral_metrics` may legitimately differ from the local version in defaults. The workflow is:
1. Run the test, get a diff.
2. Inspect it — is the change explainable (different default, different sampling order)?
3. If yes: either match the package by adjusting arguments, or accept the new behaviour and regenerate the fixture, recording why in the commit message.
4. If no: investigate before merging.

---

# Phase 4 notes: the `{targets}` pipeline as built

Phase 4 (commits `f0fecd1` → `e0bdc60`) replaced the two top-level scripts with a `targets` workflow so expensive steps (n=999 rarefaction, 20-seed RF clustering) are cached and only re-run when inputs change. The legacy scripts are retained as fallbacks until Phase 5.3 lands `report.Rmd`. What follows is design rationale and operational notes for maintenance — not a forward plan.

## Why we adopted targets

- **Expensive, deterministic steps.** Rarefaction in `calculate_cv`/`calculate_chv_nopca` and the 20-seed clustering both take a long time and produce reproducible outputs given a seed — exactly what `targets` caches well.
- **Branching is natural.** The analysis is a cross-product over `site × taxonomic_metric × spectral_metric` (and `× band_combination`, `× seed`). `tar_map()` and `tar_map_rep()` express this directly.
- **Large file inputs.** `tar_target(format = "file")` on the Zenodo download lets `targets` skip re-downloading if files are unchanged.
- **Two scripts share inputs.** `funx.R::calculate_field_diversity` is called from both; in a pipeline it's computed once.

## Phase 4a refactors (done)

For history: Phase 4a extracted pure functions out of the side-effectful analysis scripts before any `_targets.R` existed. Concretely: model-fitting loops became `fit_spectral_biodiversity_model()` and `fit_cv_band_model()`; the per-seed RF clustering became `spectral_species_one_seed()`; the `write.csv` side effect in the seed loop was removed; the site-prefix logic (`E`/`G`/`S`/`C`) was consolidated into `add_site_prefix()`. The figure helpers (item 6 in the original plan: `make_figure_5`, `make_figure_6`) were **not** extracted — they're still in `figures-for-publication.R` pending Phase 5.3.

## Target graph

The live DAG (84 targets after Phase 4d) is best inspected with:

- `targets::tar_visnetwork()` — interactive HTML widget, colour-coded by status.
- `targets::tar_mermaid()` — Mermaid string suitable for embedding.
- A hand-simplified Mermaid version (collapsing per-site and per-pair fan-out into single annotated boxes) lives in the README's *Workflow* section — that's the canonical "what readers see" diagram.

## Sub-phases (as built)

4a → 4d landed in order: 4a (function extraction; no `_targets.R` yet), 4b (minimal pipeline with file inputs + taxonomic diversity + per-site pixel extraction), 4c (CV/SV/5D-CHV + 12 continuous-metric models), 4d (20-seed spectral-species `tar_map` + 3 spectral-species models). Each sub-phase ended with `tar_make()` re-running cleanly and `testthat::test_dir("tests")` passing. The rendered manuscript-figure report is implemented in Phase 5.3 as a `tar_render()` target within this pipeline.

**Runtime capture.** `tar_meta(fields = "seconds")` records per-target wall-clock time automatically — no instrumentation needed. The README's runtime table is sourced from a snapshot of this. To regenerate after a future cold run: `targets::tar_meta(fields = "seconds")` and pick out the dominant rows.

**Observed runtimes from the full Phase 4 cold run (2026-05-22 → 2026-05-24).** All numbers from `tar_meta(fields = "seconds")` on Apple M-series silicon. Total wall-clock: **~44 h** across 84 targets.
- Per-site `pixel_values_<site>` extraction: 16–27 s each.
- Combined `pixel_values` + `min_points`: ~12 s.
- `cv_subset_red.edge_nir` (CV only, 2 bands): ~16 min.
- `cv_subset_green_red.edge_nir` (CV only, 3 bands): ~17 min.
- `cv_subset_green_red_red.edge_nir` (CV only, 3 bands): ~19 min.
- `spectral_metrics` (CV + SV + **5D CHV**, all bands, n=999): **~3 h 5 min** (11,075 s on this cold run; an earlier mid-Phase-4c estimate was ≥50 min — see "When to start worrying" below). cv_subsets and `spectral_metrics` are not in the same cost regime: cv_subsets do CV only, while `spectral_metrics` adds the 5D convex hull, and `geometry::convhulln` cost grows roughly as O(n^(d/2)) with dimension. Comparing them as "how much extra do bands cost" is misleading; the dominant factor is the CHV dimensionality, not the band count.
- 20 × `spectral_species_seed_<s>`: 112–138 min each (mean ~120 min), 40.1 h cumulative serial.
- 9 of the 20 seeds emit a benign `randomForest` "did not converge in 10 iterations" warning — same warning the legacy `spectral_species_analysis.R` produces. All seeds completed.
- Downstream of `spectral_metrics` / `mean_spectral_species`: the `spectral_taxonomic_diversity` join, both `tar_combine`s, and the 24 glmmTMB/lm fits + summaries are all fast (seconds each on ~100-row data).

**Why `spectral_metrics` can't be parallelised away.** It's tempting to look at this hour-long monolithic target and split it per site for parallelism. Don't. Design decision #2 in `_targets.R` fixes `RAREFACTION_SEED = 42` and pins `spectral_metrics` to a *single* call so the rarefaction RNG sequence is reproducible — splitting per site would change the published numbers. The cost is the price of bit-exact reproducibility, and it's the right tradeoff. This target is also the strongest motivating example for the "Why `targets`" section in the Phase 5 README: hour-long, deterministic given a seed, reused by every downstream model and figure — exactly what content-addressed caching is for.

**When to start worrying.** The first full cold-run measurement (184 min) is well above the ≥90 min "worry" threshold from the Phase 4c mid-run estimate. Two possibilities, not yet disambiguated: (a) the original ≥50 min figure was timed on a partial dataset and 184 min is the real cost of all four sites at n=999, or (b) one site's 5D hull hits a Qhull degenerate-point edge case (which can slow `convhulln` dramatically). Phase 5 housekeeping (5.8) should bisect by site via `tar_workspace()` to nail this down before publication. In the meantime, anything over ~4 h on a future run definitely warrants investigation.

## Working conventions during the transition

- Don't add new top-level scripts with side effects. Logic goes in `funx.R` (or a future `R/` directory) as named functions.
- Keep functions pure where possible: take data + parameters, return data. No `write.csv` or `dir.create` inside metric/model functions.
- Prefer adding a new function over modifying a working one — keeps existing scripts runnable as a fallback during the refactor.
- The Zenodo downloader is the canonical raster source. Don't reintroduce the legacy `data_out/combined_rasters/masked/2024/` path.

---

# Plan: publication prep (Phase 5)

Final-mile items before the code accompanies the manuscript. Most depend on the targets pipeline being in place so "run the analysis" is a one-liner.

## Status and dependency graph

| # | Item | Status |
|---|---|---|
| 5.1 | Lock dependencies with `renv` | ✓ done (commit `d741860`) |
| 5.2 | Rewrite README | ✓ done (commit `1ace498`) |
| 5.3 | Rendered `report.Rmd` via `tar_render()` | pending |
| 5.4 | Dockerfile | pending |
| 5.5 | CI for Tier 1 tests | ✓ done (commits `02f650c`, `3530e10`) |
| 5.6 | Output checksums as release artifacts | pending |
| 5.7 | Code DOI via Zenodo | pending |
| 5.8 | Repo housekeeping pass | partially in flight (drift fixes have been folded into the relevant doc-update commits as they're spotted) |

Hard ordering (from item-1 dependencies): 5.1 must precede 5.2, 5.4, 5.5 (those three reference `renv::restore()` or the locked environment).

Soft ordering (release shape): 5.7 should be last among the artifact items because the v1.0.0 release that mints the code DOI should be a complete, documented, reproducible artifact. 5.8 is the final reconciliation sweep.

```
 5.1 (renv) ─┬─→ 5.2 (README)     ─┐
            ├─→ 5.4 (Dockerfile)  ─┤
            └─→ 5.5 (CI)          ─┤
                                   ├─→ 5.7 (DOI) ──→ 5.8 (housekeeping)
 5.3 (report.Rmd) ─────────────────┤
 5.6 (checksums) ──────────────────┘
```

5.2, 5.3, 5.4, 5.5, 5.6 can proceed in any order or in parallel once 5.1 is in.

## Items

1. **Lock dependencies with `renv`.** ✓ done in commit `d741860`. `renv::init()` captured R 4.5.3 and the full pipeline dependency set, with `saltbush` pinned to commit `aece6a18`. Tier 1 tests pass against the locked environment.
2. **Rewrite README.** ✓ done in commit `1ace498` (+ saltbush hex callout added in commit `bd9f577`). Replaced the transitional README with a real entry point. Sections:
   - **Prerequisites** — R version, GDAL/PROJ, C++17 for `geometry`.
   - **Install** — `renv::restore()`. List the full pinned dependency set with versions (pulled from `renv.lock` after item 1 lands), grouped by role: spatial (`sf`, `terra`), modelling (`glmmTMB`, `performance`, `vegan`, `geometry`), clustering (`randomForest`, `cluster`), pipeline (`targets`, `tarchetypes`, optionally `crew`), figures (`ggnewscale`, `ggh4x`, `ggraph`, `tidygraph`), project-specific (`saltbush` at the pinned commit).
   - **Quick start** — `targets::tar_make()`.
   - **Expected outputs and runtime** — table of per-target wall-clock times sourced from `tar_meta()` (captured at the end of Phase 4; see "Runtime capture" above). Include total end-to-end runtime, the expensive sub-totals (rarefaction at n=999, 20-seed RF clustering), and the ~4.4 GB Zenodo download as a one-time cost.
   - **Workflow diagram** — reviewer asked for this explicitly and it ships in **two places** with slightly different content (not just rendering):
     - **README** — `targets::tar_mermaid()` output embedded inline as Mermaid; renders natively on GitHub, no local toolchain. Audience is users reproducing the code: the preprocessing block (per-band stacking + ROC-thresholded masking) is shown inside a **dashed boundary** with "done once; outputs archived on Zenodo" so readers know they download the masked rasters and skip preprocessing entirely. Mention `targets::tar_visnetwork()` as the interactive alternative.
     - **Manuscript** — rendered with `ggraph` from `targets::tar_network()`'s `vertices` / `edges` data frames; shares fonts, theme, and palette with Figures 5 and 6. Audience is paper readers: preprocessing is shown as part of the methods narrative (**no dashed boundary**), since the reader is following the analysis flow, not deciding what to skip. Add `make_figure_workflow(network)` to the figure helpers alongside `make_figure_5` / `make_figure_6`, and wire it in as its own `tar_target` (returning a ggplot, persisted to PDF via `format = "file"` exactly like the other figures).
     **Simplification rules applied** (both versions): collapse per-site fan-out, per-pair model branches, and per-band-subset CV into single annotated boxes; drop `format = "file"` file targets and `tar_combine` plumbing; keep rarefaction parameters (`n=999, seed=42`) and 20-seed averaging visible as box annotations since these are scientific choices reviewers will ask about.
     **Content** (both versions, top to bottom): preprocessing inputs (drone tiles per site + visually classified training pixels on veg/non-veg **and** shadow/non-shadow axes) → Pix4D orthomosaic + illumination correction + 5-band stacking → ROC-optimised per-site thresholds (**NIR + NDVI only**; a Red filter exists in `saltbush` but was not applied to this dataset) → masked rasters → analysis inputs (fishnets, AusPlots survey) → three parallel tracks (continuous spectral metrics CV/SV/5D-CHV; spectral species via RF + k-means k=40 averaged across 20 seeds; taxonomic diversity S/exp(H')/1/D/J') → mixed models with singular-fit refit convention → Figures 5 + 6.
     Earlier drafts (post-4b skeleton, post-4c with the metric cross-product) are usable for PR review but the canonical version is the post-4d one that includes the 20-seed spectral-species fan-out.
   - **Why `targets`** — short rationale for readers who'll wonder why a paper repo isn't just two scripts. Expand the "Why targets fits" notes from the Phase 4 plan into reader-facing prose, grounded in the runtime numbers from the section above: (a) the expensive deterministic steps (rarefaction, 20-seed clustering) are cached and skipped on re-runs; (b) the `site × taxonomic_metric × spectral_metric × band_combination × seed` cross-product is expressed declaratively via `tar_map`/`tar_map_rep` rather than nested loops; (c) `tar_visnetwork()` gives reviewers a live DAG of the analysis; (d) input changes (a new raster, a corrected survey row) invalidate only the affected branches, so partial reruns are cheap. The argument is: for a multi-hour, multi-input, cross-product analysis like this one, a Make-style DAG with content-addressed caching is the smallest tool that actually delivers "re-run only what changed" — which is what reproducibility means in practice once the analysis is larger than a single script.
   - **Data and code citations** — Zenodo data DOI (10.5281/zenodo.17089161) and the code DOI from item 7.
   - **License.**
3. **Rendered report.** `report.Rmd` that loads pipeline outputs via `tar_read()` and produces the manuscript's tables and figures. Wired into the pipeline as a `tar_render()` target so it rebuilds whenever upstream targets change. This is the artifact that ties the code in the repo to the claims in the paper.
4. **Dockerfile.** Base off `rocker/geospatial:<R version>` (ships GDAL/PROJ pre-installed); call `renv::restore()` at build time; entrypoint runs `targets::tar_make()`. Lets reviewers and readers reproduce the analysis end-to-end without configuring a local R environment.
5. **CI for Tier 1 tests.** ✓ done in commits `02f650c` (workflow) + `3530e10` (sysdeps fix). `.github/workflows/tests.yml` runs on push and PR against `main`, plus `workflow_dispatch`. ubuntu-latest, pinned R 4.5.3, `renv::restore()`, then `testthat::test_dir("tests")`. Tier 2/3 deliberately not run in CI — the rasters are 4.6 GB and GitHub-hosted runners won't tolerate it. First cold run was ~4 min; subsequent runs are faster from cache.
6. **Output checksums as release artifacts.** Ship the Tier 1 baseline RDS (and optionally one Tier 2 fixture) alongside the GitHub release so readers can `digest::digest()` their replication output and confirm an exact match.
7. **Code DOI via Zenodo–GitHub integration.** Enable the integration; cut a `v1.0.0` release; Zenodo mints a DOI for the code, separate from the data DOI (10.5281/zenodo.17089161). Cite both in the manuscript. Update CITATION.cff with the code DOI once it exists.
8. **Repo housekeeping pass.** Partially in flight — drift fixes are being folded into the relevant doc-update commits as they're spotted (e.g. the `reports/` gitignored claim was fixed in `1ace498`; the workflow's NDVI/NIR/Red claim corrected to NIR+NDVI only in the same commit; the Phase-4-as-plan-vs-built reframing in the staleness sweep). Outstanding for a dedicated sweep before the v1.0.0 release: decide the fate of `reports/interim_progress.Rmd` once `report.Rmd` (Phase 5.3) lands; check whether the rendered HTML stays tracked as a published artifact; audit `interim_progress.Rmd` itself for staleness; one more grep pass over CLAUDE.md and the README for any remaining mismatches.

Order matters: `renv` must land before items 2 (README references `renv::restore()`), 4 (Dockerfile uses `renv::restore()`), and 5 (CI runs against the locked environment).

Tier 1 + Tier 2 tests should pass against the locked environment after each item lands.

---

# Potential future improvements (post-publication)

Not on the Phase 5 critical path; record here so they're easy to pick up later.

## Make Zenodo data-version bumps a one-line change

`download_zenodo_rasters()` in `funx.R` currently hardcodes four `https://zenodo.org/records/17089161/files/<name>.tif?download=1` URLs and uses **file-size within 5%** as the "already downloaded?" check. Consequences if the underlying dataset is ever republished as a new Zenodo version:

- URLs in `funx.R` need editing in four places (would only need one if the record ID were factored out).
- The `size = ...` expected sizes need updating.
- Anyone with the *old* rasters already cached in `data/raster_images/` will silently skip the re-download — `targets` won't detect the version change because the input is a regular value target, not a `format = "file"` target.
- The data DOI also needs updating in `README.md` (×2), `CITATION.cff` (×2), and `CLAUDE.md` (×3).

Suggested refactor (~30 lines):

1. Centralise the Zenodo record ID and per-file SHA-256s at the top of `funx.R` (or in `_targets.R`):
   ```r
   ZENODO_RECORD <- "17089161"
   ZENODO_RASTERS <- list(
     list(name = "NSABHC0009_masked.tif", sha256 = "...", size = 1.6e9),
     ...
   )
   ```
2. In `download_zenodo_rasters()`, replace the size-check with `tools::sha256sum(dest) == f$sha256`. Skip if the local hash matches; download (and re-verify) otherwise. Hard-fail on hash mismatch after download.
3. Optionally wrap each per-site `raster_<site>` target in a `format = "file"` check that depends on `ZENODO_RECORD` so the DAG invalidates on bump.

Net result: bumping the data version becomes "edit one constant + paste new SHA-256s, push, let `tar_make()` recompute."

Unlikely to be needed for this paper — the dataset is what it is — but cheap to implement and a worthwhile reproducibility upgrade if a follow-up paper uses the same downloader pattern.

## Parallelise the spectral-species seeds

The 20 spectral-species seeds run serially in the current pipeline (~40 h cumulative). They're embarrassingly parallel: each seed is an independent `tar_target`. Wiring in [`{crew}`](https://wlandau.github.io/crew/) (a `tar_make()`-aware backend) would knock the wall-clock to ~5 h on an 8-core machine without changing any of the science. Two lines in `_targets.R`:

```r
library(crew)
tar_option_set(controller = crew_controller_local(workers = 8))
```

Held off during Phase 4d because the weekend serial run was good enough and adding `{crew}` to the lockfile right before a critical run was higher risk than benefit. Worth doing if anyone needs to re-run the full pipeline regularly.

## Investigate the `spectral_metrics` 3 h 5 min runtime

See "When to start worrying" above. The first full cold-run measurement (184 min) is well above the original ≥90 min estimate from the partial Phase 4c run. Could be the real cost of all four sites at n=999, or a Qhull degenerate-point edge case in one site's 5D hull. `tar_workspace(spectral_metrics)` + per-site bisection would resolve it. Not blocking publication; would shave hours off any future cold run if a site turns out to be pathological.
