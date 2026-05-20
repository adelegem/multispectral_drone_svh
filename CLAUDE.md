# multispectral_drone_svh

Code supporting Gemmell et al., "Applying the spectral variability hypothesis to arid shrublands, using multispectral drone imagery."

The analysis tests whether spectral heterogeneity from drone-borne multispectral imagery (5 bands: blue, green, red, red-edge, NIR) predicts taxonomic plant diversity across four AusPlot sites in arid NSW (NSABHC0009–0012), each surveyed as a 5×5 grid of 20 m subplots.

## Repo layout

- `continuous_metrics_analysis.R` — extracts pixel values per subplot, computes coefficient of variation (CV), spectral variance (SV), and convex hull volume (CHV) using rarefaction (n=999); fits per-taxonomic-metric × per-spectral-metric mixed models (`glmmTMB`, site as random effect). Also computes CV across band subsets (red-edge+NIR, etc.).
- `spectral_species_analysis.R` — random-forest + k-means clustering of pixels into "spectral species" (k=40), repeated over 20 seeds; computes spectral richness / Shannon / Simpson per subplot; fits mixed models against field diversity.
- `funx.R` — all the reusable functions: raster I/O, masking, pixel extraction, the three spectral-metric calculators, `calculate_field_diversity` (taxonomic diversity from AusPlots transect hits), and `download_zenodo_rasters` (added recently — auto-fetches the 4 GeoTIFFs from Zenodo).
- `data/ausplots_march_24.csv` — field survey hits.
- `data/fishnets/NSABHC00{09..12}_fishnet.shp` — the 5×5 subplot grids per site.
- `data/raster_images/` — gitignored; populated on first run by `download_zenodo_rasters()` from Zenodo 10.5281/zenodo.17089161 (~4.4 GB, 4 files: 0009 1.5 GB, 0010 554 MB, 0011 957 MB, 0012 1.5 GB).
- `data_out/` — gitignored; intermediate outputs from the analysis scripts:
  - `spectral_taxonomic_diversity.rds`, `cv_values.rds` — extracted pixel-derived metrics joined to taxonomic diversity.
  - `spectral_biodiversity_model_results.rds`, `cv_biodiversity_model_results.rds` — original 12-model and CV-band-combo loop output, written directly by `continuous_metrics_analysis.R`.
  - `model_checks/` — diagnostic per-model `.rds` (mixed-effect specification only) plus `model_checks.log`. Use these to inspect random-effect variances before deciding whether a refit is warranted.
  - `model_fits/` — final fits used for reporting. One `.rds` per taxonomic × spectral metric, suffixed `_mixed` or `_fixed` depending on whether the singular-fit refit (see Conventions) kicked in, plus `model_summaries.csv` and `top_significant_models.csv`.
- `reports/` — gitignored; rendered HTML reports + their sources (currently `interim_progress.Rmd` / `.html`). Read this for the current running state of the analysis.

## Required R packages

`sf`, `terra`, `tidyverse`, `vegan`, `data.table`, `performance`, `glmmTMB`, `geometry` (CHV), `randomForest` + `cluster` (spectral species), `pROC` (mask threshold finder in `funx.R`).

## Conventions

- Per-subplot IDs are `"row_col"` (e.g. `"3_2"`); after joining across sites the analysis prefixes them with a site letter (`E`/`G`/`S`/`C`) to keep them unique.
- Bands must be stacked in wavelength order (blue, green, red, red-edge, NIR) — several functions assume this implicitly.
- Site names are the AusPlot IDs (`NSABHC0009`..`0012`) and are extracted from filenames via `str_extract(basename(x), "^[^_]+")`.
- **Rarefaction seed:** `RAREFACTION_SEED <- 42` at the top of `continuous_metrics_analysis.R` is threaded into `calculate_spectral_metrics()` and `calculate_coefficient_of_variance()` so CV/CHV are bit-exact reproducible across runs. The spectral-species analysis loops over seeds 1–20 internally and does not need a single global seed.
- **Singular-fit refit:** the three `pielou_evenness` mixed models (CV, SV, log.CHV) consistently produce site random-intercept std.dev ≈ 1e-6 (`performance::check_singularity() == TRUE`). The convention is to refit those models as fixed-effect linear models (no `(1 | site)`) and persist them under `data_out/model_fits/` with the suffix `_fixed`. The original mixed-model diagnostic fits stay at `data_out/model_checks/` for comparison. The refit logic is currently ad-hoc and needs to be folded back into `continuous_metrics_analysis.R` — until then, the script's own `*_model_results.rds` outputs contain `NA` for the singular fits' conditional R² and the post-hoc CSVs in `model_fits/` are the canonical results.

## Current analysis state

For the running state of the analysis (which scripts have been run end-to-end, current model results, decisions awaiting co-author review) see `reports/interim_progress.Rmd`. That report is the canonical place for transient state. CLAUDE.md only documents durable structure, conventions, and the phase plan.

---

# Recommended order

Six phases. The two function refactors are bracketed by reproducibility work — pin dependencies before swapping anything, polish for release after the structure stabilises.

0. **Phase 0 — test scaffolding** (no behavioural change). Capture baseline outputs and seed the rarefaction. Without this every later phase risks silent output drift; with it each phase is verified by a one-command snapshot diff.
1. **Phase 1 — `raster` → `terra` migration.** Smallest scope, isolated to one function (`create_masked_raster()`) not called by the analysis scripts. Good first refactor — proves the test harness works on a near-zero-risk change.
2. **Phase 2 — reproducibility groundwork.** Wire the rarefaction seed into the analysis scripts (so end-user runs reproduce), add LICENSE and CITATION.cff, document R version and system deps. Done *before* the function refactors so any later output drift is attributable to the function swap, not to seed addition. Dependency pinning (`renv`) is deliberately deferred to Phase 5 — locking before saltbush / targets means constantly re-snapshotting through the refactors.
3. **Phase 3 — `funx.R` → `saltbush` migration.** Highest behavioural risk but well bounded: one function at a time, verify each swap against the baseline snapshot before moving on. Easiest first: `calculate_field_diversity` → `extract_pixel_values` → `calculate_spectral_metrics`.
4. **Phase 4 — `targets` pipeline.** Purely structural by this point. Phase 4a (factor out model loops, drop CSV side effects) verified by re-running snapshots; Phase 4b (`_targets.R`) is plumbing.
5. **Phase 5 — publication prep.** Lock dependencies with `renv` now that the package set is stable, comprehensive README, rendered report (Rmd → HTML/PDF) that reproduces manuscript figures, Dockerfile based on `rocker/geospatial`, CI for Tier 1 tests, Zenodo–GitHub integration for a code DOI, output checksums as release artifacts.

Each phase ends with `testthat::test_dir("tests")` passing.

---

# Testing strategy

No tests existed on `main` before this work. We're adding them as Phase 0 so subsequent refactors can be verified mechanically.

## Tiers

- **Tier 1 — every change (~1 s).** Taxonomic-diversity snapshot. Uses only `data/ausplots_march_24.csv`, no rasters. Exercises `calculate_field_diversity()`, which both analysis scripts call. Fixture: `tests/fixtures/taxonomic_diversity_baseline.rds`.
- **Tier 2 — per phase, manual (~minutes).** One-site `extract_pixel_values` + `calculate_spectral_metrics` snapshots on NSABHC0010 (smallest raster at 554 MB). Uses `n = 99` for rarefaction in tests, `n = 999` in the real analysis. Fixture generated once, re-checked after each function swap. Lives in `tests/test-spectral-metrics.R` and is **gated by an env var** — skipped by default; opt in with `RUN_TIER2=true Rscript -e 'testthat::test_dir("tests")'`. The test also skips automatically if the Zenodo raster isn't downloaded or if `geometry` isn't installed, so adding it doesn't break Tier 1 runs on a fresh checkout.
- **Tier 3 — release / merge to main (~hours).** Full end-to-end on all four sites. Acceptance check, not a per-commit gate.

## Tooling

- `testthat` for structure; `waldo::compare()` for readable diffs.
- Snapshots checked into the repo as `.rds` files in `tests/fixtures/` — small (summary tibbles, not rasters).
- Project is not an R package (no `DESCRIPTION`), so tests live in flat `tests/` (no `tests/testthat/` subdir). Run with `Rscript -e 'testthat::test_dir("tests")'`.

## Reproducibility prerequisite

`calculate_cv()` and `calculate_chv_nopca()` currently call `sample()` un-seeded, so identical inputs produce non-identical outputs. Phase 0 adds an optional `seed = NULL` argument to both (and to the `calculate_spectral_metrics` / `calculate_coefficient_of_variance` wrappers). Default `NULL` preserves current behaviour; tests pass `seed = 42` or similar to get exact equality. Phase 2 wires that seed into the analysis scripts themselves.

## When a snapshot diff appears

A refactor changing outputs isn't automatically wrong — `saltbush::calculate_spectral_metrics` may legitimately differ from the local version in defaults. The workflow is:
1. Run the test, get a diff.
2. Inspect it — is the change explainable (different default, different sampling order)?
3. If yes: either match the package by adjusting arguments, or accept the new behaviour and regenerate the fixture, recording why in the commit message.
4. If no: investigate before merging.

---

# Plan: `raster` → `terra` migration (Phase 1)

Medium-term goal: drop the `raster` package as a dependency. Almost all of the codebase already uses `terra::` (`rast`, `mask`, `crop`, `vect`, `extract`, `predict`, `spatSample`, `writeRaster`). One holdout remains.

## Remaining `raster::` usage

- **`funx.R:182`** in `create_masked_raster()` — `raster_data <- stack(file)` returns a `RasterStack`, which is then passed to `terra::mask()` (line 192). This works through implicit class coercion in `terra` but is fragile and silently slow. Replace with `terra::rast(file)`.
- **`funx.R:199`** in the same function — `terra::writeRaster(..., format = "GTiff", overwrite = TRUE)`. `format` is a `raster::writeRaster` argument; `terra::writeRaster` uses `filetype`. Currently this argument is silently ignored. Replace with `filetype = "GTiff"` (or omit, since `.tif` extension already implies GTiff).

The two analysis scripts (`continuous_metrics_analysis.R`, `spectral_species_analysis.R`) already use `terra::` exclusively.

## Rollout

1. Edit `create_masked_raster()` to use `terra::rast()` and fix the `writeRaster` argument.
2. Drop any `library(raster)` / `requireNamespace("raster")` calls. Currently there are none in the analysis scripts, but verify.
3. Remove `raster` from any future `DESCRIPTION` / package-loading bootstrap.

## Interaction with the other plans

- Independent of the saltbush migration — `saltbush` already uses `terra::` internally based on its DESCRIPTION, so swapping to `saltbush::create_multiband_image` effectively migrates that function for free. The local `create_masked_raster` is the only function this plan acts on.
- Should land before targets refactor: pipeline targets that produce or consume `SpatRaster` objects shouldn't have to handle mixed `RasterStack`/`SpatRaster` types.

---

# Plan: reproducibility groundwork (Phase 2)

Four cheap-but-load-bearing items between Phase 1 and Phase 3. Independent of the function refactors; address the things most likely to break a reader's replication six months after publication.

Dependency pinning with `renv` is *deliberately not* in this phase — locking before saltbush/targets means constantly re-snapshotting as packages change. It lives in Phase 5 instead, once the dependency set has stabilised.

1. **Wire the rarefaction seed into the analysis scripts.** Phase 0 added `seed = NULL` to `calculate_cv()` / `calculate_chv_nopca()` / `calculate_spectral_metrics()` / `calculate_coefficient_of_variance()`, but neither analysis script passes a value — so end-user runs of `continuous_metrics_analysis.R` still produce slightly different CV/CHV every time. Pick a documented seed (e.g. `42`), pass it at each call site, mention it in the README.
2. **Add a LICENSE.** Conventional split: CC-BY-4.0 for `data/` (matches the Zenodo deposit) and MIT or Apache-2.0 for code. A single root-level LICENSE file is fine if the data/code split is documented in the README.
3. **Add CITATION.cff.** GitHub renders it as a "Cite this repository" button, and reference managers (Zotero, Mendeley) parse it. Cross-link the eventual code DOI when Phase 5 mints it.
4. **Document R version + system deps** in a SETUP section of the README (even a stub README is fine here; the full rewrite is Phase 5). `sf` and `terra` need GDAL/PROJ; `geometry` needs C++17.

Tier 1 tests pass at the end of this phase. Tier 2 passes too if the rasters are downloaded.

---

# Plan: migrate from `funx.R` to the `saltbush` package (Phase 3)

Most of the reusable functions in `funx.R` have been generalized into the `saltbush` R package (`traitecoevo/saltbush` on GitHub). The analysis should call `saltbush::` functions instead of carrying local copies.

## Mapping

`saltbush` currently exports (at least) these four functions, all with the same names as their `funx.R` counterparts:

| funx.R                          | saltbush                          | call sites                                          |
| ------------------------------- | --------------------------------- | --------------------------------------------------- |
| `create_multiband_image()`      | `saltbush::create_multiband_image()`      | (none in the two analysis scripts; used in preprocessing) |
| `extract_pixel_values()`        | `saltbush::extract_pixel_values()`        | `continuous_metrics_analysis.R` line 22             |
| `calculate_spectral_metrics()`  | `saltbush::calculate_spectral_metrics()`  | `continuous_metrics_analysis.R` line 34             |
| `calculate_field_diversity()`   | `saltbush::calculate_field_diversity()`   | both scripts                                        |

## What stays in `funx.R` for now

`saltbush` doesn't (as of writing) expose everything the analysis uses. These need to stay local until either added upstream or kept as analysis-specific helpers:

- `calculate_coefficient_of_variance()` — band-subset CV for the red-edge/NIR combos (`continuous_metrics_analysis.R` lines 134–143). Possibly a candidate for upstreaming to `saltbush` since it's a small extension of `calculate_cv`.
- `find_optimum_thresholds()` and `create_masked_raster()` — used in the preprocessing step (not in the two analysis scripts shipped here).
- `download_zenodo_rasters()` — analysis-specific; stays local.
- The low-level metric helpers (`calculate_cv`, `calculate_sv`, `calculate_chv_nopca`) — keep until verified that `saltbush::calculate_spectral_metrics` exposes the same `rarefaction = TRUE, n = 999, min_points = ...` knobs.

## Rollout

1. **Install + verify.** `remotes::install_github("traitecoevo/saltbush")`. Diff each function signature against the `funx.R` version (`args(saltbush::extract_pixel_values)` vs `args(extract_pixel_values)`). Note any differences in argument names, defaults, or return shapes.
2. **Swap one call site at a time** and confirm the output is identical. Suggested order, easiest first:
   1. `calculate_field_diversity` (pure data, no rasters needed — fastest to validate)
   2. `extract_pixel_values`
   3. `calculate_spectral_metrics`
3. **Slim `funx.R`** — delete functions that now have a `saltbush` equivalent in active use. Leave `# moved to saltbush` only if a comment is needed to explain something non-obvious; otherwise just delete.
4. **Decide on upstreaming.** If `calculate_coefficient_of_variance` is genuinely just `calculate_cv` over a band subset, consider proposing it upstream instead of keeping a local fork.

## Interaction with the other plans

- The targets pipeline (Phase 4) should call `saltbush::` functions directly, not `funx.R` copies of the same logic — so this phase must land before that one.
- Phase 4a of the targets plan ("extract model loops into functions, drop side effects") is unaffected — those are analysis-specific and won't live in `saltbush`.
- Swapping function sources is harder once a `_targets.R` is referencing them; do the dedup first.

---

# Plan: refactor to a {targets} pipeline (Phase 4)

Replace the two top-level scripts with a `targets` workflow so expensive steps (n=999 rarefaction, 20-seed RF clustering) are cached and only re-run when inputs change.

## Why targets fits

- **Expensive, deterministic steps.** Rarefaction in `calculate_cv`/`calculate_chv_nopca` and the 20-seed clustering both take a long time and produce reproducible outputs given a seed — exactly what `targets` caches well.
- **Branching is natural.** The analysis is a cross-product over `site × taxonomic_metric × spectral_metric` (and `× band_combination`, `× seed`). `tar_map()` and `tar_map_rep()` express this directly.
- **Large file inputs.** `tar_target(format = "file")` on the Zenodo download lets `targets` skip re-downloading if files are unchanged.
- **Two scripts share inputs.** `funx.R::calculate_field_diversity` is called from both; in a pipeline it's computed once.

## Prep work (before any `_targets.R` is written)

The existing functions in `funx.R` are mostly pure and reusable. The two analysis scripts mix computation with side effects — that needs untangling first. None of these should change behavior; they're pure refactors:

1. **Extract the model-fitting loops** in `continuous_metrics_analysis.R` (lines ~72–127 and ~176–229) into functions: `fit_spectral_biodiversity_model(data, tax_metric, spec_metric)` and `fit_cv_band_model(data, tax_metric, band_combo)`, each returning a one-row tibble of coefficients / CIs / R². The driver code then becomes a `purrr::pmap_dfr` (and eventually a `tar_map`).
2. **Remove the CSV write side effect** from `get_spectral_species()` in `spectral_species_analysis.R` (the `write.csv` inside the per-seed loop). Return the tibble; let `targets` handle persistence.
3. **Make spectral-species per-seed callable.** Split `get_spectral_species()` into `spectral_species_one_seed(rasters, fishnets, seed, k, sample_size)` returning one tibble. The 20-seed iteration becomes a `tar_map_rep(seeds = 1:20)`.
4. **Move every shared function into `funx.R`.** Both scripts should `source('funx.R')` (already done in both as of the Zenodo download work).
5. **Consolidate the site-prefix logic** (`E`/`G`/`S`/`C`) into a small helper so it's not duplicated.

## Target graph (sketch)

```
file targets:
  zenodo_rasters   (format="file")  -> download_zenodo_rasters()
  fishnet_files    (format="file")
  survey_csv       (format="file")

per-site (tar_map over sites = c("NSABHC0009"..0012)):
  pixel_values_<site>          <- extract_pixel_values(...)
  cv_<site>                    <- calculate_cv(...,  rarefaction = TRUE, n = 999)
  sv_<site>                    <- calculate_sv(...)
  chv_<site>                   <- calculate_chv_nopca(...)
  cv_<bands>_<site>            <- calculate_cv(...) for each of 3 band subsets

combined:
  spectral_metrics             <- bind_rows of cv/sv/chv across sites
  cv_band_combos               <- bind_rows of band-subset CVs across sites
  taxonomic_diversity          <- calculate_field_diversity(survey_csv)
  spectral_taxonomic           <- left_join(spectral_metrics, taxonomic_diversity)

models (tar_map over taxonomic × spectral metric):
  model_<tax>_<spec>           <- fit_spectral_biodiversity_model(...)
  model_results                <- bind_rows of model_<tax>_<spec>

spectral species (tar_map_rep over seeds = 1:20):
  spectral_species_seed_<s>    <- spectral_species_one_seed(...)
  mean_spectral_species        <- aggregate across seeds
  ss_models                    <- glmmTMB fits

reporting (implemented in Phase 5):
  report                       <- tar_render("report.Rmd")
```

## Sub-phases

1. **Phase 4a — prep refactors** (no targets yet). Items 1–5 above. Verify both scripts still produce identical outputs.
2. **Phase 4b — minimal pipeline.** `_targets.R` with the file inputs, taxonomic diversity, and per-site pixel extraction. Confirms the plumbing works end-to-end on the cheap stuff.
3. **Phase 4c — spectral metrics + continuous-metrics models.** Add CV/SV/CHV targets and the `tar_map` over taxonomic × spectral metrics.
4. **Phase 4d — spectral species.** Add the 20-seed branching; consider `crew`/`future` for parallel execution.

The rendered manuscript-figure report is implemented in Phase 5 (publication prep) as a `tar_render()` target within this pipeline.

## Working conventions during the transition

- Don't add new top-level scripts with side effects. Logic goes in `funx.R` (or a future `R/` directory) as named functions.
- Keep functions pure where possible: take data + parameters, return data. No `write.csv` or `dir.create` inside metric/model functions.
- Prefer adding a new function over modifying a working one — keeps existing scripts runnable as a fallback during the refactor.
- The Zenodo downloader is the canonical raster source. Don't reintroduce the legacy `data_out/combined_rasters/masked/2024/` path.

---

# Plan: publication prep (Phase 5)

Final-mile items before the code accompanies the manuscript. Most depend on the targets pipeline being in place so "run the analysis" is a one-liner.

1. **Lock dependencies with `renv`.** `renv::init()` from the project root, commit `renv.lock`, `renv/activate.R`, and the `.gitignore` lines `renv` writes. The package set is stable by this point (saltbush in, the targets pipeline driving installs) so a single snapshot captures the final state — no re-snapshotting churn. This is the single biggest reproducibility risk in an R analysis; doing it here means the lockfile that ships with the paper is the lockfile we tested against.
2. **Rewrite README.** Replace the transitional README with a real entry point: Prerequisites (R version, GDAL/PROJ, C++17 for `geometry`), Install (`renv::restore()`), Quick start (`targets::tar_make()`), Expected outputs and runtime, Data citation (Zenodo data DOI) + code citation (the new code DOI from item 7 below), License.
3. **Rendered report.** `report.Rmd` that loads pipeline outputs via `tar_read()` and produces the manuscript's tables and figures. Wired into the pipeline as a `tar_render()` target so it rebuilds whenever upstream targets change. This is the artifact that ties the code in the repo to the claims in the paper.
4. **Dockerfile.** Base off `rocker/geospatial:<R version>` (ships GDAL/PROJ pre-installed); call `renv::restore()` at build time; entrypoint runs `targets::tar_make()`. Lets reviewers and readers reproduce the analysis end-to-end without configuring a local R environment.
5. **CI for Tier 1 tests.** GitHub Actions on push: spin up R against the locked `renv` environment, `Rscript -e 'testthat::test_dir("tests")'`. Catches dependency drift early. Don't try to run Tier 2/3 in CI — the rasters are 4.6 GB and GitHub-hosted runners won't tolerate it.
6. **Output checksums as release artifacts.** Ship the Tier 1 baseline RDS (and optionally one Tier 2 fixture) alongside the GitHub release so readers can `digest::digest()` their replication output and confirm an exact match.
7. **Code DOI via Zenodo–GitHub integration.** Enable the integration; cut a `v1.0.0` release; Zenodo mints a DOI for the code, separate from the data DOI (10.5281/zenodo.17089161). Cite both in the manuscript. Update CITATION.cff with the code DOI once it exists.

Order matters: `renv` must land before items 2 (README references `renv::restore()`), 4 (Dockerfile uses `renv::restore()`), and 5 (CI runs against the locked environment).

Tier 1 + Tier 2 tests should pass against the locked environment after each item lands.
