# multispectral_drone_svh

Code supporting Gemmell et al., "Applying the spectral variability hypothesis to arid shrublands, using multispectral drone imagery."

The analysis tests whether spectral heterogeneity from drone-borne multispectral imagery (5 bands: blue, green, red, red-edge, NIR) predicts taxonomic plant diversity across four AusPlot sites in arid NSW (NSABHC0009–0012), each surveyed as a 5×5 grid of 20 m subplots.

## Repo layout

- `continuous_metrics_analysis.R` — extracts pixel values per subplot, computes coefficient of variation (CV), spectral variance (SV), and convex hull volume (CHV) using rarefaction (n=999); fits per-taxonomic-metric × per-spectral-metric mixed models (`glmmTMB`, site as random effect). Also computes CV across band subsets (red-edge+NIR, etc.).
- `spectral_species_analysis.R` — random-forest + k-means clustering of pixels into "spectral species" (k=40), repeated over 20 seeds; computes spectral richness / Shannon / Simpson per subplot; fits mixed models against field diversity.
- `funx.R` — all the reusable functions: raster I/O, masking, pixel extraction, the three spectral-metric calculators, `calculate_field_diversity` (taxonomic diversity from AusPlots transect hits), and `download_zenodo_rasters` (added recently — auto-fetches the 4 GeoTIFFs from Zenodo).
- `data/ausplots_march_24.csv` — field survey hits.
- `data/fishnets/NSABHC00{09..12}_fishnet.shp` — the 5×5 subplot grids per site.
- `data/raster_images/` — gitignored; populated on first run by `download_zenodo_rasters()` from Zenodo 10.5281/zenodo.17089161 (~4.6 GB, 4 files).
- `data_out/` — gitignored; intermediate outputs.

## Required R packages

`sf`, `terra`, `tidyverse`, `vegan`, `data.table`, `performance`, `glmmTMB`, `geometry` (CHV), `randomForest` + `cluster` (spectral species), `pROC` (mask threshold finder in `funx.R`).

## Conventions

- Per-subplot IDs are `"row_col"` (e.g. `"3_2"`); after joining across sites the analysis prefixes them with a site letter (`E`/`G`/`S`/`C`) to keep them unique.
- Bands must be stacked in wavelength order (blue, green, red, red-edge, NIR) — several functions assume this implicitly.
- Site names are the AusPlot IDs (`NSABHC0009`..`0012`) and are extracted from filenames via `str_extract(basename(x), "^[^_]+")`.

---

# Plan: migrate from `funx.R` to the `saltbush` package

Medium-term goal: most of the reusable functions in `funx.R` have been generalized into the `saltbush` R package (`traitecoevo/saltbush` on GitHub). The analysis should call `saltbush::` functions instead of carrying local copies.

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

## Phased rollout

1. **Install + verify.** `remotes::install_github("traitecoevo/saltbush")`. Diff each function signature against the `funx.R` version (`args(saltbush::extract_pixel_values)` vs `args(extract_pixel_values)`). Note any differences in argument names, defaults, or return shapes.
2. **Swap one call site at a time** and confirm the output is identical. Suggested order, easiest first:
   1. `calculate_field_diversity` (pure data, no rasters needed — fastest to validate)
   2. `extract_pixel_values`
   3. `calculate_spectral_metrics`
3. **Slim `funx.R`** — delete functions that now have a `saltbush` equivalent in active use. Leave `# moved to saltbush` only if a comment is needed to explain something non-obvious; otherwise just delete.
4. **Decide on upstreaming.** If `calculate_coefficient_of_variance` is genuinely just `calculate_cv` over a band subset, consider proposing it upstream instead of keeping a local fork.

## Interaction with the targets plan

Do the saltbush migration before (or as the first phase of) the targets refactor. Reasons:

- The targets pipeline should call `saltbush::` functions directly, not `funx.R` copies of the same logic.
- Phase 0 of the targets plan ("extract model loops into functions, drop side effects") is unaffected — those are analysis-specific and won't live in `saltbush`.
- Swapping function sources is harder once a `_targets.R` is referencing them; do the dedup first.

---

# Recommended order

Three medium-term refactors are planned. Do them in this order:

0. **Phase 0 — test scaffolding** (no behavioural change). Capture baseline outputs and seed the rarefaction. Without this every later phase risks silent output drift; with it each phase is verified by a one-command snapshot diff.
1. **`raster` → `terra` migration.** Smallest scope, isolated to one function (`create_masked_raster()`) not called by the analysis scripts. Good first refactor — proves the test harness works on a near-zero-risk change.
2. **`funx.R` → `saltbush` migration.** Highest behavioural risk but well bounded: one function at a time, verify each swap against the baseline snapshot before moving on. Easiest first: `calculate_field_diversity` → `extract_pixel_values` → `calculate_spectral_metrics`.
3. **`targets` pipeline.** Purely structural by this point. Phase 3a (factor out model loops, drop CSV side effects) verified by re-running snapshots; Phase 3b (`_targets.R`) is plumbing.

Each phase ends with `testthat::test_dir("tests")` passing.

---

# Testing strategy

No tests exist on `main` today. We're adding them as the first refactor so the next three can be verified mechanically.

## Tiers

- **Tier 1 — every change (~1 s).** Taxonomic-diversity snapshot. Uses only `data/ausplots_march_24.csv`, no rasters. Exercises `calculate_field_diversity()`, which both analysis scripts call. Fixture: `tests/fixtures/taxonomic_diversity_baseline.rds`.
- **Tier 2 — per phase, manual (~minutes).** One-site spectral-metric snapshot on NSABHC0010 (smallest raster at 554 MB). Uses `n = 99` for rarefaction in tests, `n = 999` in the real analysis. Fixture generated once, re-checked after each function swap.
- **Tier 3 — release / merge to main (~hours).** Full end-to-end on all four sites. Acceptance check, not a per-commit gate.

## Tooling

- `testthat` for structure; `waldo::compare()` for readable diffs.
- Snapshots checked into the repo as `.rds` files in `tests/fixtures/` — small (summary tibbles, not rasters).
- Project is not an R package (no `DESCRIPTION`), so tests live in flat `tests/` (no `tests/testthat/` subdir). Run with `Rscript -e 'testthat::test_dir("tests")'`.

## Reproducibility prerequisite

`calculate_cv()` and `calculate_chv_nopca()` currently call `sample()` un-seeded, so identical inputs produce non-identical outputs. Phase 0 adds an optional `seed = NULL` argument to both (and to the `calculate_spectral_metrics` / `calculate_coefficient_of_variance` wrappers). Default `NULL` preserves current behaviour; tests pass `seed = 42` or similar to get exact equality.

## When a snapshot diff appears

A refactor changing outputs isn't automatically wrong — `saltbush::calculate_spectral_metrics` may legitimately differ from the local version in defaults. The workflow is:
1. Run the test, get a diff.
2. Inspect it — is the change explainable (different default, different sampling order)?
3. If yes: either match the package by adjusting arguments, or accept the new behaviour and regenerate the fixture, recording why in the commit message.
4. If no: investigate before merging.

---

# Plan: complete the `raster` → `terra` migration

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

# Plan: refactor to a {targets} pipeline

Medium-term goal: replace the two top-level scripts with a `targets` workflow so expensive steps (n=999 rarefaction, 20-seed RF clustering) are cached and only re-run when inputs change.

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

reporting:
  report                       <- tar_render("report.Rmd")
```

## Phased rollout

1. **Phase 0 — prep refactors** (no targets yet). Items 1–5 above. Verify both scripts still produce identical outputs.
2. **Phase 1 — minimal pipeline.** `_targets.R` with the file inputs, taxonomic diversity, and per-site pixel extraction. Confirms the plumbing works end-to-end on the cheap stuff.
3. **Phase 2 — spectral metrics + continuous-metrics models.** Add CV/SV/CHV targets and the `tar_map` over taxonomic × spectral metrics.
4. **Phase 3 — spectral species.** Add the 20-seed branching; consider `crew`/`future` for parallel execution.
5. **Phase 4 — report.** `tar_render` an Rmd that produces the manuscript's tables and figures.

## Working conventions during the transition

- Don't add new top-level scripts with side effects. Logic goes in `funx.R` (or a future `R/` directory) as named functions.
- Keep functions pure where possible: take data + parameters, return data. No `write.csv` or `dir.create` inside metric/model functions.
- Prefer adding a new function over modifying a working one — keeps existing scripts runnable as a fallback during the refactor.
- The Zenodo downloader is the canonical raster source. Don't reintroduce the legacy `data_out/combined_rasters/masked/2024/` path.
