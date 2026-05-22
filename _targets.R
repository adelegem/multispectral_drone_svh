# =============================================================================
# _targets.R — pipeline definition for multispectral_drone_svh
# =============================================================================
#
#   ##############################################################
#   #                                                            #
#   #   DRAFT SKELETON — not yet wired, sharing for co-author    #
#   #   feedback before Phase 4b begins. Nothing here runs       #
#   #   end-to-end. Section markers below are the CLAUDE.md      #
#   #   phase plan.                                              #
#   #                                                            #
#   ##############################################################
#
# Why targets at all
#   The analysis is a cross-product of (site × taxonomic_metric × spectral_metric
#   × band_combination × seed). Two expensive steps dominate runtime:
#     - rarefaction at n=999 for CV/CHV
#     - 20-seed RF + k-means + per-pixel classification for spectral species
#   Both are deterministic given a seed. `targets` caches each step keyed on
#   (function body × inputs × dependencies × packages), so a one-line tweak in
#   `funx.R` invalidates only the targets that actually use the changed code.
#
# Why this skeleton has 4b wired and 4c/4d/5 stubbed
#   Phase 4b is the cheapest plumbing — file inputs, per-site pixel extraction,
#   taxonomic diversity. If that works end-to-end against the four-site dataset,
#   the rest is mechanical: each subsequent phase is "drop in another tar_map
#   over (X × Y)". Better to confirm the cheap plumbing first than design the
#   whole graph at once.
#
# To run when wired
#   install.packages(c("targets", "tarchetypes"))   # not currently installed
#   library(targets)
#   tar_make()
#   tar_visnetwork()                                # render the live DAG
#   tar_read(spectral_taxonomic_diversity)          # pull any intermediate
#
# =============================================================================
# DESIGN DECISIONS FOR CO-AUTHOR REVIEW
# =============================================================================
#
# 1. Per-site branching (tar_map over the 4 sites) for pixel extraction.
#      + each site caches independently; redownloading or remasking one raster
#        doesn't invalidate the other three
#      + obvious parallel boundary if/when we add `crew`
#      - adds 4 nearly-identical targets per phase rather than one combined
#      Decision: per-site. Keeps the masking step revisitable for a single
#      site without invalidating the other three.
#
# 2. Seed handling.
#      Two independent seeds, both reproducible — the name overlap caused
#      confusion in co-author review, so worth spelling out:
#        - RAREFACTION_SEED = 42 is a single constant used ONLY inside
#          calculate_spectral_metrics() for the n=999 rarefaction draw that
#          produces CV / SV / CHV_nopca. One seed, one draw.
#        - The 20 spectral-species seeds (1..20) drive 20 independent k-means
#          clusterings whose richness/Shannon/Simpson outputs are averaged
#          per subplot. These remain 20 distinct seeds — a single seed across
#          all 20 would defeat the averaging — expanded via
#          tar_map_rep(seeds = 1:20) in Phase 4d.
#      Decision: RAREFACTION_SEED = 42 ships with the paper; spectral-species
#      averaging stays at 20 distinct seeds (1..20).
#
# 3. Parallelism via `crew`.
#      20 spectral-species seeds × ~2 h/seed = ~40 h serial.
#      With crew_controller_local(workers = 4): ~10 h. Same memory pressure as
#      one process per worker (full-raster terra::predict is the bottleneck),
#      so workers = 4 may already pressure a 32 GB machine. Worth measuring.
#      Decision: workers = 3 on Will's 32 GB box (split the difference between
#      2-safe and 4-fast). Back off to 2 if RAM pressure shows up during the
#      first run; bump to 4 if it doesn't.
#
# 4. Where do figures live?
#      Option A: standalone tar_target(fig_5_masked_24, make_figure_5(...))
#      Option B: inline in report.Rmd, rebuilt only via tar_render() (Phase 5)
#      Decision: (B) — inline in report.Rmd. Figures and narrative move
#      together in revisions, so coupling them to the same render is cleaner
#      than maintaining standalone figure targets.
#
# 5. Singular-fit refit convention.
#      "Singular" = a variance component (here the site random-intercept
#      variance) is estimated at zero, the boundary of parameter space. For
#      pielou_evenness across our 4 sites, between-site variation is ~nil so
#      (1 | site) has nothing to estimate → glmmTMB returns a degenerate fit
#      flagged by performance::check_singularity(). Clean fix: drop the
#      random effect, refit as fixed-effect lm.
#      The 3 pielou_evenness mixed models trigger this and are currently
#      refit by hand under data_out/model_fits/*_fixed.rds. Options for the
#      targets graph:
#        - a single fit_spectral_biodiversity_model() that internally decides
#          mixed vs fixed (returns the chosen model)
#        - two separate targets (model_mixed, model_fixed) with mixed-first
#          fallback logic
#      Decision: fold into fit_spectral_biodiversity_model() — the singular
#      check + refit happens inside, caller gets one model object and a
#      flag column indicating which kind it got.
#
# 6. Manuscript report.
#      report.Rmd under tar_render() is Phase 5. The interim_progress.Rmd in
#      reports/ is the running-state log, not the manuscript artifact.
#      Decision: interim_progress.Rmd stays hand-rendered (it tracks transient
#      analysis state, not pipeline outputs); only the Phase 5 manuscript
#      report.Rmd is wired into tar_render().
#
# =============================================================================

library(targets)
library(tarchetypes)

# Packages every target sees on dispatch. Listing here (rather than relying on
# library() at the top of funx.R) means {targets} fingerprints package versions
# alongside function bodies.
tar_option_set(
  packages = c(
    "tidyverse", "data.table", "sf", "terra",
    "vegan", "performance", "glmmTMB",
    "randomForest", "saltbush"
  ),
  format = "rds"
)

# Project constants — referenced by multiple targets below.
RAREFACTION_SEED <- 42
SITES            <- c("NSABHC0009", "NSABHC0010", "NSABHC0011", "NSABHC0012")
WAVELENGTHS      <- c("blue", "green", "red", "red_edge", "nir")
BAND_SUBSETS <- list(
  red.edge_nir           = c("red_edge", "nir"),
  green_red.edge_nir     = c("green", "red_edge", "nir"),
  green_red_red.edge_nir = c("green", "red", "red_edge", "nir")
)

# Functions live in funx.R. {targets} tracks function bodies, so editing one
# function only invalidates the targets that actually use it.
source("funx.R")


# =============================================================================
# PHASE 4b — file inputs + per-site pixel extraction + taxonomic diversity
# =============================================================================
# Minimal plumbing. If `tar_make()` survives this section end-to-end on the four
# sites, the rest of the graph is mechanical.

# Per-site lookup table: drives tar_map below. Keeping the file paths explicit
# here (rather than glob-matching inside each target) means each branch's input
# files are a clean tar_target dependency.
sites_lookup <- tibble::tibble(
  site         = SITES,
  raster_file  = file.path("data/raster_images", paste0(SITES, "_masked.tif")),
  fishnet_file = file.path("data/fishnets",       paste0(SITES, "_fishnet.shp"))
)

# Per-site file + extraction tar_map.
#
# unlist = FALSE keeps the structured list of branches (accessible as
# phase_4b_branches$pixel_values, etc.) so Phase 4c's tar_combine can
# fan back in over $pixel_values without re-enumerating site names.
#
# raster_<site>:  depends on raster_files (the download gate — ensures
#                 the .tif exists on a fresh checkout) then content-tracks
#                 the single per-site .tif. A change to one raster
#                 invalidates only that site's pixel_values.
# fishnet_<site>: shapefiles are tracked in git; just content-track each
#                 .shp directly (sidecar .dbf/.shx/.prj changes are
#                 unusual enough to not warrant per-sidecar tracking).
# pixel_values_<site>: standard extraction; the saltbush wrapper handles
#                      one site at a time.
phase_4b_branches <- tar_map(
  values = sites_lookup,
  names  = site,
  unlist = FALSE,
  tar_target(raster,  { raster_files; raster_file }, format = "file"),
  tar_target(fishnet, fishnet_file,                  format = "file"),
  tar_target(pixel_values,
             extract_pixel_values(raster, fishnet, WAVELENGTHS))
)

phase_4b_targets <- list(

  # ---- file inputs --------------------------------------------------------
  # `format = "file"` tracks a file by content hash + mtime. Subsequent
  # tar_make() calls skip if unchanged.
  tar_target(survey_csv, "data/ausplots_march_24.csv", format = "file"),

  # Zenodo downloader gates the per-site raster file targets above. It's
  # idempotent (skips files already present at expected size) and returns
  # the four paths. NOT format = "file" on purpose — content-tracking the
  # 4-file bundle would invalidate every per-site branch on a single
  # raster change. Per-site `raster_<site>` file targets handle that.
  tar_target(raster_files,
             download_zenodo_rasters("data/raster_images")),

  # ---- taxonomic diversity -----------------------------------------------
  # Single fast target — no per-site branching needed. Reads the CSV at the
  # boundary so `survey_csv` is the file dependency.
  tar_target(taxonomic_diversity,
             calculate_field_diversity(
               readr::read_csv(survey_csv, show_col_types = FALSE)
             )$final_results),

  phase_4b_branches
)


# =============================================================================
# PHASE 4c — spectral metrics + continuous-metrics models
# =============================================================================
# Mirrors continuous_metrics_analysis.R: rarefied spectral metrics, CV across
# band subsets, joined with taxonomic diversity, then the 12 + 12 fixed/mixed
# models. spectral_metrics is a single target rather than per-site because the
# n = 999 rarefaction draw is RNG-order-sensitive — splitting per site would
# change the published numbers.
#
# Per-pair model targets hold the fitted glmmTMB / lm object so Phase 5
# figures can pull cached fits via tar_read() for predicted lines. Per-pair
# summary targets are tar_combine'd into the model_results tibbles.

# Per-band-subset CV tar_map. Uses a list-column for the wavelength vector;
# `band_vec` substitutes as the literal R expression, `bands_label` as a
# string — avoids the name clash that occurs when a tar_map column name
# matches a function arg or output column name.
cv_branches <- tar_map(
  values = tibble::tibble(
    bands_label = names(BAND_SUBSETS),
    band_vec    = unname(BAND_SUBSETS)
  ),
  names  = bands_label,
  unlist = FALSE,
  tar_target(
    cv_subset,
    calculate_coefficient_of_variance(
      pixel_values,
      wavelengths = band_vec,
      min_points  = min_points,
      n           = 999,
      rarefaction = TRUE,
      seed        = RAREFACTION_SEED
    ) %>% dplyr::mutate(bands = bands_label)
  )
)

# Per-pair model + summary tar_map for spectral metrics models.
# tar_map names: model_<tax>_<spec>, summary_<tax>_<spec>.
model_branches <- tar_map(
  values = tidyr::expand_grid(
    tax_metric  = c("species_richness", "exp_shannon",
                    "inv_simpson", "pielou_evenness"),
    spec_metric = c("CV", "SV", "log.CHV")
  ),
  names  = c(tax_metric, spec_metric),
  unlist = FALSE,
  tar_target(model,
             fit_spectral_biodiversity_model(
               spectral_taxonomic_diversity, tax_metric, spec_metric)),
  tar_target(summary,
             summarise_spectral_biodiversity_model(
               model, tax_metric, spec_metric))
)

# Per-pair model + summary tar_map for CV band-combination models.
# tar_map names: cv_model_<tax>_<bands>, cv_summary_<tax>_<bands>.
cv_model_branches <- tar_map(
  values = tidyr::expand_grid(
    tax_metric = c("species_richness", "exp_shannon",
                   "inv_simpson", "pielou_evenness"),
    band_combo = names(BAND_SUBSETS)
  ),
  names  = c(tax_metric, band_combo),
  unlist = FALSE,
  tar_target(cv_model,
             fit_cv_band_model(cv_values, tax_metric, band_combo)),
  tar_target(cv_summary,
             summarise_cv_band_model(cv_model, tax_metric, band_combo))
)

phase_4c_targets <- list(

  # ---- combined pixel_values --------------------------------------------
  # Bind the four per-site pixel_values into one tibble. spectral_metrics
  # below consumes this in a single rarefaction call so the published RNG
  # sequence is preserved.
  tar_combine(pixel_values,
              phase_4b_branches$pixel_values,
              command = dplyr::bind_rows(!!!.x)),

  # min_points: smallest non-NA pixel count across all (site, subplot)
  # combinations — drives the rarefaction sample size. Previously hardcoded
  # at 133609 in the script; derived from data here.
  tar_target(min_points,
             pixel_values %>%
               na.omit() %>%
               dplyr::count(site, subplot_id) %>%
               dplyr::pull(n) %>%
               min()),

  # ---- spectral metrics (CV / SV / CHV_nopca / log.CHV) ------------------
  # Single call to preserve the script's rarefaction RNG sequence. Adds
  # log.CHV for the manuscript scale. subplot_id stays unprefixed until
  # after the taxonomic join — taxonomic_diversity uses unprefixed
  # subplot_id, and prefixing here would break the join.
  tar_target(
    spectral_metrics,
    calculate_spectral_metrics(
      pixel_values,
      wavelengths = WAVELENGTHS,
      min_points  = min_points,
      n           = 999,
      rarefaction = TRUE,
      seed        = RAREFACTION_SEED
    ) %>%
      dplyr::mutate(log.CHV = log(CHV_nopca))
  ),

  # ---- CV per band subset + combined cv_values ---------------------------
  # cv_band_combos gets prefixed here (no further taxonomic join precedes
  # the rbind into cv_values below).
  cv_branches,
  tar_combine(cv_band_combos,
              cv_branches$cv_subset,
              command = dplyr::bind_rows(!!!.x) %>% add_site_prefix()),

  # ---- spectral × taxonomic join ----------------------------------------
  # Join unprefixed, then prefix the result — matches the legacy script
  # (continuous_metrics_analysis.R lines 59–62).
  tar_target(spectral_taxonomic_diversity,
             dplyr::left_join(spectral_metrics, taxonomic_diversity,
                              by = c("site", "subplot_id")) %>%
               add_site_prefix()),

  # ---- cv_values: band-subset CVs + all-bands CV + taxonomic join -------
  # Matches the script's `cv_values` structure (one row per site × subplot ×
  # bands with taxonomic metrics joined). Used by the CV-band models below.
  tar_target(
    cv_values,
    {
      all_bands <- spectral_taxonomic_diversity %>%
        dplyr::select(site, subplot_id, CV) %>%
        dplyr::mutate(bands = "all_bands")

      combined <- dplyr::bind_rows(all_bands, cv_band_combos)

      spectral_taxonomic_diversity %>%
        dplyr::select(species_richness, exp_shannon, inv_simpson,
                      pielou_evenness, site, subplot_id) %>%
        dplyr::left_join(combined, by = c("site", "subplot_id"))
    }
  ),

  # ---- per-pair model fits + tar_combine'd summary tables ---------------
  model_branches,
  tar_combine(spectral_biodiversity_model_results,
              model_branches$summary,
              command = dplyr::bind_rows(!!!.x)),

  cv_model_branches,
  tar_combine(cv_biodiversity_model_results,
              cv_model_branches$cv_summary,
              command = dplyr::bind_rows(!!!.x))
)


# =============================================================================
# PHASE 4d — 20-seed spectral species
# =============================================================================
# spectral_species_one_seed is RNG-deterministic given a seed, so per-seed
# branches are independent and trivially cacheable. Serial runtime is ~40 h
# (~2 h per seed, measured from data_out/spectral_species_seeds checkpoints
# May 20-21 2026). With crew_controller_local(workers = N) the seeds run in
# parallel; not configured here because each terra::predict on a 1.5 GB
# raster has real RAM cost and the right N depends on the host.
#
# The SS models (richness, Shannon, Simpson) pair each spectral metric with
# its matching taxonomic metric — 3 models, not the 4x3 cross-product of
# Phase 4c — and use raw (unscaled) predictors to match
# spectral_species_analysis.R lines 82, 87, 92.

# Per-seed tar_map. Each branch reads the per-site raster + fishnet file
# paths from the Phase 4b file targets — terra::SpatRasters don't serialise
# across target boundaries, so spectral_species_one_seed() opens them
# itself given paths.
ss_seed_branches <- tar_map(
  values = tibble::tibble(seed = 1:20),
  names  = seed,
  unlist = FALSE,
  tar_target(
    spectral_species_seed,
    spectral_species_one_seed(
      image_paths   = ss_image_paths,
      fishnet_paths = ss_fishnet_paths,
      seed          = seed,
      k_clusters    = 40,
      sample_size   = 1250
    )
  )
)

# Per-pair model + summary tar_map for the three SS models. tar_map names:
# ss_model_<tax>_<ss>, ss_summary_<tax>_<ss>. scale_predictor = FALSE matches
# the legacy spectral_species_analysis.R (no scale() on the predictor).
ss_model_branches <- tar_map(
  values = tibble::tribble(
    ~tax_metric,         ~ss_metric,
    "species_richness",  "spectral_richness",
    "shannon_diversity", "shannon_spectral",
    "simpson_diversity", "simpson_spectral"
  ),
  names  = c(tax_metric, ss_metric),
  unlist = FALSE,
  tar_target(ss_model,
             fit_spectral_biodiversity_model(
               mean_spectral_species, tax_metric, ss_metric,
               scale_predictor = FALSE)),
  tar_target(ss_summary,
             summarise_spectral_biodiversity_model(
               ss_model, tax_metric, ss_metric,
               scale_predictor = FALSE))
)

phase_4d_targets <- list(

  # ---- file-path vectors fed to spectral_species_one_seed --------------
  # tar_combine over the Phase 4b file targets so this stays in sync with
  # the per-site raster/fishnet content tracking — any raster change
  # invalidates downstream seeds correctly.
  tar_combine(ss_image_paths,
              phase_4b_branches$raster,
              command = c(!!!.x)),
  tar_combine(ss_fishnet_paths,
              phase_4b_branches$fishnet,
              command = c(!!!.x)),

  # ---- 20 per-seed spectral-species tibbles ----------------------------
  ss_seed_branches,

  # ---- combined long table across seeds --------------------------------
  tar_combine(spectral_species,
              ss_seed_branches$spectral_species_seed,
              command = dplyr::bind_rows(!!!.x)),

  # ---- per-subplot mean across the 20 seeds + taxonomic join -----------
  # Matches spectral_species_analysis.R lines 60-78: group, summarise,
  # left_join with taxonomic_diversity, then add the E/G/S/C prefix
  # (same convention as spectral_taxonomic_diversity above).
  tar_target(
    mean_spectral_species,
    spectral_species %>%
      dplyr::group_by(subplot_id, site) %>%
      dplyr::summarise(
        spectral_richness = mean(spectral_richness),
        shannon_spectral  = mean(shannon_spectral),
        simpson_spectral  = mean(simpson_spectral),
        .groups = "drop"
      ) %>%
      dplyr::left_join(taxonomic_diversity,
                       by = c("site", "subplot_id")) %>%
      add_site_prefix()
  ),

  # ---- per-pair model fits + tar_combine'd summary table --------------
  ss_model_branches,
  tar_combine(spectral_species_model_results,
              ss_model_branches$ss_summary,
              command = dplyr::bind_rows(!!!.x))
)


# =============================================================================
# PHASE 5 — figures + manuscript report  (STUB)
# =============================================================================
#
# Likely inlined into report.Rmd rather than standalone targets (see decision
# #4). report.Rmd uses tar_read() to pull spectral_taxonomic_diversity,
# model_results, cv_band_results, model_fits (a list of the per-tax × per-spec
# fitted objects, for Figure 5's predicted lines), mean_spectral_species, and
# ss_models. Reformatted into the manuscript's Tables 1–2 and Figures 5–6.
#
#   report                       tar_render("report.Rmd")


# =============================================================================
# Return value — `targets` reads the list at the bottom of _targets.R.
# Currently Phases 4b + 4c + 4d.
# =============================================================================
list(phase_4b_targets, phase_4c_targets, phase_4d_targets)
