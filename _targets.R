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

phase_4b_targets <- list(

  # ---- file inputs --------------------------------------------------------
  # `format = "file"` tracks a file by content hash + mtime. Subsequent
  # tar_make() calls skip if unchanged.

  tar_target(survey_csv, "data/ausplots_march_24.csv", format = "file"),

  # Zenodo downloader gates the per-site raster file targets below. It's
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

  # ---- per-site file + pixel extraction ----------------------------------
  # tar_map expands this block into one set of targets per row of
  # sites_lookup:
  #   raster_NSABHC0009, fishnet_NSABHC0009, pixel_values_NSABHC0009, ...
  #
  # raster_<site>:  depends on raster_files (the download gate — ensures
  #                 the .tif exists on a fresh checkout) then content-tracks
  #                 the single per-site .tif. A change to one raster
  #                 invalidates only that site's pixel_values.
  # fishnet_<site>: shapefiles are tracked in git; just content-track each
  #                 .shp directly (sidecar .dbf/.shx/.prj changes are
  #                 unusual enough to not warrant per-sidecar tracking).
  # pixel_values_<site>: standard extraction; saltbush::extract_pixel_values
  #                      under the wrapper handles one site at a time.
  tar_map(
    values = sites_lookup,
    names  = site,        # use only `site` for branch suffixes
    tar_target(raster,
               { raster_files; raster_file },
               format = "file"),
    tar_target(fishnet, fishnet_file, format = "file"),
    tar_target(pixel_values,
               extract_pixel_values(raster, fishnet, WAVELENGTHS))
  )

  # tar_combine() to bind the four pixel_values_<site> targets back into one
  # tibble would slot in here, once we're ready to consume them in 4c.
)


# =============================================================================
# PHASE 4c — spectral metrics + continuous-metrics models  (STUB)
# =============================================================================
# Wires the existing continuous_metrics_analysis.R flow into targets:
#
#   spectral_metrics_<site>      tar_map over SITES; calls
#                                calculate_spectral_metrics(..., seed = 42,
#                                n = 999) returning CV/SV/CHV_nopca
#
#   spectral_metrics             tar_combine across sites + log.CHV mutation
#                                + add_site_prefix()
#
#   cv_band_combos_<bands>_<site>
#                                tar_map over BAND_SUBSETS × SITES; calls
#                                calculate_coefficient_of_variance()
#
#   cv_band_combos               tar_combine across sites and band subsets
#
#   spectral_taxonomic_diversity left_join(spectral_metrics, taxonomic_diversity)
#
#   model_<tax>_<spec>           tar_map over (taxonomic_metric × spectral_metric)
#                                calling fit_spectral_biodiversity_model();
#                                each branch holds the fitted glmmTMB / lm
#                                object (for figures + post-hoc diagnostics)
#
#   model_results                tar_combine; one-row-per-fit summary tibble
#
#   cv_band_results_<bands>_<tax> tar_map over (band_combo × taxonomic_metric)
#                                 calling fit_cv_band_model()
#
#   cv_band_results              tar_combine to one tibble
#
# Open: see design decision #5 (where the singular-refit lives).


# =============================================================================
# PHASE 4d — 20-seed spectral species  (STUB)
# =============================================================================
#
#   spectral_species_seed_<s>    tar_map_rep(seeds = 1:20) calling
#                                spectral_species_one_seed(image_paths,
#                                fishnet_paths, seed, k_clusters = 40,
#                                sample_size = 1250)
#
#   spectral_species             tar_combine across seeds
#
#   mean_spectral_species        group_by(subplot_id, site) %>% summarise(mean)
#                                joined with taxonomic_diversity
#
#   ss_models (sr_model, sh_model, si_model)
#                                three glmmTMB fits against mean_spectral_species
#
# Cost: ~40 h serial. With crew_controller_local(workers = 4): ~10 h, modulo
# RAM (single-process peak is dominated by terra::predict on a 1.5 GB raster).


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
# Right now this is Phase 4b only.
# =============================================================================
phase_4b_targets
