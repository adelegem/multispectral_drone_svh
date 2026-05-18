library(testthat)

# Tier 2 — Zenodo-raster processing for one site (NSABHC0010, smallest).
# Skipped by default; opt in with RUN_TIER2=true.
# Also skips if the raster hasn't been downloaded or if `geometry` is missing.

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(data.table)
  library(sf)
  library(terra)
})

project_root <- normalizePath(file.path(dirname(testthat::test_path()), ".."))
source(file.path(project_root, "funx.R"))

site         <- "NSABHC0010"
raster_path  <- file.path(project_root, "data", "raster_images",
                          paste0(site, "_masked.tif"))
fishnet_path <- file.path(project_root, "data", "fishnets",
                          paste0(site, "_fishnet.shp"))
wavelengths  <- c("blue", "green", "red", "red_edge", "nir")

skip_unless_tier2_ready <- function() {
  if (Sys.getenv("RUN_TIER2", unset = "false") != "true") {
    skip("Tier 2 test (set RUN_TIER2=true to run)")
  }
  if (!file.exists(raster_path)) {
    skip(paste("raster not downloaded:", raster_path))
  }
  if (!file.exists(fishnet_path)) {
    skip(paste("fishnet missing:", fishnet_path))
  }
}

baseline_path <- function(name) {
  file.path(project_root, "tests", "fixtures",
            paste0(name, "_baseline_", site, ".rds"))
}

# Cache the (slow) pixel extraction across tests in the same run
pixels_cache <- NULL
get_pixels <- function() {
  if (is.null(pixels_cache)) {
    pixels_cache <<- extract_pixel_values(raster_path, fishnet_path, wavelengths)
  }
  pixels_cache
}

test_that("extract_pixel_values summary is stable (NSABHC0010)", {
  skip_unless_tier2_ready()

  pixels <- get_pixels()

  # Snapshot summary stats per subplot rather than the millions of raw pixel
  # rows: pixel count + per-band mean. Catches drift in spatial clipping or
  # band ordering without a huge fixture.
  pixel_summary <- pixels %>%
    dplyr::group_by(site, subplot_id) %>%
    dplyr::summarise(
      n_pixels = dplyr::n(),
      dplyr::across(dplyr::all_of(wavelengths), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(subplot_id)

  fixture <- baseline_path("pixel_summary")
  if (!file.exists(fixture)) {
    saveRDS(pixel_summary, fixture)
    skip("baseline created; rerun to compare")
  }
  expect_equal(pixel_summary, readRDS(fixture))
})

test_that("calculate_spectral_metrics output is stable (NSABHC0010)", {
  skip_unless_tier2_ready()
  skip_if_not_installed("geometry")

  pixels <- get_pixels()

  # Use a modest min_points and small n for test speed (real analysis uses
  # 133609 and n=999). Seeded so rarefaction is reproducible.
  min_pixels_per_subplot <- pixels |>
    dplyr::group_by(subplot_id) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::pull(n) |>
    min()

  metrics <- calculate_spectral_metrics(
    pixels,
    wavelengths = wavelengths,
    min_points  = min(min_pixels_per_subplot, 5000L),
    n           = 99,
    rarefaction = TRUE,
    seed        = 42
  ) %>%
    dplyr::arrange(subplot_id)

  fixture <- baseline_path("spectral_metrics")
  if (!file.exists(fixture)) {
    saveRDS(metrics, fixture)
    skip("baseline created; rerun to compare")
  }
  expect_equal(metrics, readRDS(fixture))
})
