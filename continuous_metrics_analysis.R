# load libraries
library(sf)
library(terra)
library(tidyverse)
library(vegan)
library(data.table)
library(performance)
library(glmmTMB)
source('funx.R')

# Seed used for all rarefaction draws below. Documented so readers can
# reproduce the published CV / CHV values exactly.
RAREFACTION_SEED <- 42

# CALCULATE SPECTRAL METRICS

# load subplot files
subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

# load files

# raster images are too large to be loaded to github (gitignored locally).
# download_zenodo_rasters() fetches them from Zenodo (10.5281/zenodo.17089161)
# into data/raster_images/, skipping any that are already present.
raster_files <- download_zenodo_rasters("data/raster_images")

pixel_values <- extract_pixel_values(raster_files, subplot_files, c('blue', 'green', 'red', 'red_edge', 'nir'))

# calculate min pixels per subplot for rarefraction
pixel_per_plot <- pixel_values %>%
  na.omit() %>%
  group_by(site, subplot_id) %>%
  summarise(count = n())

min(pixel_per_plot$count)
# 133609

# calculate spectral metrics
spectral_metrics <- calculate_spectral_metrics(pixel_values,
                                         c('blue', 'green', 'red', 'red_edge', 'nir'),
                                         rarefaction = T,
                                         n = 999,
                                         min_points = 133609,
                                         seed = RAREFACTION_SEED)

# log transform CHV to meet model reqs
spectral_metrics$log.CHV <- log(spectral_metrics$CHV_nopca)

# CALCULATE TAXONOMIC METRICS

# read in survey data
survey_data <- read_csv('data/ausplots_march_24.csv')


# calculate taxonomic diversity
taxonomic_diversity <- calculate_field_diversity(survey_data)$final_results

# join taxonomic and spectral diversity
spectral_taxonomic_diversity <- left_join(spectral_metrics, taxonomic_diversity, by = c('site', 'subplot_id'))

# unique subplot_id across sites (prefixes "row_col" with E/G/S/C)
spectral_taxonomic_diversity <- add_site_prefix(spectral_taxonomic_diversity)

# define spectral and taxonomic metrics
taxonomic_metrics <- c("species_richness", "exp_shannon", "inv_simpson", "pielou_evenness")
spectral_metrics <- c("CV", "SV", "log.CHV")

# linear mixed-effect models: <tax> ~ scale(<spec>) + (1 | site).
# expand_grid varies the last column fastest, so this iterates spec
# within tax, matching the prior nested-loop row order.
spectral_biodiversity_model_results <- expand_grid(
  tax_metric  = taxonomic_metrics,
  spec_metric = spectral_metrics
) %>%
  pmap_dfr(fit_spectral_biodiversity_model, data = spectral_taxonomic_diversity)

# CALCULATE CV FOR BAND COMBINATIONS + NDVI

# calculate CV based on band combinations

# red edge and nir
red_edge_nir_cv <- calculate_coefficient_of_variance(pixel_values, wavelengths = c('red_edge', 'nir'), min_points = 133609, rarefaction = T, seed = RAREFACTION_SEED)
red_edge_nir_cv$bands <- 'red.edge_nir'

# green, red-edge, nir
green_re_nir_cv <- calculate_coefficient_of_variance(pixel_values, wavelengths = c('green', 'red_edge', 'nir'), min_points = 133609, rarefaction = T, seed = RAREFACTION_SEED)
green_re_nir_cv$bands <- 'green_red.edge_nir'

# green, red, red-edge, nir
veg_bands_cv <- calculate_coefficient_of_variance(pixel_values, wavelengths = c('green', 'red', 'red_edge', 'nir'), min_points = 133609, rarefaction = T, seed = RAREFACTION_SEED)
veg_bands_cv$bands <- 'green_red_red.edge_nir'

cv_values <- rbind(red_edge_nir_cv, veg_bands_cv, green_re_nir_cv)

# unique subplot id across sites (E/G/S/C prefix)
cv_values <- add_site_prefix(cv_values)

# add five band dataset CV values
cv_values <- spectral_taxonomic_diversity %>%
  select(site, subplot_id, CV) %>%
  mutate(bands = 'all_bands') %>%
  rbind(cv_values)

# add taxonomic metrics
cv_values <- spectral_taxonomic_diversity %>%
  select(species_richness, exp_shannon, inv_simpson, pielou_evenness, site, subplot_id) %>%
  left_join(cv_values, by = c('site', 'subplot_id'))

# MODELS

taxonomic_metrics <- c("species_richness", "exp_shannon", "inv_simpson", "pielou_evenness")
band_combinations <- c("red.edge_nir", "green_red.edge_nir", "green_red_red.edge_nir")

# Outer = band_combo, inner = tax_metric — last column varies fastest under
# expand_grid, matching the prior nested-loop row order.
cv_biodiversity_model_results <- expand_grid(
  band_combo = band_combinations,
  tax_metric = taxonomic_metrics
) %>%
  pmap_dfr(fit_cv_band_model, data = cv_values)

# Persist outputs for downstream reporting. data_out/ is gitignored.
if (!dir.exists("data_out")) dir.create("data_out", recursive = TRUE)
saveRDS(spectral_taxonomic_diversity,        "data_out/spectral_taxonomic_diversity.rds")
saveRDS(spectral_biodiversity_model_results, "data_out/spectral_biodiversity_model_results.rds")
saveRDS(cv_values,                           "data_out/cv_values.rds")
saveRDS(cv_biodiversity_model_results,       "data_out/cv_biodiversity_model_results.rds")

