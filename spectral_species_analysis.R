# SPECTRAL SPECIES ANALYSIS

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(vegan)  # for diversity()
library(randomForest)
library(cluster)
library(glmmTMB)
library(performance)
source('funx.R')

# Survey data drives the taxonomic-diversity join near the end of this script.
# Loaded up-front so the script runs stand-alone (not just after sourcing
# continuous_metrics_analysis.R in the same session).
survey_data <- read_csv('data/ausplots_march_24.csv', show_col_types = FALSE)

# site names
sites <- c("NSABHC0009", "NSABHC0010", "NSABHC0011", "NSABHC0012")

# paths — rasters are gitignored; download_zenodo_rasters() fetches them from
# Zenodo (10.5281/zenodo.17089161) into data/raster_images/ if not already present.
download_zenodo_rasters("data/raster_images")
image_paths <- paste0("data/raster_images/", sites, "_masked.tif")

# fishnet paths
fishnet_paths <- paste0("data/fishnets/", sites, "_fishnet.shp")

# 20 seeds
twentyseeds <- c(1:20)

# output dir to save each iteration in case of failing
output_dir <- 'data_out/spectral_species_seeds'

# Per-seed checkpoint dir created on first write. Phase 4d will replace this
# loop with tar_map_rep(seeds = 1:20); per-seed persistence then comes from
# tar_target caching, and this manual checkpoint goes away.
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

all_results <- list()
for (seed in twentyseeds) {
  cat("Running seed:", seed, "\n")

  seed_result <- spectral_species_one_seed(
    image_paths, fishnet_paths, seed,
    k_clusters = 40, sample_size = 1250
  )

  # Crash-recovery checkpoint per seed (resumeable if R dies mid-run)
  write.csv(seed_result,
            file = file.path(output_dir,
                             paste0("spectral_species_seed_", seed, ".csv")),
            row.names = FALSE)

  all_results[[as.character(seed)]] <- seed_result
}

spectral_species <- do.call(rbind, all_results)

# mean spectral species per subplot
mean_spectral_species <- spectral_species %>%
  group_by(subplot_id, site) %>%
  summarise(spectral_richness = mean(spectral_richness),
            shannon_spectral = mean(shannon_spectral),
            simpson_spectral = mean(simpson_spectral))

# assign site
mean_spectral_species$site <- substr(mean_spectral_species$site, 1, 10)

# load taxonomic diersity
taxonomic_diversity <- calculate_field_diversity(survey_data)$final_results %>%
  select(species_richness, shannon_diversity, simpson_diversity, subplot_id, site)

# join
mean_spectral_species <- mean_spectral_species %>%
  left_join(taxonomic_diversity, by = c('site', 'subplot_id'))

# MODELS
# spectral-species richness
sr_model <- glmmTMB(species_richness ~ spectral_richness + (1 | site), data = mean_spectral_species)
summary(sr_model)
r2_nakagawa(sr_model)
rmse(sr_model)

# spectral shannons diversity
sh_model <- glmmTMB(shannon_diversity ~ shannon_spectral + (1 | site), data = mean_spectral_species)
summary(sh_model)
r2_nakagawa(sh_model)
rmse(sh_model)

# spectral simpsons diversity
si_model <- glmmTMB(simpson_diversity ~ simpson_spectral + (1 | site), data = mean_spectral_species)
summary(si_model)
r2_nakagawa(si_model)
rmse(si_model)

# Persist outputs for downstream reporting. data_out/ is gitignored.
if (!dir.exists("data_out")) dir.create("data_out", recursive = TRUE)
saveRDS(spectral_species,      "data_out/spectral_species.rds")
saveRDS(mean_spectral_species, "data_out/mean_spectral_species.rds")
saveRDS(sr_model,              "data_out/sr_model.rds")
saveRDS(sh_model,              "data_out/sh_model.rds")
saveRDS(si_model,              "data_out/si_model.rds")

spectral_species_rmse <- tibble::tibble(
  model    = c("species_richness", "shannon_diversity", "simpson_diversity"),
  response = c("spectral_richness", "shannon_spectral", "simpson_spectral"),
  rmse     = c(rmse(sr_model), rmse(sh_model), rmse(si_model))
)
saveRDS(spectral_species_rmse, "data_out/spectral_species_rmse.rds")

