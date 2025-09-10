# SPECTRAL SPECIES ANALYSIS

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(vegan)  # for diversity()
library(randomForest)
library(cluster)
library(glmmTMB)
library(performance)

# site names
sites <- c("NSABHC0009", "NSABHC0010", "NSABHC0011", "NSABHC0012")

# paths
image_paths <- paste0("data_out/combined_rasters/masked/2024/", sites, "_combined_image_masked.tif")

# fishnet paths
fishnet_paths <- paste0("data/fishnets/", sites, "_fishnet.shp")

# 20 seeds
twentyseeds <- c(1:20)

# output dir to save each iteration in case of failing
output_dir <- 'data_out/spectral_species_seeds'

# spectral species function
get_spectral_species <- function(image_paths, fishnet_paths, seeds, k_clusters = 40, sample_size = 2500, output_dir){

  rasters <- lapply(image_paths, terra::rast)

  all_results <- list()

  for (seed in seeds){
    set.seed(seed)

    cat("Running seed:", seed, "\n")

    # sample pixels per raster, combine
    sample_list <- lapply(rasters, function(r) {
      spatSample(r, size = sample_size, method = 'random', na.rm = TRUE, as.df = TRUE)
    })

    combined_samples <- do.call(rbind, sample_list)

    # train RF with proximity matrix
    rf_model <- randomForest(x = combined_samples, proximity = TRUE, ntree = 500)
    prox_matrix <- rf_model$proximity

    # k-means clustering
    km <- kmeans(prox_matrix, centers = k_clusters, nstart = 10)
    combined_samples$cluster <- km$cluster

    # train RF classifier to classify all pixels based on clusters
    rf_cluster_model <- randomForest(x = combined_samples[, 1:(ncol(combined_samples)-1)],
                                     y = as.factor(combined_samples$cluster), ntree = 500)

    seed_results <- list()

    for (i in seq_along(rasters)){
      r <- rasters[[i]]
      fishnet_path <- fishnet_paths[[i]]

      # predict spectral species for full raster
      r_cluster <- terra::predict(r, rf_cluster_model, type = "response", na.rm = TRUE)

      # read subplot shapefile and convert to SpatVector
      subplot <- read_sf(fishnet_path) %>%
        dplyr::select(geometry) %>%
        dplyr::mutate(subplot_id = unlist(lapply(1:5, function(x) paste(x, 1:5, sep = "_")))) %>%
        vect()

      # extract cluster classification values for subplots
      cluster_vals <- terra::extract(r_cluster, subplot)

      # confirm name of 2 is cluster
      colnames(cluster_vals)[2] <- "cluster"

      # assign subplot_id to each extracted pixel by matching ID to subplot features
      cluster_vals$subplot_id <- subplot$subplot_id[cluster_vals$ID]

      # filter out NA clusters
      cluster_vals <- cluster_vals %>% filter(!is.na(cluster))

      # spectral species richness per subplot
      spectral_richness <- cluster_vals %>%
        group_by(subplot_id) %>%
        summarise(spectral_species_richness = n_distinct(cluster), .groups = "drop")

      # community matrix for Shannon and Simpson diversity
      community_matrix <- cluster_vals %>%
        group_by(subplot_id, cluster) %>%
        summarise(count = n(), .groups = "drop") %>%
        pivot_wider(names_from = cluster, values_from = count, values_fill = 0) %>%
        tibble::column_to_rownames("subplot_id")

      result_df <- data.frame(
        subplot_id = rownames(community_matrix),
        spectral_richness = spectral_richness$spectral_species_richness[match(rownames(community_matrix), spectral_richness$subplot_id)],
        shannon_spectral = diversity(community_matrix, index = "shannon"),
        simpson_spectral = diversity(community_matrix, index = "simpson"),
        seed = seed,
        site = substr(basename(fishnet_path), 1, 10)
      )

      seed_results[[i]] <- result_df
    }

    # combine results for all rasters for this seed
    combined_seed_results <- do.call(rbind, seed_results)

    # save the results to CSV
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write.csv(combined_seed_results, file = file.path(output_dir, paste0("spectral_species_seed_", seed, ".csv")), row.names = FALSE)

    # store in all_results
    all_results[[as.character(seed)]] <- combined_seed_results

  }

  # combine all seeds results
  do.call(rbind, all_results)
}

# get spectral species
spectral_species <- get_spectral_species(image_paths, fishnet_paths, twentyseeds, 40, 1250, output_dir)

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
summary(model)
r2_nakagawa(model)

# spectral shannons diversity
sh_model <- glmmTMB(shannon_diversity ~ shannon_spectral+ (1 | site), data = mean_spectral_species)
summary(model)
r2_nakagawa(model)

# spectral simpsons diversity
si_model <- glmmTMB(simpson_diversity ~ simpson_spectral + (1 | site), data = mean_spectral_species)
summary(model)
r2_nakagawa(model)

