# load libraries
library(sf)
library(terra)
library(tidyverse)
library(vegan)
library(data.table)
library(performance)
library(glmmTMB)
source('funx.R')

# CALCULATE SPECTRAL METRICS

# load subplot files
subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

# load files

# raster images are too large to be loaded to github but can be found on zenodo
# 10.5281/zenodo.17089161
raster_files <- list.files('data/raster_images', pattern = '_combined_image_masked.tif$', full.names = TRUE)

pixel_values <- extract_pixel_values(raster_files, subplot_files, c('blue', 'green', 'red', 'red_edge', 'nir'))

# calculate min pixels per subplot for rarefraction
pixel_per_plot <- pixel_values %>%
  na.omit() %>%
  group_by(identifier, subplot_id) %>%
  summarise(count = n())

min(pixel_per_plot$count)
# 133609

# calculate spectral metrics
spectral_metrics <- calculate_spectral_metrics(pixel_values,
                                         c('blue', 'green', 'red', 'red_edge', 'nir'),
                                         rarefaction = T,
                                         n = 999,
                                         min_points = 133609)

# log transform CHV to meet model reqs
spectral_metrics$log.CHV <- log(spectral_metrics$CHV)

# CALCULATE TAXONOMIC METRICS

# read in survey data
survey_data <- read_csv('data/ausplots_march_24.csv')


# calculate taxonomic diversity
taxonomic_diversity <- calculate_field_diversity(survey_data)$final_results

# join taxonomic and spectral diversity
spectral_taxonomic_diversity <- left_join(spectral_metrics, field_diversity, by = c('site', 'subplot_id'))

# unique subplot_id for different sites
spectral_taxonomic_diversity <- spectral_taxonomic_diversity %>%
  mutate(subplot_id = case_when(
  site == 'NSABHC0009' ~ paste0('E',subplot_id),
  site == 'NSABHC0010' ~ paste0('G',subplot_id),
  site == 'NSABHC0011' ~ paste0('S',subplot_id),
  site == 'NSABHC0012' ~ paste0('C',subplot_id)
))

# define spectral and taxonomic metrics
taxonomic_metrics <- c("species_richness", "exp_shannon", "inv_simpson", "pielou_evenness")
spectral_metrics <- c("CV", "SV", "log.CHV")

# store empty data frame
spectral_biodiversity_model_results <- data.frame()

# linear mixed-effect models
for (tax in taxonomic_metrics) {
  for (spec in spectral_metrics) {

    # fixed effect = spectral metric, random effect = site
    # model = taxonomic_metric ~ scale(spectral_metric) + (1 + site)
    model <- try(glmmTMB(
      formula = as.formula(paste(tax, "~ scale(", spec, ") +(1 | site)")),
      data = spectral_taxonomic_diversity
    ), silent = TRUE)

    # confirm convergence + check singularity for model where random effect negligible (pielous evenness)
    if (!inherits(model, "try-error")) {

      converged <- tryCatch({
        model$sdr$pdHess
      }, error = function(e) FALSE)

      is_singular <- performance::check_singularity(model)

      # provide r2 value using r2_nakagawa
      r2_values <- r2_nakagawa(model)

      model_summary <- summary(model)
      intercept <- model_summary$coefficients$cond["(Intercept)", "Estimate"]

      # p values and beta coeffs
      spec_term <- paste0("scale(", spec, ")")
      p_value <- model_summary$coefficients$cond[spec_term, "Pr(>|z|)"]
      beta <- model_summary$coefficients$cond[spec_term, "Estimate"]

      # confidence intervals (Wald for speed/stability)
      conf_int <- confint(model, method = "Wald")
      beta_ci_lower <- conf_int[spec_term, 1]
      beta_ci_upper <- conf_int[spec_term, 2]

      # significance
      significance <- if_else(p_value < 0.05, 'yes', 'no')

      # store model results
      spectral_biodiversity_model_results <- bind_rows(spectral_biodiversity_model_results, data.frame(
        taxonomic_metric = tax,
        spectral_metric = spec,
        r2_marginal = r2_values$R2_marginal,
        r2_conditional = r2_values$R2_conditional,
        beta = beta,
        beta_ci_lower = beta_ci_lower,
        beta_ci_upper = beta_ci_upper,
        intercept = intercept,
        p_value = p_value,
        significance = significance,
        converged = converged,
        is_singular = is_singular
      ))
    }
  }
}

# CALCULATE CV FOR BAND COMBINATIONS + NDVI

# calculate CV based on band combinations

# red edge and nir
red_edge_nir_cv <- calculate_coefficient_of_variance(pixel_values, wavelengths = c('red_edge', 'nir'), min_points = 133609, rarefaction = T)
red_edge_nir_cv$bands <- 'red.edge_nir'

# green, red-edge, nir
green_re_nir_cv <- calculate_coefficient_of_variance(pixel_values, wavelengths = c('green', 'red_edge', 'nir'), min_points = 133609, rarefaction = T)
green_re_nir_cv$bands <- 'green_red.edge_nir'

# green, red, red-edge, nir
veg_bands_cv <- calculate_coefficient_of_variance(pixel_values, wavelengths = c('green', 'red', 'red_edge', 'nir'), min_points = 133609, rarefaction = T)
veg_bands_cv$bands <- 'green_red_red.edge_nir'

# CV of NDVI

# read in ndvi files - can be created from raster images available at 10.5281/zenodo.17089161
ndvi_files <- list.files('data/raster_images', pattern = '_combined_image_masked_NDVI.tif$', full.names = TRUE)

# extract ndvi values
ndvi_values <- extract_pixel_values_ndvi(ndvi_files, subplot_files, 'NDVI')

ndvi_cv <- calculate_coefficient_of_variance(ndvi_values, 'NDVI', min_points = 133609, rarefaction = T)
ndvi_cv$bands <- 'NDVI'

cv_values <- rbind(ndvi_cv, red_edge_nir_cv, veg_bands_cv, green_re_nir_cv)

# unique subplot id for each site
cv_values <- cv_values %>%
  mutate(subplot_id = case_when(
  site == 'NSABHC0009' ~ paste0('E',subplot_id),
  site == 'NSABHC0010' ~ paste0('G',subplot_id),
  site == 'NSABHC0011' ~ paste0('S',subplot_id),
  site == 'NSABHC0012' ~ paste0('C',subplot_id)
))

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
band_combinations <- c("NDVI", "red.edge_nir", "green_red.edge_nir", "green_red_red.edge_nir")

# empty df to store results
cv_biodiversity_model_results <- data.frame()


for (band in band_combinations) {
  df_filtered <- cv_values %>%
    filter(bands == band)

  for (tax in taxonomic_metrics) {

    # build model formula
    model_formula <- as.formula(
      paste0(tax, " ~ scale(CV) + (1 | site)")
    )

    # fit model
    model <- try(glmmTMB(model_formula, data = df_filtered))

    if (!inherits(model, "try-error")) {

      # get model summary and R²
      model_summary <- summary(model)
      r2_values <- r2_nakagawa(model)

      intercept <- model_summary$coefficients$cond["(Intercept)", "Estimate"]

      # extract spectral metric values
      spec_term <- paste0("scale(CV)")
      p_spec <- model_summary$coefficients$cond[spec_term, "Pr(>|z|)"]
      beta_spec <- model_summary$coefficients$cond[spec_term, "Estimate"]

      # confidence intervals
      conf_int <- confint(model, method = "Wald")
      beta_spec_ci_lower <- conf_int[spec_term, 1]
      beta_spec_ci_upper <- conf_int[spec_term, 2]

      # significance
      sig_spec <- if_else(p_spec < 0.05, "yes", "no")


      # store results
      cv_biodiversity_model_results <- bind_rows(cv_biodiversity_model_results, data.frame(
        bands = band,
        taxonomic_metric = tax,
        spectral_metric = 'CV',
        r2_marginal = r2_values$R2_marginal,
        r2_conditional = r2_values$R2_conditional,
        beta_spec = beta_spec,
        beta_spec_ci_lower = beta_spec_ci_lower,
        beta_spec_ci_upper = beta_spec_ci_upper,
        p_spec = p_spec,
        sig_spec = sig_spec,
        intercept = intercept
      ))

    }
  }
}

