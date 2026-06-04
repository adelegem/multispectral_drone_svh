# Figures for publication — manuscript PNG exporter.
#
# Produces 600-dpi PNGs at the dimensions the manuscript submission needs.
# Inline figure viewing happens in reports/report.html (Phase 5.3); this
# script exists alongside it solely to write the high-resolution PNGs.
#
# Prerequisite: legacy data_out/ inputs from continuous_metrics_analysis.R.
# A future cleanup would migrate readRDS() → tar_read() so this script
# pulls from the targets cache directly; not blocking publication.
#
# Figure 5: spectral vs taxonomic diversity scatter plots (masked_24_plot)
# Figure 6: CV beta coefficients across spectral band combinations (cv_bands_plot)
#
# All required packages are in renv.lock.

library(tidyverse)
library(glmmTMB)        # for predict.glmmTMB dispatch on cached _mixed fits
library(ggnewscale)
library(ggh4x)
library(patchwork)
source("funx.R")

metrics         <- readRDS("data_out/spectral_taxonomic_diversity.rds")
model_results   <- readRDS("data_out/spectral_biodiversity_model_results.rds")
cv_band_results <- readRDS("data_out/cv_biodiversity_model_results.rds")

# Cached fits from continuous_metrics_analysis.R (data_out/model_fits/), one per
# tax × spec combo. Mixed-effects models use the _mixed suffix; the three
# singular pielou_evenness models are refit as fixed-effect lm and use _fixed.
fit_files <- list.files("data_out/model_fits",
                        pattern = "^model_.*[.]rds$", full.names = TRUE)
fits <- setNames(
  lapply(fit_files, readRDS),
  fit_files |>
    basename() |>
    str_remove("^model_") |>
    str_remove("_(mixed|fixed)[.]rds$")
)

masked_24_plot <- make_figure_5(metrics, model_results, fits)
cv_bands_plot  <- make_figure_6(model_results, cv_band_results)

masked_24_plot
cv_bands_plot

if (!dir.exists("maps_graphs")) dir.create("maps_graphs")
ggsave("maps_graphs/masked_24_plot.png",          masked_24_plot, width = 9,  height = 8, dpi = 600)
ggsave("maps_graphs/cv_band_combinations_plot.png", cv_bands_plot, width = 12, height = 4, dpi = 600)

# Spectral-species figure
survey_data  <- read_csv("data/ausplots_march_24.csv", show_col_types = FALSE)
taxonomic_ss <- calculate_field_diversity(survey_data)$final_results %>%
  select(subplot_id, site, species_richness, shannon_diversity, simpson_diversity)
ss_spectral  <- readRDS("data_out/mean_spectral_species.rds") %>%
  select(subplot_id, site, spectral_richness, shannon_spectral, simpson_spectral)
ss_data <- ss_spectral %>%
  left_join(taxonomic_ss, by = c("subplot_id", "site")) %>%
  add_site_prefix()

ss_figure <- make_figure_ss(ss_data)
ss_figure
ggsave("maps_graphs/ss_figure.png", ss_figure, width = 12, height = 4, dpi = 600)
