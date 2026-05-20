# Figures for publication
#
# Prerequisite: run continuous_metrics_analysis.R first.
# Reads from data_out/ — all files written by that script.
#
# Figure 5: spectral vs taxonomic diversity scatter plots (masked_24_plot)
# Figure 6: CV beta coefficients across spectral band combinations (cv_bands_plot)
#
# Additional packages beyond those listed in the README:
#   install.packages(c("ggnewscale", "ggh4x"))

library(tidyverse)
library(glmmTMB)        # for predict.glmmTMB dispatch on cached _mixed fits
library(ggnewscale)
library(ggh4x)
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
ggsave("maps_graphs/cv_band_combinations_plot.png", cv_bands_plot, width = 10, height = 4, dpi = 600)
