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
library(glmmTMB)
library(ggnewscale)
library(ggh4x)

metrics      <- readRDS("data_out/spectral_taxonomic_diversity.rds")
model_results <- readRDS("data_out/spectral_biodiversity_model_results.rds")

taxonomic_labels <- c(
  species_richness = "Species\nRichness",
  exp_shannon      = "Exponential\nShannon's",
  inv_simpson      = "Inverse\nSimpson's",
  pielou_evenness  = "Pielou's\nEvenness"
)

spectral_labels <- c(
  CV      = "CV",
  SV      = "SV",
  log.CHV = "log(CHV)"
)

df_long <- metrics %>%
  pivot_longer(
    cols      = c(species_richness, exp_shannon, inv_simpson, pielou_evenness),
    names_to  = "taxonomic_metric",
    values_to = "taxonomic_value"
  ) %>%
  pivot_longer(
    cols      = c(CV, SV, log.CHV),
    names_to  = "spectral_metric",
    values_to = "spectral_value"
  ) %>%
  mutate(taxonomic_metric = factor(taxonomic_metric, levels = names(taxonomic_labels)))

# Re-fit models to generate predicted lines (fast — no rarefaction).
# Significance comes from the saved model_results; only significant
# relationships get a line drawn.
df_long$predicted_values <- NA_real_

for (tax in levels(df_long$taxonomic_metric)) {
  for (spec in unique(df_long$spectral_metric)) {

    model <- try(glmmTMB(
      formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
      data    = metrics
    ), silent = TRUE)

    if (!inherits(model, "try-error")) {
      preds <- predict(model, newdata = metrics, type = "response")
      rows  <- df_long$taxonomic_metric == tax & df_long$spectral_metric == spec
      df_long$predicted_values[rows] <- preds
    }
  }
}

df_long <- df_long %>%
  left_join(
    model_results %>% dplyr::select(taxonomic_metric, spectral_metric, significance),
    by = c("taxonomic_metric", "spectral_metric")
  ) %>%
  mutate(taxonomic_metric = factor(taxonomic_metric, levels = c(
    "species_richness", "exp_shannon", "inv_simpson", "pielou_evenness"
  )))

masked_24_plot <- df_long %>%
  ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
  geom_point(alpha = 0.7) +
  facet_grid(
    taxonomic_metric ~ spectral_metric,
    scales   = "free",
    labeller = labeller(taxonomic_metric = taxonomic_labels,
                        spectral_metric  = spectral_labels),
    switch   = "both"
  ) +
  labs(title = "", color = "Site") +
  theme_minimal() +
  scale_color_manual(values = c("darkgreen", "saddlebrown", "navajowhite2", "darkseagreen")) +
  new_scale_colour() +
  geom_line(
    data = df_long %>% filter(significance == "yes"),
    aes(x = spectral_value, y = predicted_values, color = site),
    linetype = "solid", linewidth = 1, show.legend = FALSE
  ) +
  scale_color_manual(values = rep("black", 4)) +
  facetted_pos_scales(
    y = list(
      taxonomic_metric == "species_richness" ~ scale_y_continuous(breaks = c(5, 10, 15),   limits = c(0, 20)),
      taxonomic_metric == "exp_shannon"      ~ scale_y_continuous(breaks = c(5, 10, 15),   limits = c(0.7, 16)),
      taxonomic_metric == "inv_simpson"      ~ scale_y_continuous(breaks = c(4, 8, 12),    limits = c(0, 12.5)),
      taxonomic_metric == "pielou_evenness"  ~ scale_y_continuous(breaks = c(0.8, 0.9, 1), limits = c(0.65, 1))
    )
  ) +
  facetted_pos_scales(
    x = list(
      spectral_metric == "CV"      ~ scale_x_continuous(breaks = c(0.15, 0.25, 0.35),         limits = c(0.1, 0.4)),
      spectral_metric == "log.CHV" ~ scale_x_continuous(breaks = c(-20.0, -18.0, -16.0),       limits = c(-20, -15.5)),
      spectral_metric == "SV"      ~ scale_x_continuous(breaks = c(0.0005, 0.0015, 0.0025),    limits = c(0, 0.0027))
    )
  ) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x        = element_text(angle = 90, hjust = 1, size = 12),
    axis.title.x       = element_text(size = 16),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 14),
    strip.text.x       = element_text(size = 14),
    strip.text.y       = element_text(size = 14),
    axis.title.x.bottom = element_blank(),
    axis.title.y       = element_blank(),
    strip.placement    = "outside",
    legend.position    = "bottom",
    legend.text        = element_text(size = 13),
    legend.title       = element_text(size = 13),
    plot.background    = element_rect(fill = "white", color = "white")
  )

masked_24_plot

if (!dir.exists("maps_graphs")) dir.create("maps_graphs")
ggsave("maps_graphs/masked_24_plot.png", masked_24_plot, width = 9, height = 8, dpi = 600)


# =============================================================================
# Figure 6: CV beta coefficients across spectral band combinations
# =============================================================================
# Reads data_out/cv_biodiversity_model_results.rds  (band-combo models)
# and reuses model_results (already loaded above) for the all-bands CV baseline.

cv_band_results <- readRDS("data_out/cv_biodiversity_model_results.rds")

# Pull the all-bands CV rows from the main model results and align column names
all_bands_cv <- model_results %>%
  filter(spectral_metric == "CV") %>%
  transmute(
    bands              = "all_bands",
    taxonomic_metric,
    beta_spec          = beta,
    beta_spec_ci_lower = beta_ci_lower,
    beta_spec_ci_upper = beta_ci_upper,
    sig_spec           = significance
  )

all_cv_results <- bind_rows(all_bands_cv, cv_band_results) %>%
  filter(taxonomic_metric != "pielou_evenness") %>%
  mutate(
    bands = factor(bands, levels = c(
      "red.edge_nir",
      "green_red.edge_nir",
      "green_red_red.edge_nir",
      "all_bands"
    )),
    taxonomic_metric = factor(taxonomic_metric, levels = c(
      "species_richness", "exp_shannon", "inv_simpson"
    ))
  )

band_labels <- c(
  red.edge_nir           = "Red edge and NIR",
  green_red.edge_nir     = "Green, red edge, NIR",
  green_red_red.edge_nir = "Vegetation bands",
  all_bands              = "All Bands"
)

band_colours <- c(
  all_bands              = "black",
  green_red_red.edge_nir = "#B31441",
  green_red.edge_nir     = "#D95F02",
  red.edge_nir           = "#E6AB02"
)

taxonomic_labels_fig6 <- c(
  species_richness = "Species Richness",
  exp_shannon      = "Exponential Shannon",
  inv_simpson      = "Inverse Simpson"
)

cv_bands_plot <- all_cv_results %>%
  ggplot(aes(x = beta_spec, y = bands, color = bands)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = beta_spec_ci_lower, xmax = beta_spec_ci_upper),
    height = 0.1, linewidth = 0.8
  ) +
  geom_point(size = 3.5, stroke = 1) +
  scale_color_manual(values = band_colours) +
  scale_y_discrete(labels = band_labels) +
  facet_wrap(
    ~ taxonomic_metric, nrow = 1,
    labeller = labeller(taxonomic_metric = taxonomic_labels_fig6)
  ) +
  labs(x = "Beta co-efficient", y = NULL) +
  guides(color = "none") +
  theme_minimal() +
  theme(
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(2, "lines"),
    text              = element_text(size = 20),
    strip.text        = element_text(size = 18, face = "bold"),
    axis.text.y       = element_text(size = 16),
    axis.text.x       = element_text(size = 14),
    axis.title.x      = element_text(size = 18),
    axis.title.y      = element_blank(),
    plot.background   = element_rect(fill = "white", color = "white")
  )

cv_bands_plot

ggsave("maps_graphs/cv_band_combinations_plot.png", cv_bands_plot, width = 10, height = 4, dpi = 600)
