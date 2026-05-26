# ROC-based threshold selection for NDVI and NIR masking
#
# Fits ROC curves (Youden index) to visually classified training points to find
# the optimal pixel-value threshold for separating vegetation from bare ground
# (NDVI) and non-shadow from shadow (NIR). These thresholds were then passed to
# create_masked_raster() (funx.R) to produce the masked rasters archived on
# Zenodo (10.5281/zenodo.17089161).
#
# Run from the project root:
#   Rscript preprocessing/02_roc_thresholds.R
#
# Inputs (tracked in this repo):
#   data/ndvi_reference_point_values.csv  — 100 pts/site, classes: veg | non-veg
#   data/nir_reference_point_values.csv   — 100 pts/site, classes: shadow | non-shadow
#
# Outputs:
#   data_out/ndvi_thresholds_2024.csv     — optimal threshold + diagnostics per site
#   data_out/nir_thresholds.csv
#   maps_graphs/roc_curves.png            — ROC curves for all sites x both metrics
#
# Expected thresholds (from original analysis):
#   NDVI: NSABHC0009 = 0.02198, NSABHC0010 = 0.02310, NSABHC0011 = 0.06775, NSABHC0012 = 0.04794
#   NIR:  NSABHC0009 = 0.03323, NSABHC0010 = 0.05510, NSABHC0011 = 0.04225, NSABHC0012 = 0.03722

library(pROC)
library(tidyverse)

if (!dir.exists("data_out"))    dir.create("data_out")
if (!dir.exists("maps_graphs")) dir.create("maps_graphs")

ndvi_pts <- read_csv("data/ndvi_reference_point_values.csv", show_col_types = FALSE)
nir_pts  <- read_csv("data/nir_reference_point_values.csv",  show_col_types = FALSE)


# =============================================================================
# Core function: ROC analysis for one metric across all sites
# =============================================================================
# positive_class: the class coded as 1 (target of detection)
# Returns a list with one element per site, each containing:
#   threshold, auc, sensitivity, specificity, n_pos, n_neg, roc object

run_roc_analysis <- function(df, value_col, class_col, positive_class,
                             site_col = "site") {
  lapply(unique(df[[site_col]]), function(s) {
    d <- df[df[[site_col]] == s, ]
    binary <- as.integer(d[[class_col]] == positive_class)

    roc_obj  <- roc(binary, d[[value_col]], quiet = TRUE, direction = "<")
    best     <- coords(roc_obj, "best",
                       ret = c("threshold", "sensitivity", "specificity"),
                       best.method = "youden")

    # coords() can return multiple rows if several thresholds tie — take the first
    best <- best[1, ]

    list(
      site        = s,
      roc         = roc_obj,
      threshold   = best$threshold,
      auc         = as.numeric(auc(roc_obj)),
      auc_ci      = as.numeric(ci.auc(roc_obj)),   # lower, AUC, upper (DeLong 95% CI)
      sensitivity = best$sensitivity,
      specificity = best$specificity,
      n_pos       = sum(binary == 1),
      n_neg       = sum(binary == 0)
    )
  }) |> setNames(unique(df[[site_col]]))
}


# =============================================================================
# Run analyses
# =============================================================================

ndvi_roc <- run_roc_analysis(ndvi_pts,
                             value_col      = "ndvi",
                             class_col      = "class",
                             positive_class = "veg")

nir_roc  <- run_roc_analysis(nir_pts,
                             value_col      = "nir",
                             class_col      = "class",
                             positive_class = "non-shadow")


# =============================================================================
# Extract summary tables
# =============================================================================

summarise_roc <- function(roc_list, value_label) {
  map_dfr(roc_list, function(x) {
    tibble(
      site          = x$site,
      metric        = value_label,
      threshold     = round(x$threshold,   6),
      auc           = round(x$auc,         4),
      auc_ci_lower  = round(x$auc_ci[1],   4),
      auc_ci_upper  = round(x$auc_ci[3],   4),
      sensitivity   = round(x$sensitivity, 3),
      specificity   = round(x$specificity, 3),
      n_veg_pos     = x$n_pos,
      n_nonveg_neg  = x$n_neg
    )
  })
}

ndvi_summary <- summarise_roc(ndvi_roc, "NDVI")
nir_summary  <- summarise_roc(nir_roc,  "NIR")

write_csv(ndvi_summary, "data_out/ndvi_thresholds_2024.csv")
write_csv(nir_summary,  "data_out/nir_thresholds.csv")

cat("=== NDVI thresholds (2024) ===\n")
print(ndvi_summary, width = 120)
cat("\n=== NIR thresholds ===\n")
print(nir_summary, width = 120)


# =============================================================================
# ROC curve plots
# =============================================================================

roc_plot_data <- function(roc_list, metric_label) {
  map_dfr(roc_list, function(x) {
    tibble(
      site        = x$site,
      metric      = metric_label,
      specificity = x$roc$specificities,
      sensitivity = x$roc$sensitivities
    )
  })
}

optimal_points <- function(summary_tbl) {
  summary_tbl %>%
    transmute(site, metric,
              x = 1 - specificity,
              y = sensitivity)
}

# One annotation block per facet, placed in the lower-right empty space.
# Each row is one site; columns are tab-aligned with sprintf padding.
panel_annotations <- bind_rows(ndvi_summary, nir_summary) %>%
  arrange(metric, site) %>%
  group_by(metric) %>%
  summarise(
    label = paste(
      sprintf("%-12s  t = %.4f   AUC = %.3f   Sens = %.3f   Spec = %.3f",
              site, threshold, auc, sensitivity, specificity),
      collapse = "\n"
    ),
    .groups = "drop"
  )

site_colours <- c(
  NSABHC0009 = "darkgreen",
  NSABHC0010 = "saddlebrown",
  NSABHC0011 = "navajowhite4",
  NSABHC0012 = "darkseagreen4"
)

plot_df <- bind_rows(
  roc_plot_data(ndvi_roc, "NDVI"),
  roc_plot_data(nir_roc,  "NIR")
)

opt_df <- bind_rows(
  optimal_points(ndvi_summary),
  optimal_points(nir_summary)
)

roc_plot <- ggplot(plot_df,
                   aes(x = 1 - specificity, y = sensitivity, colour = site)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = opt_df, aes(x = x, y = y, colour = site),
             size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
  geom_label(data = panel_annotations,
             aes(x = 0.97, y = 0.08, label = label),
             hjust = 1, vjust = 0, size = 2.6,
             inherit.aes = FALSE,
             label.size = 0.3, fill = "grey97", colour = "grey20",
             family = "mono") +
  facet_wrap(~ metric, nrow = 1) +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(x = "False positive rate (1 – Specificity)",
       y = "True positive rate (Sensitivity)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 13),
    legend.position   = "bottom",
    plot.background   = element_rect(fill = "white", colour = "white")
  )

roc_plot

ggsave("maps_graphs/roc_curves.png",      roc_plot, width = 12, height = 6, dpi = 300)
ggsave("preprocessing/roc_curves.png", roc_plot, width = 12, height = 6, dpi = 120)
