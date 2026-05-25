

# DOWNLOAD RASTERS FROM ZENODO -----------------------------------------------
# Downloads the 4 masked multispectral rasters from Zenodo record 10.5281/zenodo.17089161
# into `dest_dir`. Files already present (matching expected size within 1%) are skipped.
# Total download is ~4.6 GB.
download_zenodo_rasters <- function(dest_dir = "data/raster_images") {
  files <- list(
    list(name = "NSABHC0009_masked.tif",
         url  = "https://zenodo.org/records/17089161/files/NSABHC0009_masked.tif?download=1",
         size = 1.6e9),
    list(name = "NSABHC0010_masked.tif",
         url  = "https://zenodo.org/records/17089161/files/NSABHC0010_masked.tif?download=1",
         size = 553.9e6),
    list(name = "NSABHC0011_masked.tif",
         url  = "https://zenodo.org/records/17089161/files/NSABHC0011_masked.tif?download=1",
         size = 956.5e6),
    list(name = "NSABHC0012_masked.tif",
         url  = "https://zenodo.org/records/17089161/files/NSABHC0012_masked.tif?download=1",
         size = 1.5e9)
  )

  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (f in files) {
    dest <- file.path(dest_dir, f$name)
    if (file.exists(dest) && file.info(dest)$size > 0.95 * f$size) {
      message("Already present, skipping: ", dest)
      next
    }
    message("Downloading ", f$name, " (~", round(f$size / 1e6), " MB) ...")
    utils::download.file(f$url, destfile = dest, mode = "wb", quiet = FALSE)
  }

  invisible(list.files(dest_dir, pattern = "_masked\\.tif$", full.names = TRUE))
}

# COMBINE TIFS FUNCTION (FOR CREATING MULTI BAND IMAGE) -----------------------------------------------------------------------

create_multiband_image <- function(mosaics_dir, desired_band_order){
# folder list | recursive = won't pick folders within folders
folders <- list.dirs(mosaics_dir, full.names = FALSE, recursive = FALSE)
# need to make this part an argument e.g. option to exclude certain folders or include certain folders
folders <- folders[folders != "point_clouds"]

## NOTE: spectral band image tif file names must be named after their band (e.g., blue, nir, etc),
#  otherwise change 'desired_band_order' to match file names
#  should be combined in wavelength order, esp for biodivmapR processes (i.e. as above)

# loop thru each folder
for (folder in folders) {
  # create path
  folder_path <- file.path(mosaics_dir, folder)

  # list of tif files | \\. represents . (dots need to be escaped w \, \ need to be escaped with  \). $means at end of file name/string
  tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

  # load as raster
  rasters <- lapply(tif_files, terra::rast)

  # extract band names from file names
  band_names <- tools::file_path_sans_ext(basename(tif_files))

  # stack rasters and assign band names
  combined_image <- terra::rast(rasters)
  names(combined_image) <- band_names

  # reorder the bands based on the desired band order
  combined_image <- combined_image[[match(desired_band_order, band_names)]]

  #create output directory folder if it doesn't exist
  output_dir <- file.path("data_out/combined_rasters", substr(basename(mosaics_dir), 1, 4))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # create output file as .tif and as .envi
  output_filename <- file.path(output_dir, paste0(folder, "_combined_image"))
  # write .tif file
  terra::writeRaster(combined_image, filename = paste0(output_filename, '.tif'),
                     filetype = "GTiff", gdal = c("INTERLEAVE=BAND"), overwrite = TRUE)
  plot(combined_image)
}
}

# FIND OPTIMUM THRESHOLDS FUNCTION
# class = class column for classifications (e.g. veg, ground etc) (col name)
# value = ndvi, nir - what threshold are you seeking (col name)
# site = plot reference (col name)

find_optimum_thresholds <- function(df, class, value, site, class_value) {
  # empty df
  threshold_df <- data.frame(site = character(), threshold = numeric())

  # iterate over sites
  for (site_name in unique(df[[site]])) {

    # filter current location
    site_data <- subset(df, df[[site]] == site_name)

    # binary outcome variable for veg and non-veg
    site_data$binary_class <- ifelse(site_data[[class]] == class_value, 1, 0)

    # ROC curve
    roc_result <- pROC::roc(site_data$binary_class, site_data[[value]])

    # find optimum threshold
    best_threshold <- coords(roc_result, 'best')$threshold

    # append to the result data frame
    threshold_df <- rbind(threshold_df, data.frame(site = site_name, threshold = best_threshold))
  }

  return(threshold_df)
}

# CREATE_MASKED_RASTER FUNCTION
#input can be directory with a number of files, a single file, or string of files.
#ndvi and nir thresholds can be provided as a df, if there are diff optimum values per site
# or as a single value for all sites
# this function assumes that layers are stacked in WAVELENGTH ORDER
#think about how you can make this more general for users - e.g. it requires the plot id to
# be in the file name currently - think about usability
create_masked_raster <- function(input, output_dir,
                                 band_names,
                                 NDVI_Thresh = 0.2, NIR_Thresh = 0.2, Red_Thresh= 0.7,
                                 NDVI_Thresh_df = NULL, NIR_Thresh_df = NULL, Red_Thresh_df = NULL,
                                 red_band_index = 3, NIR_band_index = 5) {

  if (dir.exists(input)) {
    # list all ENVI or TIF files in the directory
    files <- list.files(file.path(input), pattern = '\\.(envi|tif)$', full.names = TRUE)
  } else if (file.exists(input) || is.character(input)) {
    # single file input or string of files
    files <- input
  } else {
    stop("Invalid input provided.")
  }

  print(paste("Files found:", files))

  if (length(files) == 0) {
    stop("No files found.")
  }

  for (file in files) {
    # extract the site identifier from file name
    file_id <- strsplit(basename(file), "_")[[1]][1]

    # check if NDVI_Thresh_df is provided and extract the relevant threshold values
    if (!is.null(NDVI_Thresh_df)) {
      if (file_id %in% NDVI_Thresh_df[[1]]) {
        NDVI_Thresh <- NDVI_Thresh_df[[2]][NDVI_Thresh_df[[1]] == file_id]
      } else {
        stop(paste("No NDVI threshold values found for file", file_id))
      }
    }

    # check if NIR_Thresh_df is provided and extract the relevant threshold values
    if (!is.null(NIR_Thresh_df)) {
      if (file_id %in% NIR_Thresh_df[[1]]) {
        NIR_Thresh <- NIR_Thresh_df[[2]][NIR_Thresh_df[[1]] == file_id]
      } else {
        stop(paste("No NIR threshold values found for site", file_id))
      }
    }

    # check if NIR_Thresh_df is provided and extract the relevant threshold values
    if (!is.null(Red_Thresh_df)) {
      if (file_id %in% Red_Thresh_df[[1]]) {
        Red_Thresh <- Red_Thresh_df[[2]][Red_Thresh_df[[1]] == file_id]
      } else {
        stop(paste("No NIR threshold values found for site", file_id))
      }
    }

    # Read the multi-band raster
    raster_data <- terra::rast(file)

    # Identify the bands for Red and NIR
    red <- raster_data[[red_band_index]]
    nir <- raster_data[[NIR_band_index]]

    # Calculate NDVI
    ndvi <- (nir - red) / (nir + red)

    # Create a mask based on NDVI and NIR thresholds
    mask <- (ndvi < NDVI_Thresh) | (nir < NIR_Thresh) | (red > Red_Thresh)

    # Apply the mask to the raster data
    raster_data_masked <- terra::mask(raster_data, mask, maskvalue = TRUE, updatevalue = NA)

    # Save the masked raster
    masked_filename <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file)), '_masked.tif'))
    terra::writeRaster(raster_data_masked, filename = masked_filename, filetype = "GTiff", overwrite = TRUE)

    print(paste("Masked raster saved to:", masked_filename))
  }
}




## EXTRACT PIXEL VALUES FUNCTION

# Wraps saltbush::extract_pixel_values to preserve the analysis-specific
# "row_col" subplot_id and `site` column name. Fishnets are read in row-major
# order, so aoi_id 1..25 maps 1:1 to subplot_id 1_1, 1_2, ..., 5_5.
extract_pixel_values <- function(raster_files, subplot_files, wavelength_names) {
  subplot_ids <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep = "_")))

  saltbush::extract_pixel_values(
    raster_files     = raster_files,
    aoi_files        = subplot_files,
    wavelength_names = wavelength_names
  ) %>%
    mutate(subplot_id = subplot_ids[aoi_id]) %>%
    dplyr::select(-aoi_id) %>%
    rename(site = site_name)
}


## EXTRACT PIXEL VALUES NDVI

extract_pixel_values_ndvi <- function(raster_files, subplot_files, wavelength_names) {

  all_pixel_values_list <- list()

  for (raster_file in raster_files) {

    # identify the string that represents the site name
    site <- str_extract(basename(raster_file), "^[^_]+")

    # choose matching subplot shapefile
    subplot_file <- subplot_files[grep(paste0('^', site), basename(subplot_files))]

    if (length(subplot_file) == 0) next  # skip if no matching subplot

    # read subplot shapefile
    subplots <- read_sf(subplot_file) %>%
      dplyr::select('geometry')

    subplots$subplot_id <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep = "_")))

    # read raster
    raster_data <- rast(raster_file)
    names(raster_data) <- wavelength_names

    pixel_values_list <- list()

    for (i in 1:nrow(subplots)) {
      subplot <- subplots[i, ]
      subplot_id <- subplot$subplot_id
      subplot_sp <- as(subplot, "Spatial")

      subplot_vect <- vect(subplot)  # convert sf to SpatVector

      cropped_raster <- crop(raster_data, subplot_vect)
      masked_raster <- mask(cropped_raster, subplot_vect)

      pixel_values <- as.data.frame(masked_raster)
      colnames(pixel_values) <- wavelength_names
      pixel_values$subplot_id <- subplot_id
      pixel_values_list[[i]] <- pixel_values
    }

    # Combine and add metadata
    all_pixel_values <- bind_rows(pixel_values_list) %>%
      na.omit()

    # add to overall list with all raster data pixel values
    all_pixel_values_list[[site]] <- all_pixel_values
  }

  combined_values <- bind_rows(all_pixel_values_list, .id = 'site')

  return(combined_values)
}

## CALCULATE SPECTRAL METRICS FUNCTIONS

#rarefraction cv function
calculate_cv <- function(pixel_values_df,
                         wavelengths,
                         rarefaction = TRUE,
                         min_points = NULL,
                         n = 999,
                         seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Convert the dataframe to a data.table for efficiency
  setDT(pixel_values_df)

  if (rarefaction) {
    # Initialize a list to store CV values for each replication
    cv_list <- vector("list", n)

    # Resample and calculate CV for each iteration
    for (i in seq_len(n)) {
      # Sample to the minimum number of points per subplot
      sampled_df <- pixel_values_df[, .SD[sample(.N, min_points)], by = subplot_id, .SDcols = wavelengths]

      # Calculate CV for each wavelength within each subplot
      cv_data <- sampled_df[, lapply(.SD, function(x) sd(x) / abs(mean(x, na.rm = TRUE))), by = subplot_id]

      # Sum across wavelengths and normalize by the number of bands (ignoring NAs)
      cv_data[, CV := rowSums(.SD, na.rm = TRUE) / (length(wavelengths) - rowSums(is.na(.SD))), .SDcols = wavelengths]

      # Store CV values along with subplots
      cv_list[[i]] <- cv_data[, .(subplot_id = subplot_id, CV), by = subplot_id]  # Ensure subplot_id is included
    }

    # Collapse the list of CV data tables into a single data table and calculate the average CV for each subplot
    cv <- rbindlist(cv_list)[, .(CV = mean(CV, na.rm = TRUE)), by = subplot_id]

  } else {
    # If rarefaction is FALSE, directly calculate the CV without resampling
    cv <- pixel_values_df[, lapply(.SD, function(x) sd(x) / abs(mean(x, na.rm = TRUE))), by = subplot_id, .SDcols = wavelengths]

    # Sum across wavelengths and normalize by the number of bands (ignoring NAs)
    cv[, CV := rowSums(.SD, na.rm = TRUE) / (length(wavelengths) - rowSums(is.na(.SD))), .SDcols = wavelengths]

    # Ensure subplot_id is included in the output
    cv <- cv[, .(subplot_id = subplot_id, CV)]
  }

  return(cv)
}



# Wraps saltbush::calculate_spectral_metrics to preserve the analysis-specific
# column names (site/subplot_id/CHV_nopca) used by the rest of the pipeline.
calculate_spectral_metrics <- function(pixel_values_df,
                                       wavelengths,
                                       min_points,
                                       n = 999,
                                       rarefaction = TRUE,
                                       seed = NULL) {
  renamed <- pixel_values_df %>%
    rename(site_name = site, aoi_id = subplot_id)

  out <- saltbush::calculate_spectral_metrics(
    renamed,
    wavelengths = wavelengths,
    rarefaction = rarefaction,
    min_points  = min_points,
    n           = n,
    seed        = seed
  )

  out %>%
    rename(subplot_id = aoi_id, CHV_nopca = CHV) %>%
    dplyr::select(site, subplot_id, CV, SV, CHV_nopca)
}

# CALCULATE CV ONLY
## FUNCTION FOR CALCULATING ALL METRICS
calculate_coefficient_of_variance <- function(pixel_values_df,
                                       wavelengths,
                                       min_points,
                                       n = 999,
                                       rarefaction = TRUE,
                                       seed = NULL) {  # Add rarefaction here
  results <- list()

  for (site in unique(pixel_values_df$site)) {
    site_pixel_values <- pixel_values_df %>% filter(site == !!site)

    # Calculate CV, pass rarefaction where needed
    cv <- calculate_cv(site_pixel_values, wavelengths = wavelengths, rarefaction = rarefaction, n = n, min_points = min_points, seed = seed)
    results[[site]] <- list(CV = cv)
  }

  combined_cv <- bind_rows(lapply(results, function(x) x$CV), .id = 'site')


  # create a data frame for cv only
  combined_metrics <- combined_cv


  return(combined_metrics)
}



# CALCULATE_FIELD_DIVERSITY -----------------------------------------------

# Reconstruct subplot_id (5x5 grid) from AusPlots transect/point columns.
# Returns survey_data with added subplot_id (+ helper columns), filtered to drop
# point_number2 == 100 which would otherwise straddle subplot boundaries.
bin_survey_subplots <- function(survey_data) {
  survey_data %>%
    mutate(
      transect_direction = gsub('[[:digit:]]+', '', transect),
      transect_number = as.numeric(gsub(".*?([0-9]+).*", "\\1", transect)),
      point_number2 = case_when(
        transect_direction == "E-W" ~ 100 - point_number,
        transect_direction == "N-S" ~ 100 - point_number,
        TRUE ~ point_number
      ),
      transect_direction2 = case_when(
        transect_direction %in% c("W-E", "E-W") ~ "W-E",
        transect_direction %in% c("N-S", "S-N") ~ "S-N"
      ),
      X_plot = case_when(
        transect_direction2 == "W-E" ~ point_number2,
        transect_direction2 == "S-N" & transect_number == 1 ~ 10,
        transect_direction2 == "S-N" & transect_number == 2 ~ 30,
        transect_direction2 == "S-N" & transect_number == 3 ~ 50,
        transect_direction2 == "S-N" & transect_number == 4 ~ 70,
        transect_direction2 == "S-N" & transect_number == 5 ~ 90
      ),
      Y_plot = case_when(
        transect_direction2 == "S-N" ~ point_number2,
        transect_direction2 == "W-E" & transect_number == 1 ~ 10,
        transect_direction2 == "W-E" & transect_number == 2 ~ 30,
        transect_direction2 == "W-E" & transect_number == 3 ~ 50,
        transect_direction2 == "W-E" & transect_number == 4 ~ 70,
        transect_direction2 == "W-E" & transect_number == 5 ~ 90
      )
    ) %>%
    filter(point_number2 != 100) %>%
    mutate(
      subplot_row = pmin(ceiling((Y_plot + 1) / 20), 5),
      subplot_col = pmin(ceiling((X_plot + 1) / 20), 5),
      subplot_id  = paste(subplot_row, subplot_col, sep = "_")
    )
}

# Wraps saltbush::calculate_field_diversity to compute per-subplot diversity.
# Subplot binning (transect/point -> 5x5 grid) is AusPlots-specific and stays
# local; the per-group diversity math is delegated to saltbush.
calculate_field_diversity <- function(survey_data) {
  binned <- bin_survey_subplots(survey_data)

  out <- saltbush::calculate_field_diversity(
    binned,
    group_by_cols = c("site_location_name", "subplot_id")
  )

  final_results <- out$taxonomic_diversity %>%
    rename(site = site_location_name) %>%
    select(subplot_id, species_richness, shannon_diversity,
           simpson_diversity, pielou_evenness, exp_shannon,
           inv_simpson, site) %>%
    arrange(site, subplot_id)

  list(
    final_results = final_results,
    community_matrices = out$community_matrices
  )
}


# SITE-PREFIX HELPER ----------------------------------------------------------

# Prepend a single-letter code (E/G/S/C) derived from the AusPlots site name to
# `subplot_id`, so subplot ids remain unique after binding the four sites
# (each site has its own 5×5 grid with overlapping "row_col" ids).
add_site_prefix <- function(data) {
  data %>%
    mutate(subplot_id = case_when(
      site == "NSABHC0009" ~ paste0("E", subplot_id),
      site == "NSABHC0010" ~ paste0("G", subplot_id),
      site == "NSABHC0011" ~ paste0("S", subplot_id),
      site == "NSABHC0012" ~ paste0("C", subplot_id)
    ))
}


# MIXED-MODEL FITTING ---------------------------------------------------------

# Fit one mixed model: <tax_metric> ~ scale(<spec_metric>) + (1 | site).
# If the mixed fit is singular (site random-intercept variance ≈ 0, which
# pielou_evenness consistently triggers across our four sites), automatically
# refit as a fixed-effect lm without (1 | site). Returns the fitted model
# (class glmmTMB or lm). Caller doesn't need to know which kind it got —
# summarise_spectral_biodiversity_model() handles both.
fit_spectral_biodiversity_model <- function(data, tax_metric, spec_metric,
                                             scale_predictor = TRUE) {
  pred_expr <- if (scale_predictor) paste0("scale(", spec_metric, ")") else spec_metric
  mixed_formula <- as.formula(paste(tax_metric, "~", pred_expr, "+ (1 | site)"))
  mixed <- glmmTMB::glmmTMB(mixed_formula, data = data)

  if (isTRUE(performance::check_singularity(mixed))) {
    fixed_formula <- as.formula(paste(tax_metric, "~", pred_expr))
    return(stats::lm(fixed_formula, data = data))
  }

  mixed
}

# One-row summary tibble from a fitted spectral-biodiversity model. Accepts
# either glmmTMB (mixed) or lm (singular-refit fallback). `model_kind`
# records which was used; `is_singular` is TRUE only when the mixed fit was
# singular and the refit kicked in.
summarise_spectral_biodiversity_model <- function(model, tax_metric, spec_metric,
                                                  scale_predictor = TRUE) {
  spec_term <- if (scale_predictor) paste0("scale(", spec_metric, ")") else spec_metric

  if (inherits(model, "glmmTMB")) {
    coefs     <- summary(model)$coefficients$cond
    conf_int  <- confint(model, method = "Wald")
    r2_values <- performance::r2_nakagawa(model)
    p_value   <- coefs[spec_term, "Pr(>|z|)"]

    tibble(
      taxonomic_metric = tax_metric,
      spectral_metric  = spec_metric,
      model_kind       = "mixed",
      r2_marginal      = unname(r2_values$R2_marginal),
      r2_conditional   = unname(r2_values$R2_conditional),
      beta             = coefs[spec_term, "Estimate"],
      beta_ci_lower    = conf_int[spec_term, 1],
      beta_ci_upper    = conf_int[spec_term, 2],
      intercept        = coefs["(Intercept)", "Estimate"],
      p_value          = p_value,
      significance     = if_else(p_value < 0.05, "yes", "no"),
      converged        = tryCatch(model$sdr$pdHess, error = function(e) FALSE),
      is_singular      = FALSE
    )
  } else if (inherits(model, "lm")) {
    coefs    <- summary(model)$coefficients
    conf_int <- confint(model)
    r2       <- summary(model)$r.squared
    p_value  <- coefs[spec_term, "Pr(>|t|)"]

    tibble(
      taxonomic_metric = tax_metric,
      spectral_metric  = spec_metric,
      model_kind       = "fixed",
      r2_marginal      = r2,
      r2_conditional   = r2,                # no random effect; marginal == conditional
      beta             = coefs[spec_term, "Estimate"],
      beta_ci_lower    = conf_int[spec_term, 1],
      beta_ci_upper    = conf_int[spec_term, 2],
      intercept        = coefs["(Intercept)", "Estimate"],
      p_value          = p_value,
      significance     = if_else(p_value < 0.05, "yes", "no"),
      converged        = TRUE,
      is_singular      = TRUE               # marks "mixed was singular → refit"
    )
  } else {
    stop("Unexpected model class: ", paste(class(model), collapse = ", "))
  }
}

# Fit one CV-vs-taxonomic model on the band-combination subset of cv_values:
# <tax_metric> ~ scale(CV) + (1 | site), restricted to rows where bands == band_combo.
# Same singular-refit fallback as fit_spectral_biodiversity_model().
fit_cv_band_model <- function(data, tax_metric, band_combo) {
  df_filtered   <- data %>% filter(bands == band_combo)
  mixed_formula <- as.formula(paste0(tax_metric, " ~ scale(CV) + (1 | site)"))
  mixed <- glmmTMB::glmmTMB(mixed_formula, data = df_filtered)

  if (isTRUE(performance::check_singularity(mixed))) {
    fixed_formula <- as.formula(paste0(tax_metric, " ~ scale(CV)"))
    return(stats::lm(fixed_formula, data = df_filtered))
  }

  mixed
}

# One-row summary tibble from a fitted CV-band model. Same handling as
# summarise_spectral_biodiversity_model() but with the cv-band column naming
# (bands, beta_spec, p_spec, sig_spec) the manuscript figure uses.
summarise_cv_band_model <- function(model, tax_metric, band_combo) {
  spec_term <- "scale(CV)"

  if (inherits(model, "glmmTMB")) {
    coefs     <- summary(model)$coefficients$cond
    conf_int  <- confint(model, method = "Wald")
    r2_values <- performance::r2_nakagawa(model)
    p_spec    <- coefs[spec_term, "Pr(>|z|)"]

    tibble(
      bands              = band_combo,
      taxonomic_metric   = tax_metric,
      spectral_metric    = "CV",
      model_kind         = "mixed",
      r2_marginal        = unname(r2_values$R2_marginal),
      r2_conditional     = unname(r2_values$R2_conditional),
      beta_spec          = coefs[spec_term, "Estimate"],
      beta_spec_ci_lower = conf_int[spec_term, 1],
      beta_spec_ci_upper = conf_int[spec_term, 2],
      p_spec             = p_spec,
      sig_spec           = if_else(p_spec < 0.05, "yes", "no"),
      intercept          = coefs["(Intercept)", "Estimate"]
    )
  } else if (inherits(model, "lm")) {
    coefs    <- summary(model)$coefficients
    conf_int <- confint(model)
    r2       <- summary(model)$r.squared
    p_spec   <- coefs[spec_term, "Pr(>|t|)"]

    tibble(
      bands              = band_combo,
      taxonomic_metric   = tax_metric,
      spectral_metric    = "CV",
      model_kind         = "fixed",
      r2_marginal        = r2,
      r2_conditional     = r2,
      beta_spec          = coefs[spec_term, "Estimate"],
      beta_spec_ci_lower = conf_int[spec_term, 1],
      beta_spec_ci_upper = conf_int[spec_term, 2],
      p_spec             = p_spec,
      sig_spec           = if_else(p_spec < 0.05, "yes", "no"),
      intercept          = coefs["(Intercept)", "Estimate"]
    )
  } else {
    stop("Unexpected model class: ", paste(class(model), collapse = ", "))
  }
}


# PUBLICATION FIGURES ---------------------------------------------------------

# Figure 5 — per spec × tax scatter (4×3 facet) with model-prediction lines drawn
# only where `model_results$significance == "yes"`. Predictions come from `fits`
# (cached fits from data_out/model_fits/, named "<tax>_<spec>"), so this
# function does no model fitting itself. Requires ggnewscale + ggh4x loaded.
make_figure_5 <- function(metrics, model_results, fits) {
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
    mutate(
      taxonomic_metric = factor(taxonomic_metric, levels = names(taxonomic_labels)),
      predicted_values = NA_real_
    )

  for (tax in levels(df_long$taxonomic_metric)) {
    for (spec in unique(df_long$spectral_metric)) {
      fit <- fits[[paste0(tax, "_", spec)]]
      if (is.null(fit)) next
      preds <- predict(fit, newdata = metrics, type = "response")
      rows  <- df_long$taxonomic_metric == tax & df_long$spectral_metric == spec
      df_long$predicted_values[rows] <- preds
    }
  }

  df_long <- df_long %>%
    left_join(
      model_results %>% dplyr::select(taxonomic_metric, spectral_metric, significance),
      by = c("taxonomic_metric", "spectral_metric")
    )

  df_long %>%
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
    ggnewscale::new_scale_colour() +
    geom_line(
      data = df_long %>% filter(significance == "yes"),
      aes(x = spectral_value, y = predicted_values, color = site),
      linetype = "solid", linewidth = 1, show.legend = FALSE
    ) +
    scale_color_manual(values = rep("black", 4)) +
    ggh4x::facetted_pos_scales(
      y = list(
        taxonomic_metric == "species_richness" ~ scale_y_continuous(breaks = c(5, 10, 15),   limits = c(0, 20)),
        taxonomic_metric == "exp_shannon"      ~ scale_y_continuous(breaks = c(5, 10, 15),   limits = c(0.7, 16)),
        taxonomic_metric == "inv_simpson"      ~ scale_y_continuous(breaks = c(4, 8, 12),    limits = c(0, 12.5)),
        taxonomic_metric == "pielou_evenness"  ~ scale_y_continuous(breaks = c(0.8, 0.9, 1), limits = c(0.65, 1))
      )
    ) +
    ggh4x::facetted_pos_scales(
      x = list(
        spectral_metric == "CV"      ~ scale_x_continuous(breaks = c(0.15, 0.25, 0.35),      limits = c(0.1, 0.4)),
        spectral_metric == "log.CHV" ~ scale_x_continuous(breaks = c(-20.0, -18.0, -16.0),   limits = c(-20, -15.5)),
        spectral_metric == "SV"      ~ scale_x_continuous(breaks = c(0.0005, 0.0015, 0.0025), limits = c(0, 0.0027))
      )
    ) +
    theme(
      plot.title          = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x         = element_text(angle = 90, hjust = 1, size = 12),
      axis.title.x        = element_text(size = 16),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      axis.text.y         = element_text(size = 14),
      strip.text.x        = element_text(size = 14),
      strip.text.y        = element_text(size = 14),
      axis.title.x.bottom = element_blank(),
      axis.title.y        = element_blank(),
      strip.placement     = "outside",
      legend.position     = "bottom",
      legend.text         = element_text(size = 13),
      legend.title        = element_text(size = 13),
      plot.background     = element_rect(fill = "white", color = "white")
    )
}

# Figure 6 — CV beta coefficients (Wald CIs) across band-combination models,
# plus the all-bands CV row pulled from `model_results`. pielou_evenness is
# dropped per the manuscript's figure design.
make_figure_6 <- function(model_results, cv_band_results) {
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

  all_cv_results %>%
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
    labs(x = expression(beta * " co-efficient"), y = NULL) +
    guides(color = "none") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(2, "lines"),
      text             = element_text(size = 20),
      strip.text       = element_text(size = 18, face = "bold"),
      axis.text.y      = element_text(size = 16),
      axis.text.x      = element_text(size = 14),
      axis.title.x     = element_text(size = 18),
      axis.title.y     = element_blank(),
      plot.background  = element_rect(fill = "white", color = "white")
    )
}


# Workflow diagram — hand-built spec (not auto-generated from tar_network()).
#
# The pipeline DAG has ~80 nodes; the paper figure needs ~13. Every conceptual
# box here is a manual collapse of one or more `tar_target`s, so the figure
# stays the same shape even as 4c/4d expand the implementation graph. The
# anti-drift property comes from the consistency check at the bottom: if
# `_targets.R` adds a target that doesn't fit any existing conceptual box, the
# function warns so future-you decides whether the figure needs updating.
#
# `version = "manuscript"` — preprocessing flows into the analysis without a
# boundary; reads as part of the methods narrative.
# `version = "readme"`     — preprocessing wrapped in a dashed rectangle with
# "done once; outputs archived on Zenodo" so users reproducing the code know
# they download masked rasters and skip preprocessing entirely.
#
# Requires ggraph + tidygraph + ggforce loaded.
make_figure_workflow <- function(version = c("manuscript", "readme")) {
  version <- match.arg(version)

  nodes <- tibble::tribble(
    ~name,              ~label,                                                                             ~x,  ~y, ~group,
    "drone_tiles",      "Drone tiles\n(4 sites)",                                                          -2,  7,  "preprocess",
    "training",         "Visually classified\ntraining pixels\n(veg / non-veg,\nshadow / non-shadow)",      2,  7,  "preprocess",
    "preprocess",       "Pre-processing\n(Pix4D ortho +\nillumination correction +\nband stacking)",      -2,  6,  "preprocess",
    "thresholds",       "ROC curves →\nper-site thresholds\n(NIR, NDVI)",                                   2,  6,  "preprocess",
    "masked_rasters",   "Masked multiband\nraster\n(4 sites, 5 bands)",                                     0,  5,  "preprocess",
    "fishnets",         "Fishnets\n(5×5 grid\nper site)",                                                -2.5,  5,  "input",
    "survey",           "AusPlots survey\n(March 2024)",                                                  2.5,  5,  "input",
    "spectral_metrics", "Spectral metrics\nCV / SV / 5D CHV\nrarefied n=999, seed=42",                   -2.5,  4,  "analysis",
    "spectral_species", "Spectral species\nRF + k-means k=40\n×20 seeds averaged →\nrichness, H', 1/D",     0,  4,  "analysis",
    "taxonomic",        "Taxonomic diversity\nS, exp(H'), 1/D, J'",                                       2.5,  4,  "analysis",
    "cv_subsets",       "CV across\nband subsets (×3)",                                                  -2.5,  3,  "analysis",
    "models",           "Mixed models per pair\ny ~ x + (1|site)\nsingular-fit refit to lm\nwhen site variance ≈ 0", 0, 2, "model",
    "figures",          "Figure 5: spectral × taxonomic\nFigure 6: CV β across bands",                      0,  1,  "output"
  )

  edges <- tibble::tribble(
    ~from,              ~to,
    "drone_tiles",      "preprocess",
    "training",         "thresholds",
    "preprocess",       "masked_rasters",
    "thresholds",       "masked_rasters",
    "masked_rasters",   "spectral_metrics",
    "masked_rasters",   "spectral_species",
    "fishnets",         "spectral_metrics",
    "fishnets",         "spectral_species",
    "survey",           "taxonomic",
    "spectral_metrics", "cv_subsets",
    "spectral_metrics", "models",
    "cv_subsets",       "models",
    "spectral_species", "models",
    "taxonomic",        "models",
    "models",           "figures"
  )

  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  layout <- ggraph::create_layout(graph, layout = "manual",
                                  x = nodes$x, y = nodes$y)

  p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        start_cap = ggraph::label_rect(node1.label,
                                       padding  = ggplot2::margin(2.5, 3.5, 2.5, 3.5, "mm"),
                                       fontsize = 8.5),
        end_cap   = ggraph::label_rect(node2.label,
                                       padding  = ggplot2::margin(2.5, 3.5, 2.5, 3.5, "mm"),
                                       fontsize = 8.5)
      ),
      arrow      = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed"),
      colour     = "grey40",
      edge_width = 0.5
    ) +
    ggraph::geom_node_label(
      ggplot2::aes(label = label, fill = group),
      colour        = "black",
      label.padding = grid::unit(c(2, 3, 2, 3), "mm"),
      lineheight    = 0.9,
      size          = 3
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        preprocess = "grey92",
        input      = "lightyellow",
        analysis   = "#cfe6f4",
        model      = "#cce8cc",
        output     = "wheat"
      ),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 1.2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
      plot.margin     = grid::unit(c(6, 6, 6, 6), "mm")
    )

  if (version == "readme") {
    # Hand-placed dashed boundary + label. Tight rect around the four
    # preprocessing nodes at y=6 and y=7 (drone_tiles, training, preprocess,
    # thresholds) plus masked_rasters at y=5. The label sits in the
    # left margin of the rect's top so it doesn't overlap nodes.
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin     = -3.4, xmax = 3.4,
        ymin     =  4.4, ymax = 7.6,
        fill     = NA,
        colour   = "grey40",
        linetype = "dashed",
        linewidth = 0.4
      ) +
      ggplot2::annotate(
        "text",
        x      = -3.4, y = 7.6,
        hjust  = 0,    vjust = -0.2,
        label  = "Preprocessing — done once; outputs archived on Zenodo",
        colour = "grey25",
        size   = 3.2
      )
  }

  p
}

# Compute spectral-species per-subplot summaries (richness + Shannon + Simpson)
# for a single RNG seed across all supplied sites. Pure: takes paths in, returns
# one data.frame, no file I/O. The driver (or a Phase 4d tar_map_rep over seeds)
# handles persistence.
#
# Pipeline per seed:
#   1. Random-sample `sample_size` pixels per raster; combine across sites.
#   2. Train RF with proximity matrix on the combined sample (the "spectral
#      species" signature step).
#   3. k-means cluster the proximity matrix into `k_clusters` groups.
#   4. Train a second RF mapping pixel values -> cluster label, then predict the
#      cluster for every pixel of each raster.
#   5. Aggregate per-pixel cluster labels to per-subplot richness + Shannon +
#      Simpson via the fishnet polygons.
#
# Inputs are file paths (not opened SpatRaster objects) because terra objects
# don't serialise cleanly across tar_target boundaries — matching the Phase 4d
# branch signature avoids a second refactor when targets lands.
spectral_species_one_seed <- function(image_paths, fishnet_paths, seed,
                                      k_clusters = 40, sample_size = 1250) {
  set.seed(seed)

  rasters <- lapply(image_paths, terra::rast)

  # Stage 1: sample pixels per raster, train RF + proximity matrix
  sample_list <- lapply(rasters, function(r) {
    terra::spatSample(r, size = sample_size, method = "random",
                      na.rm = TRUE, as.df = TRUE)
  })
  combined_samples <- do.call(rbind, sample_list)

  rf_model <- randomForest::randomForest(x = combined_samples,
                                         proximity = TRUE, ntree = 500)
  prox_matrix <- rf_model$proximity

  # Stage 2: cluster proximity, train a per-pixel cluster classifier
  km <- kmeans(prox_matrix, centers = k_clusters, nstart = 10)
  combined_samples$cluster <- km$cluster
  rf_cluster_model <- randomForest::randomForest(
    x = combined_samples[, 1:(ncol(combined_samples) - 1)],
    y = as.factor(combined_samples$cluster),
    ntree = 500
  )

  # Stage 3+4+5: predict + aggregate per site
  per_site_results <- vector("list", length(rasters))
  for (i in seq_along(rasters)) {
    r <- rasters[[i]]
    fishnet_path <- fishnet_paths[[i]]

    r_cluster <- terra::predict(r, rf_cluster_model, type = "response",
                                na.rm = TRUE)

    subplot <- sf::read_sf(fishnet_path) %>%
      dplyr::select(geometry) %>%
      dplyr::mutate(subplot_id = unlist(lapply(1:5, function(x) paste(x, 1:5, sep = "_")))) %>%
      terra::vect()

    cluster_vals <- terra::extract(r_cluster, subplot)
    colnames(cluster_vals)[2] <- "cluster"
    cluster_vals$subplot_id <- subplot$subplot_id[cluster_vals$ID]
    cluster_vals <- cluster_vals %>% filter(!is.na(cluster))

    spectral_richness <- cluster_vals %>%
      group_by(subplot_id) %>%
      summarise(spectral_species_richness = n_distinct(cluster), .groups = "drop")

    community_matrix <- cluster_vals %>%
      group_by(subplot_id, cluster) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = cluster, values_from = count, values_fill = 0) %>%
      tibble::column_to_rownames("subplot_id")

    per_site_results[[i]] <- data.frame(
      subplot_id        = rownames(community_matrix),
      spectral_richness = spectral_richness$spectral_species_richness[
        match(rownames(community_matrix), spectral_richness$subplot_id)
      ],
      shannon_spectral  = vegan::diversity(community_matrix, index = "shannon"),
      simpson_spectral  = vegan::diversity(community_matrix, index = "simpson"),
      seed              = seed,
      site              = substr(basename(fishnet_path), 1, 10)
    )
  }

  do.call(rbind, per_site_results)
}
