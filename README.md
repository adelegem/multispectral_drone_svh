# multispectral_drone_svh

This repository supports the manuscript *"Applying the spectral variability hypothesis to arid shrublands, using multispectral drone imagery"* by Gemmell et al.

Data DOI: [10.5281/zenodo.17089161](https://doi.org/10.5281/zenodo.17089161) — the four multispectral GeoTIFFs (~4.6 GB total) are auto-downloaded by the analysis on first run.

## Setup

**R version:** 4.5 or later (developed and tested on 4.5.3).

**System dependencies:**

- GDAL and PROJ — required by `sf` and `terra`. On macOS: `brew install gdal proj`. On Debian/Ubuntu: `apt install libgdal-dev libproj-dev`.
- A C++17 compiler — required by `geometry` (convex-hull volume).

**R packages:**

```r
install.packages(c(
  "sf", "terra", "tidyverse", "vegan", "data.table", "performance",
  "glmmTMB", "geometry", "randomForest", "cluster", "pROC", "testthat"
))
```

A `renv` lockfile will be added in a follow-up release for exact dependency reproduction.

## Run

```sh
Rscript continuous_metrics_analysis.R    # CV / SV / CHV + mixed models
Rscript spectral_species_analysis.R      # spectral-species clustering + mixed models
```

Both scripts call `download_zenodo_rasters()` on startup, which fetches the four GeoTIFFs into `data/raster_images/` and skips files already present.

## Tests

```sh
Rscript -e 'testthat::test_dir("tests")'                   # Tier 1, ~1 s
RUN_TIER2=true Rscript -e 'testthat::test_dir("tests")'    # Tier 2, ~minutes; needs rasters downloaded
```

## Reproducibility

Rarefaction draws in `continuous_metrics_analysis.R` are seeded with `RAREFACTION_SEED = 42`. The spectral-species analysis loops over seeds 1–20 internally and averages across them. Both conventions are documented in the manuscript methods.

## License

- Source code: [MIT](LICENSE).
- Data files shipped under `data/` and the multispectral rasters auto-downloaded from Zenodo: [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/), matching the Zenodo deposit.

## Citation

If you use this code, please cite both the software and the dataset. See [`CITATION.cff`](CITATION.cff) for structured metadata (GitHub will render it as a "Cite this repository" button). The dataset DOI is [10.5281/zenodo.17089161](https://doi.org/10.5281/zenodo.17089161).
