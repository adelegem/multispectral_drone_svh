# `preprocessing/`

Scripts that produced the masked multispectral rasters archived on Zenodo ([10.5281/zenodo.17089161](https://doi.org/10.5281/zenodo.17089161)). **Not part of the main `{targets}` pipeline** — these scripts run *upstream* of the analysis and only need to be re-run if the raw drone tiles or training data change.

## Why is this separate from `tar_make()`?

The masked rasters are the *input* to the analysis pipeline. They're produced once per dataset, archived publicly, and downloaded by `tar_make()` via `download_zenodo_rasters()`. Anyone reproducing the analysis grabs the archived rasters directly and skips this stage entirely — see the top-level [README](../README.md)'s *Reproducing the analysis* section.

These scripts exist for **transparency and completeness**: they document how the raw drone outputs were turned into the Zenodo-archived rasters, so reviewers can inspect every step.

## Inputs (not tracked in this repo)

- Per-site **raw drone tiles** captured by the multispectral drone (4 sites: NSABHC0009–0012).
- Per-site **visually classified training pixels** on two axes: veg/non-veg **and** shadow/non-shadow.

These are large and project-specific; contact the authors if you need access for independent verification.

## Pre-step: Pix4D (outside R)

Before any R code runs, **[Pix4D](https://www.pix4d.com/)** orthomosaics each site's drone tiles into a per-band GeoTIFF and applies illumination correction. This is a closed-source GUI workflow and isn't scripted here. The output of Pix4D is one GeoTIFF per band per site, which becomes the input to `01_stack_bands.R` below.

## R scripts (placeholder)

The actual scripts will be added in a follow-up commit. Planned breakdown (3 files), each `source("../funx.R")` for the shared helpers:

1. **`01_stack_bands.R`** — stack the 5 per-band Pix4D outputs into a single multiband GeoTIFF per site, in wavelength order (blue, green, red, red-edge, NIR).
2. **`02_find_thresholds.R`** — apply `funx.R::find_optimum_thresholds()` (ROC against the training pixels) to compute the per-site **NIR + NDVI** thresholds. *A Red filter exists in `saltbush` but was deliberately not applied to this dataset.*
3. **`03_apply_mask.R`** — apply the thresholds via `funx.R::create_masked_raster()` to produce `NSABHC00{09..12}_masked.tif`. These outputs are what gets uploaded to Zenodo and consumed by `tar_make()`.

All three are independent R sessions — they don't share state. The shared helpers they use (`create_masked_raster`, `find_optimum_thresholds`) live in `../funx.R` for re-use by tests, but are not called from the targets pipeline itself.

## How to run (once scripts land)

```sh
# From the project root, with the renv environment activated:
Rscript preprocessing/01_stack_bands.R
Rscript preprocessing/02_find_thresholds.R
Rscript preprocessing/03_apply_mask.R
```

These produce the masked TIFFs locally. Uploading the results to Zenodo (and minting the data DOI) is a manual web step done by the project authors.
