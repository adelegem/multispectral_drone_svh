library(testthat)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# Resolve project root regardless of where tests are invoked from
project_root <- normalizePath(file.path(dirname(testthat::test_path()), ".."))

source(file.path(project_root, "funx.R"))

test_that("calculate_field_diversity output matches baseline", {
  survey <- read_csv(
    file.path(project_root, "data", "ausplots_march_24.csv"),
    show_col_types = FALSE
  )

  result <- calculate_field_diversity(survey)$final_results

  baseline <- readRDS(file.path(project_root, "tests", "fixtures",
                                "taxonomic_diversity_baseline.rds"))

  expect_equal(result, baseline)
})
