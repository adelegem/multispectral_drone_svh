# Print digests of the canonical pipeline outputs.
#
# Used in two places:
#   1. The README's "Verifying your reproduction" section lists these digests
#      so readers can compare against their own `digest::digest()` of the same
#      target after running tar_make().
#   2. At each new release: regenerate this output and update the README +
#      release notes.
#
# Run from the project root, after a full `tar_make()` has populated the cache:
#
#   Rscript tools/checksums.R
#
# Digests are MD5 (digest::digest()'s default) — match the README example.
# Bit-exact reproducibility is contingent on the locked renv environment
# (see README "Reproducibility notes").

library(targets)
library(digest)

canonical_targets <- c(
  "taxonomic_diversity",
  "spectral_taxonomic_diversity",
  "mean_spectral_species",
  "spectral_biodiversity_model_results",
  "cv_biodiversity_model_results",
  "spectral_species_model_results"
)

cat("# Canonical pipeline outputs (MD5 digests via digest::digest())\n\n")
for (n in canonical_targets) {
  d <- digest::digest(tar_read_raw(n))
  cat(sprintf("- `%-40s` %s\n", n, d))
}

cat("\n# Tracked test fixtures\n\n")
fixture_files <- list.files("tests/fixtures", pattern = "\\.rds$", full.names = TRUE)
for (f in fixture_files) {
  d <- digest::digest(readRDS(f))
  cat(sprintf("- `%-50s` %s\n", f, d))
}
