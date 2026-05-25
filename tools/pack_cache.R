# Package the {targets} cache for release distribution.
#
# Strategy: ship a "lightweight" cache that excludes pixel_values_* objects
# (~3.6 GB of raw per-subplot pixel matrices). Readers rebuild those in
# ~2 min of compute (~10–15 min wall-clock with R startup); everything
# downstream (the expensive 44 h of 5D CHV + 20-seed spectral-species
# clustering + models + report) is included and verified bit-for-bit
# against the rebuild via targets' content-addressed hashing.
#
# Run from the project root after a successful tar_make():
#
#   Rscript tools/pack_cache.R
#
# Produces _targets_lite.tar.gz in the project root (~2 MB, ~95 files).
# Attach to the GitHub release via:
#
#   gh release upload v1.0.0 _targets_lite.tar.gz

stopifnot(file.exists("_targets/meta/meta"))

# Refuse to package a cache that isn't fully up to date — shipping invalid
# hashes would silently break the reader's verification step.
outdated <- targets::tar_outdated()
if (length(outdated) > 0) {
  stop("Cache contains outdated targets; run tar_make() first.\nOutdated:\n  ",
       paste(outdated, collapse = "\n  "))
}

excluded <- list.files("_targets/objects/", pattern = "^pixel_values")
included <- setdiff(list.files("_targets/objects/"), excluded)

tarball <- "_targets_lite.tar.gz"
status <- system2("tar", c("-C", "_targets", "-czf", tarball,
                            "meta", paste0("objects/", included)))
if (status != 0) stop("tar exited with code ", status)

size_mb <- file.info(tarball)$size / 1024^2
cat(sprintf("Packaged %s: %.2f MB, %d objects.\n",
            tarball, size_mb, length(included)))
cat(sprintf("Excluded %d pixel_values_* objects (~3.6 GB of raw pixel matrices).\n",
            length(excluded)))
cat("Readers rebuild those in ~5 min via tar_make(); downstream stays cached.\n\n")
cat("Upload to the release with:\n")
cat(sprintf("  gh release upload v1.0.0 %s\n", tarball))
