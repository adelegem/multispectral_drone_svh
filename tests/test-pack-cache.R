# Test that tools/pack_cache.R produces a lightweight cache tarball with the
# expected shape. Skips when there's no _targets/ cache to package — so it's
# a no-op on a fresh checkout and on CI runners, and an automatic check
# every time you run `testthat::test_dir("tests")` after a tar_make().

project_root <- normalizePath(file.path(dirname(testthat::test_path()), ".."))

skip_if_no_clean_cache <- function() {
  if (!file.exists(file.path(project_root, "_targets", "meta", "meta"))) {
    testthat::skip("No _targets/ cache — run tar_make() first")
  }
  # tar_outdated() reads the cache via cwd-relative paths
  old_wd <- setwd(project_root); on.exit(setwd(old_wd), add = TRUE)
  outdated <- tryCatch(targets::tar_outdated(), error = function(e) "error")
  if (identical(outdated, "error")) {
    testthat::skip("targets metadata unreadable")
  }
  if (length(outdated) > 0) {
    testthat::skip(sprintf(
      "Cache outdated (%d targets) — pack_cache.R would refuse anyway", length(outdated)
    ))
  }
}

test_that("pack_cache.R produces _targets_lite.tar.gz with the expected shape", {
  skip_if_no_clean_cache()

  tarball_path <- file.path(project_root, "_targets_lite.tar.gz")
  preexists    <- file.exists(tarball_path)
  # Only delete on exit if we created it; don't clobber a pre-existing one.
  if (!preexists) on.exit(unlink(tarball_path), add = TRUE)

  # pack_cache.R uses cwd-relative paths so we run it from the project root.
  old_wd <- setwd(project_root); on.exit(setwd(old_wd), add = TRUE)
  capture.output(source("tools/pack_cache.R", local = new.env()))

  expect_true(file.exists(tarball_path))
  contents <- untar(tarball_path, list = TRUE)

  # Meta ledger is load-bearing — without it tar_make() can't use the
  # cached objects at all.
  expect_true("meta/meta" %in% contents,
              info = "tarball must include _targets/meta/meta")

  # The whole point of this artifact: no pixel_values_* objects (the 3.6 GB
  # of raw pixel matrices we want readers to rebuild locally).
  pixel_in_tar <- grep("^objects/pixel_values", contents, value = TRUE)
  expect_length(pixel_in_tar, 0L)

  # Every non-pixel object in the local cache must be in the tarball.
  all_objects     <- list.files("_targets/objects/")
  local_non_pixel <- all_objects[!grepl("^pixel_values", all_objects)]
  missing         <- setdiff(paste0("objects/", local_non_pixel), contents)
  expect_length(missing, 0L)

  # Size budget: well under GitHub's 2 GB release-asset limit. If this fails,
  # pixel_values_* (or something else huge) leaked into the tarball.
  size_mb <- file.info(tarball_path)$size / 1024^2
  expect_lt(size_mb, 50)
})
