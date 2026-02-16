# Coordinate mapping cache (internal) --------------------------------------

#' Default cache directory for SNPmineR
#' @noRd
coord_cache_dir <- function(cache_dir = NULL) {
  if (!is.null(cache_dir)) return(normalizePath(cache_dir, winslash = "/", mustWork = FALSE))
  # Base R; no new dependency
  d <- tools::R_user_dir("SNPmineR", which = "cache")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  normalizePath(d, winslash = "/", mustWork = FALSE)
}

#' Build a stable cache key
#' @noRd
coord_cache_key <- function(human_protein_id, target_protein_id, method = "stub") {
  # Keep it simple + filesystem-safe
  paste0(method, "__", human_protein_id, "__", target_protein_id)
}

#' Get cached object (or NULL)
#' @noRd
coord_cache_get <- function(key, cache_dir = NULL) {
  d <- coord_cache_dir(cache_dir)
  f <- file.path(d, paste0(key, ".rds"))
  if (!file.exists(f)) return(NULL)
  readRDS(f)
}

#' Save object to cache
#' @noRd
coord_cache_set <- function(key, value, cache_dir = NULL) {
  d <- coord_cache_dir(cache_dir)
  f <- file.path(d, paste0(key, ".rds"))
  saveRDS(value, f)
  invisible(f)
}
