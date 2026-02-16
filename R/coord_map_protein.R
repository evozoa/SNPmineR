# Protein coordinate mapping ------------------------------------------------

#' Map protein coordinates between human and a target species (stub)
#'
#' This is a v0.1 stub to stabilize the API. In future versions, this function
#' will align human and target protein sequences and map positions through the
#' alignment. For now, it returns `NA` positions with an explicit status field.
#'
#' Caching is supported via an RDS cache keyed on (human_protein_id, target_protein_id).
#'
#' @param human_protein_id Character scalar. Human protein identifier.
#' @param target_protein_id Character scalar. Target protein identifier.
#' @param human_protein_pos Numeric/integer vector of positions on the human protein.
#' @param method Character. Mapping method identifier (default "stub").
#' @param cache Logical. Use cache if TRUE.
#' @param cache_dir Optional cache directory; defaults to user cache via `tools::R_user_dir()`.
#' @param overwrite Logical. If TRUE, recompute even if cached (useful later).
#' @return A tibble with columns:
#'   human_protein_id, target_protein_id, human_protein_pos, target_protein_pos,
#'   confidence, status.
#' @export
coord_map_protein <- function(human_protein_id,
                              target_protein_id,
                              human_protein_pos,
                              method = "stub",
                              cache = TRUE,
                              cache_dir = NULL,
                              overwrite = FALSE) {

  if (!is.character(human_protein_id) || length(human_protein_id) != 1) {
    stop("`human_protein_id` must be a single character string.")
  }
  if (!is.character(target_protein_id) || length(target_protein_id) != 1) {
    stop("`target_protein_id` must be a single character string.")
  }
  if (!is.numeric(human_protein_pos)) {
    stop("`human_protein_pos` must be numeric/integer.")
  }

  key <- coord_cache_key(human_protein_id, target_protein_id, method = method)

  if (cache && !overwrite) {
    hit <- coord_cache_get(key, cache_dir = cache_dir)
    if (!is.null(hit)) return(hit)
  }

  # ---- STUB LOGIC ----
  out <- tibble::tibble(
    human_protein_id = human_protein_id,
    target_protein_id = target_protein_id,
    human_protein_pos = as.numeric(human_protein_pos),
    target_protein_pos = NA_real_,
    confidence = 0,
    status = "stub_unimplemented"
  )

  if (cache) coord_cache_set(key, out, cache_dir = cache_dir)
  out
}
