# Ortholog DB validation ---------------------------------------------------

#' Validate an ortholog mapping table
#'
#' @param ortholog_db A data.frame/tibble with required columns.
#' @return Invisibly returns TRUE; errors if invalid.
#' @keywords internal
validate_ortholog_db <- function(ortholog_db) {
  if (!is.data.frame(ortholog_db)) {
    stop("`ortholog_db` must be a data.frame or tibble.")
  }

  required <- c("human_gene_id", "target_species", "target_gene_id")
  missing <- setdiff(required, names(ortholog_db))
  if (length(missing) > 0) {
    stop("`ortholog_db` is missing required column(s): ", paste(missing, collapse = ", "))
  }

  # Optional but recommended columns
  # ortholog_type (e.g., '1:1', '1:many'), score (numeric), evidence/source.
  invisible(TRUE)
}
