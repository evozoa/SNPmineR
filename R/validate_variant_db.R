# Variant DB validation ----------------------------------------------------

#' Validate a human variant consequence table
#'
#' SNPmineR expects a precomputed human variant annotation table (e.g., from VEP)
#' containing rsIDs and coding consequences. This function validates that the
#' required columns exist and are of appropriate types.
#'
#' @param variant_db A data.frame/tibble with required columns.
#' @return Invisibly returns TRUE; errors if invalid.
#' @keywords internal
validate_variant_db <- function(variant_db) {
  if (!is.data.frame(variant_db)) {
    stop("`variant_db` must be a data.frame or tibble.")
  }

  required <- c(
    "rsid",
    "gene_id",
    "transcript_id",
    "protein_id",
    "protein_pos",
    "aa_ref",
    "aa_alt",
    "consequence"
  )

  missing <- setdiff(required, names(variant_db))
  if (length(missing) > 0) {
    stop("`variant_db` is missing required column(s): ", paste(missing, collapse = ", "))
  }

  # Basic type checks (tolerant but protective)
  if (!is.character(variant_db$rsid)) stop("`variant_db$rsid` must be character.")
  if (!is.character(variant_db$gene_id)) stop("`variant_db$gene_id` must be character.")
  if (!is.character(variant_db$transcript_id)) stop("`variant_db$transcript_id` must be character.")
  if (!is.character(variant_db$protein_id)) stop("`variant_db$protein_id` must be character.")
  if (!is.numeric(variant_db$protein_pos)) stop("`variant_db$protein_pos` must be numeric/integer.")
  if (!is.character(variant_db$aa_ref)) stop("`variant_db$aa_ref` must be character.")
  if (!is.character(variant_db$aa_alt)) stop("`variant_db$aa_alt` must be character.")
  if (!is.character(variant_db$consequence)) stop("`variant_db$consequence` must be character.")

  invisible(TRUE)
}
