# SNPmineR utilities ------------------------------------------------------

#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

#' Validate rsID vector
#' @param rsids Character vector of rsIDs (e.g., "rs429358")
#' @return Invisibly returns TRUE; errors if invalid.
#' @keywords internal
snp_validate_rsids <- function(rsids) {
  if (!is.character(rsids)) stop("`rsids` must be a character vector.")
  if (length(rsids) == 0) stop("`rsids` must have length > 0.")
  bad <- !grepl("^rs[0-9]+$", rsids)
  if (any(bad)) {
    stop("Invalid rsID(s): ", paste(unique(rsids[bad]), collapse = ", "))
  }
  invisible(TRUE)
}
