#' Validate rsID vector format
#'
#' @param rsids Character vector of rsIDs.
#'
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
validate_rsids <- function(rsids) {
  if (missing(rsids) || is.null(rsids) || length(rsids) == 0) {
    stop("`rsids` must be a non-empty character vector of rsIDs.", call. = FALSE)
  }
  if (!is.character(rsids)) {
    stop("`rsids` must be a character vector.", call. = FALSE)
  }

  bad <- !grepl("^rs[0-9]+$", rsids)
  if (any(is.na(rsids)) || any(bad)) {
    stop(
      "All `rsids` must match pattern '^rs[0-9]+$' and be non-missing.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_targets_or_taxa <- function(targets = NULL, taxa = NULL, context = "This operation") {
  has_targets <- !is.null(targets) && length(targets) > 0
  has_taxa <- !is.null(taxa)

  if (!has_targets && !has_taxa) {
    stop(
      paste0(
        context,
        " requires taxonomic scoping via either `targets` or `taxa`. ",
        "Searching across all NCBI taxa is not allowed by default."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
