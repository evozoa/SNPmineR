#' Create a taxonomic scope object
#'
#' `TaxonSet` objects capture an explicit set of taxonomic constraints used to
#' avoid global searches over all NCBI taxa.
#'
#' @param x A named list containing `name`, `taxids`, and/or `species_codes`.
#'
#' @return A `TaxonSet` S3 object.
#' @export
#'
#' @examples
#' taxa_set(list(name = "Mammalia", taxids = c(40674, 9685), species_codes = c("hsapiens", "cfamiliaris")))
taxa_set <- function(x) {
  if (!is.list(x)) {
    stop("`x` must be a list.", call. = FALSE)
  }

  name <- x$name %||% "custom"
  taxids <- x$taxids %||% integer(0)
  species_codes <- x$species_codes %||% character(0)
  source <- x$source %||% "user"
  date <- x$date %||% as.character(Sys.Date())

  out <- list(
    name = as.character(name)[1],
    taxids = as.integer(taxids),
    species_codes = as.character(species_codes),
    source = as.character(source)[1],
    date = as.character(date)[1]
  )
  class(out) <- c("TaxonSet", "list")
  out
}

#' Create `TaxonSet` from taxids
#'
#' @param taxids Integer taxonomic IDs.
#' @param name Optional label for the set.
#'
#' @return A `TaxonSet` S3 object.
#' @export
#'
#' @examples
#' taxa_set_from_taxids(c(40674, 7898), name = "Mammalia + Teleost")
taxa_set_from_taxids <- function(taxids, name = "custom_taxids") {
  taxa_set(list(name = name, taxids = taxids, source = "taxids"))
}

#' Create `TaxonSet` from species codes
#'
#' @param species_codes Character species codes (e.g., `"mmusculus"`).
#' @param name Optional label for the set.
#'
#' @return A `TaxonSet` S3 object.
#' @export
#'
#' @examples
#' taxa_set_from_species(c("mmusculus", "rnorvegicus"), name = "Rodents")
taxa_set_from_species <- function(species_codes, name = "custom_species") {
  taxa_set(list(name = name, species_codes = species_codes, source = "species_codes"))
}

#' @export
print.TaxonSet <- function(x, ...) {
  cat("<TaxonSet>", x$name, "\n")
  cat("  taxids:", paste(x$taxids, collapse = ", "), "\n")
  cat("  species_codes:", paste(x$species_codes, collapse = ", "), "\n")
  cat("  source:", x$source, " date:", x$date, "\n")
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
