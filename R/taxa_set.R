# TaxonSet ---------------------------------------------------------------

#' Create a TaxonSet
#'
#' A TaxonSet defines a clade- or user-specified set of taxa used to constrain
#' ortholog resolution and reference lookups. This prevents accidental global
#' searches across all taxa.
#'
#' @param name Character. A label for the set (e.g., "Mammalia").
#' @param taxids Integer vector of NCBI taxonomy IDs (optional).
#' @param species_codes Character vector of species codes (optional; e.g., "mmusculus").
#' @param source Character. Provenance of the TaxonSet (default "user").
#' @return An object of class `TaxonSet`.
#' @export
taxa_set <- function(name, taxids = NULL, species_codes = NULL, source = "user") {
  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single character string.")
  }
  if (!is.null(taxids) && !is.numeric(taxids)) {
    stop("`taxids` must be numeric/integer if provided.")
  }
  if (!is.null(species_codes) && !is.character(species_codes)) {
    stop("`species_codes` must be character if provided.")
  }

  structure(
    list(
      name = name,
      taxids = taxids,
      species_codes = species_codes,
      source = source,
      date = as.character(Sys.Date())
    ),
    class = "TaxonSet"
  )
}

#' @export
print.TaxonSet <- function(x, ...) {
  cat("<TaxonSet>\n")
  cat("  name: ", x$name, "\n", sep = "")
  cat("  taxids: ", if (is.null(x$taxids)) "NULL" else length(x$taxids), "\n", sep = "")
  cat("  species_codes: ", if (is.null(x$species_codes)) "NULL" else length(x$species_codes), "\n", sep = "")
  cat("  source: ", x$source, "\n", sep = "")
  cat("  date: ", x$date, "\n", sep = "")
  invisible(x)
}

#' Create a TaxonSet from taxids
#' @param taxids Integer vector of NCBI taxonomy IDs.
#' @param name Optional name; defaults to "taxids".
#' @return A `TaxonSet`.
#' @export
taxa_set_from_taxids <- function(taxids, name = "taxids") {
  taxa_set(name = name, taxids = taxids, species_codes = NULL, source = "taxids")
}

#' Create a TaxonSet from species codes
#' @param species_codes Character vector of species codes (e.g., "mmusculus").
#' @param name Optional name; defaults to "species".
#' @return A `TaxonSet`.
#' @export
taxa_set_from_species <- function(species_codes, name = "species") {
  taxa_set(name = name, taxids = NULL, species_codes = species_codes, source = "species_codes")
}
