#' Load a scoped reference database (stub)
#'
#' This scaffold intentionally requires explicit taxonomic constraints.
#'
#' @param taxa Optional `TaxonSet` object.
#' @param targets Optional character vector of species codes.
#' @param path Optional local path to serialized reference data.
#'
#' @return A tibble of reference records filtered to requested taxa/targets.
#' @export
#'
#' @examples
#' refdb_load(targets = "mmusculus")
refdb_load <- function(taxa = NULL, targets = NULL, path = NULL) {
  assert_targets_or_taxa(targets = targets, taxa = taxa, context = "`refdb_load()`")

  refdb <- if (!is.null(path) && file.exists(path)) {
    readRDS(path)
  } else {
    tibble::tibble(
      target_species = character(0),
      target_taxid = integer(0),
      ref_source = character(0),
      build = character(0)
    )
  }

  validate_ref_db(refdb)
  refdb_filter(refdb, taxa = taxa, targets = targets)
}

#' Filter a reference database by taxonomic scope
#'
#' @param refdb A validated reference database tibble.
#' @param taxa Optional `TaxonSet` object.
#' @param targets Optional character vector of species codes.
#'
#' @return Filtered tibble.
#' @export
#'
#' @examples
#' refdb <- tibble::tibble(
#'   target_species = c("mmusculus", "drerio"),
#'   target_taxid = c(10090L, 7955L),
#'   ref_source = c("toy", "toy"),
#'   build = c("v1", "v1")
#' )
#' refdb_filter(refdb, targets = "mmusculus")
refdb_filter <- function(refdb, taxa = NULL, targets = NULL) {
  validate_ref_db(refdb)
  assert_targets_or_taxa(targets = targets, taxa = taxa, context = "`refdb_filter()`")

  out <- refdb
  if (!is.null(targets)) {
    out <- out[out$target_species %in% targets, , drop = FALSE]
  }

  if (inherits(taxa, "TaxonSet")) {
    if (length(taxa$species_codes) > 0) {
      out <- out[out$target_species %in% taxa$species_codes, , drop = FALSE]
    }
    if (length(taxa$taxids) > 0) {
      out <- out[out$target_taxid %in% taxa$taxids, , drop = FALSE]
    }
  }

  tibble::as_tibble(out)
}
