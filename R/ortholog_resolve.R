# Ortholog resolution ------------------------------------------------------

#' Resolve orthologs for human genes in target taxa
#'
#' Returns a tidy mapping of human genes to target orthologs, with explicit
#' handling of 1:many relationships and strict taxonomic scoping.
#'
#' @param human_gene_ids Character vector of human gene identifiers.
#' @param targets Character vector of target species codes (optional).
#' @param taxa A `TaxonSet` (optional). Used to constrain target species.
#' @param ortholog_db Ortholog mapping table. Must include: human_gene_id,
#'   target_species, target_gene_id. Optional: score, ortholog_type, evidence.
#' @param one_to_many Policy for 1:many orthologs: "keep_all", "best_hit",
#'   or "drop_ambiguous".
#' @param min_score Optional numeric threshold; if provided and score column
#'   exists, filter out mappings below this.
#' @return A tibble with ortholog mappings.
#' @export
ortholog_resolve <- function(human_gene_ids,
                             targets = NULL,
                             taxa = NULL,
                             ortholog_db,
                             one_to_many = c("keep_all", "best_hit", "drop_ambiguous"),
                             min_score = NULL) {
  one_to_many <- match.arg(one_to_many)

  if (!is.character(human_gene_ids) || length(human_gene_ids) == 0) {
    stop("`human_gene_ids` must be a non-empty character vector.")
  }
  if (is.null(targets) && is.null(taxa)) {
    stop("Provide either `targets` (species codes) or `taxa` (TaxonSet) to constrain the search.")
  }
  if (!is.null(taxa) && !inherits(taxa, "TaxonSet")) {
    stop("`taxa` must be a TaxonSet object created by taxa_set().")
  }
  if (!is.null(targets) && !is.character(targets)) {
    stop("`targets` must be a character vector of species codes.")
  }

  validate_ortholog_db(ortholog_db)

  targets_resolved <- targets %||% taxa$species_codes %||% character(0)

  # If user provided only taxids (species_codes NULL), we cannot filter by species here.
  # We still enforce the "no global search" rule, but allow pass-through with warning.
  filter_by_species <- length(targets_resolved) > 0

  db <- ortholog_db

  # Filter to requested human genes
  db <- db[db$human_gene_id %in% human_gene_ids, , drop = FALSE]

  # Filter by targets if available
  if (filter_by_species) {
    db <- db[db$target_species %in% targets_resolved, , drop = FALSE]
  }

  # Optional score filtering
  if (!is.null(min_score) && "score" %in% names(db)) {
    if (!is.numeric(db$score)) {
      stop("`ortholog_db$score` must be numeric if present and `min_score` is used.")
    }
    db <- db[db$score >= min_score, , drop = FALSE]
  }

  # Convert to tibble early for consistent printing downstream
  out <- tibble::as_tibble(db)

  # Handle one-to-many relationships per (human_gene_id, target_species)
  if (nrow(out) == 0) {
    # Return a consistent empty tibble with key columns
    return(tibble::tibble(
      human_gene_id = character(0),
      target_species = character(0),
      target_gene_id = character(0)
    ))
  }

  # Determine multiplicity groups
  key <- paste(out$human_gene_id, out$target_species, sep = "||")
  counts <- table(key)
  is_ambig <- counts[key] > 1

  if (one_to_many == "keep_all") {
    return(out)
  }

  if (one_to_many == "drop_ambiguous") {
    return(out[!is_ambig, , drop = FALSE])
  }

  # one_to_many == "best_hit"
  if ("score" %in% names(out)) {
    if (!is.numeric(out$score)) {
      stop("`ortholog_db$score` must be numeric to use one_to_many = 'best_hit'.")
    }
    # Choose max score within each (human_gene_id, target_species)
    out$..key <- key
    # base R split/apply to avoid dplyr dependency here
    idx <- unlist(lapply(split(seq_len(nrow(out)), out$..key), function(ii) {
      ii[which.max(out$score[ii])]
    }), use.names = FALSE)
    out <- out[idx, , drop = FALSE]
    out$..key <- NULL
    return(out)
  }

  # If no score column, best_hit cannot be determined reliably.
  stop("one_to_many = 'best_hit' requires a numeric `score` column in `ortholog_db`.")
}
