# Core API ----------------------------------------------------------------

#' Map rsIDs to orthologous nonsynonymous variants (v0.1 scaffold)
#'
#' This is the main user-facing entry point. In v0.1 it:
#'   1) filters a user-supplied human variant annotation table to rsIDs of interest
#'   2) resolves target orthologs for implicated human genes
#'
#' IMPORTANT: Taxonomic scope must be constrained by supplying either
#' `targets` or `taxa`. Global searches are not allowed.
#'
#' @param rsids Character vector of rsIDs.
#' @param targets Character vector of target species codes (optional).
#' @param taxa A `TaxonSet` object (optional).
#' @param variant_db Human consequence table (e.g., VEP-derived). Must satisfy
#'   `validate_variant_db()` schema.
#' @param ortholog_db Ortholog mapping table. Must satisfy `validate_ortholog_db()` schema.
#' @param transcript_policy Transcript selection policy passed to `snp_rsids_to_consequence()`.
#' @param protein_coding_only Logical; passed to `snp_rsids_to_consequence()`.
#' @param one_to_many Policy for 1:many orthologs: "keep_all", "best_hit", "drop_ambiguous".
#' @param min_score Optional score threshold passed to `ortholog_resolve()` (if score column exists).
#' @param format Output format: "long" (default) or "nested".
#' @param ref_db_targets Optional reference database for target species (reserved for future use).
#' @param map_coords Logical; if TRUE, attach coordinate-mapping scaffold columns (requires `target_protein_id` for full support).
#' @param cache Logical; if TRUE, use the coordinate-mapping cache when `map_coords = TRUE`.
#' @param cache_dir Optional cache directory; defaults to a user cache directory via `tools::R_user_dir()`.
#' @return A tibble of class `SNPmineResult`.
#' @export
snp_map_ns <- function(rsids,
                       targets = NULL,
                       taxa = NULL,
                       variant_db,
                       ortholog_db,
                       ref_db_targets = NULL,
                       map_coords = FALSE,
                       cache = TRUE,
                       cache_dir = NULL,
                       transcript_policy = c("canonical", "all", "pick_longest"),
                       protein_coding_only = TRUE,
                       one_to_many = c("keep_all", "best_hit", "drop_ambiguous"),
                       min_score = NULL,
                       format = c("long", "nested")) {
  format <- match.arg(format)
  transcript_policy <- match.arg(transcript_policy)
  one_to_many <- match.arg(one_to_many)

  # Enforce scoping up front (same rule as before)
  if (is.null(targets) && is.null(taxa)) {
    stop("Provide either `targets` (species codes) or `taxa` (TaxonSet) to constrain the search.")
  }

  # Step 1: rsIDs -> human coding consequences (tidy tibble)
  cons <- snp_rsids_to_consequence(
    rsids = rsids,
    variant_db = variant_db,
    transcript_policy = transcript_policy,
    protein_coding_only = protein_coding_only
  )

  # If no consequences found, return an empty-but-consistent result
  if (nrow(cons) == 0) {
    out <- tibble::tibble(
      rsid = character(0),
      gene_id = character(0),
      transcript_id = character(0),
      protein_id = character(0),
      protein_pos = numeric(0),
      aa_ref = character(0),
      aa_alt = character(0),
      consequence = character(0),
      target_species = character(0),
      target_gene_id = character(0)
    )
    class(out) <- c("SNPmineResult", class(out))
    return(out)
  }

  # Step 2: human genes -> orthologs (tidy tibble)
  ortho <- ortholog_resolve(
    human_gene_ids = unique(cons$gene_id),
    targets = targets,
    taxa = taxa,
    ortholog_db = ortholog_db,
    one_to_many = one_to_many,
    min_score = min_score
  )

  # Join consequence rows to ortholog rows via human gene
  # (ortholog_resolve uses `human_gene_id`; cons uses `gene_id`)
  if (nrow(ortho) == 0) {
    out <- tibble::as_tibble(cons)
    out$target_species <- NA_character_
    out$target_gene_id <- NA_character_
    class(out) <- c("SNPmineResult", class(out))
    return(out)
  }

  merged <- merge(
    x = tibble::as_tibble(cons),
    y = tibble::as_tibble(ortho),
    by.x = "gene_id",
    by.y = "human_gene_id",
    all.x = TRUE
  )
  # Optional: coordinate mapping stub (protein_pos only for now)
  if (isTRUE(map_coords)) {

    # For v0.1 stub we do not require sequences; we only require protein_id presence.
    if (!("protein_id" %in% names(merged))) {
      stop("map_coords = TRUE requires `protein_id` in the consequence table (variant_db).")
    }

    # If your ortholog_db already includes target_protein_id, we can use it.
    # Otherwise, we cannot map protein coordinates yet; fill NAs with explicit status.
    if (!("target_protein_id" %in% names(merged))) {
      merged$target_protein_id <- NA_character_
      merged$target_protein_pos <- NA_real_
      merged$coord_confidence <- 0
      merged$coord_status <- "stub_no_target_protein_id"
    } else {
      # Call stub mapper rowwise (safe now; later we can vectorize)
      mapped <- mapply(
        FUN = function(hp, tp, pos) {
          if (is.na(tp) || is.na(hp)) {
            return(c(NA_real_, 0, "stub_missing_protein_id"))
          }
          tmp <- coord_map_protein(
            human_protein_id = hp,
            target_protein_id = tp,
            human_protein_pos = pos,
            method = "stub",
            cache = cache,
            cache_dir = cache_dir
          )
          c(tmp$target_protein_pos[1], tmp$confidence[1], tmp$status[1])
        },
        hp = merged$protein_id,
        tp = merged$target_protein_id,
        pos = merged$protein_pos,
        SIMPLIFY = TRUE
      )

      merged$target_protein_pos <- as.numeric(mapped[1, ])
      merged$coord_confidence <- as.numeric(mapped[2, ])
      merged$coord_status <- as.character(mapped[3, ])
    }
  }

  # Keep a clean, predictable column order
  preferred <- c(
    "rsid","gene_id","transcript_id","protein_id","protein_pos","aa_ref","aa_alt","consequence",
    "target_species","target_gene_id","score","ortholog_type","evidence","source"
  )
  keep <- intersect(preferred, names(merged))
  merged <- merged[, c(keep, setdiff(names(merged), keep)), drop = FALSE]

  if (format == "nested") {
    # one row per rsid with list-column of mappings
    nested <- split(merged, merged$rsid)
    out <- tibble::tibble(
      rsid = names(nested),
      mappings = unname(nested)
    )
    class(out) <- c("SNPmineResult", class(out))
    return(out)
  }

  class(merged) <- c("SNPmineResult", class(merged))
  merged
}
