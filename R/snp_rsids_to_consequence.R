# rsID -> consequence ------------------------------------------------------

#' Convert rsIDs to coding consequences (scaffold)
#'
#' This function filters a user-supplied human variant annotation table
#' (e.g., Ensembl VEP output that has been converted to a tibble) down to the
#' variants of interest and applies a transcript selection policy.
#'
#' SNPmineR deliberately avoids live queries; instead, users provide `variant_db`.
#'
#' @param rsids Character vector of rsIDs.
#' @param variant_db Human variant consequence table with required columns.
#'   See `validate_variant_db()` for schema.
#' @param transcript_policy Transcript selection policy:
#'   - "canonical": keep rows where `is_canonical == TRUE` if that column exists;
#'                  otherwise keep all and warn.
#'   - "all": keep all transcript consequences.
#'   - "pick_longest": placeholder for future logic; currently keeps all and warns.
#' @param protein_coding_only Logical; if TRUE and `is_protein_coding` exists, filter to TRUE.
#' @return A tibble of class `SNPmineHumanConsequence`.
#' @export
snp_rsids_to_consequence <- function(rsids,
                                     variant_db,
                                     transcript_policy = c("canonical", "all", "pick_longest"),
                                     protein_coding_only = TRUE) {
  transcript_policy <- match.arg(transcript_policy)
  snp_validate_rsids(rsids)
  validate_variant_db(variant_db)

  x <- variant_db[variant_db$rsid %in% rsids, , drop = FALSE]
  x <- tibble::as_tibble(x)

  if (nrow(x) == 0) {
    # Return consistent empty tibble with expected columns
    out <- tibble::tibble(
      rsid = character(0),
      gene_id = character(0),
      transcript_id = character(0),
      protein_id = character(0),
      protein_pos = numeric(0),
      aa_ref = character(0),
      aa_alt = character(0),
      consequence = character(0)
    )
    class(out) <- c("SNPmineHumanConsequence", class(out))
    return(out)
  }

  # Optional filter: protein coding
  if (protein_coding_only && "is_protein_coding" %in% names(x)) {
    x <- x[isTRUE(x$is_protein_coding) | x$is_protein_coding == TRUE, , drop = FALSE]
  }

  # Transcript policy handling
  if (transcript_policy == "canonical") {
    if ("is_canonical" %in% names(x)) {
      x <- x[isTRUE(x$is_canonical) | x$is_canonical == TRUE, , drop = FALSE]
    } else {
      warning("transcript_policy = 'canonical' requested but `is_canonical` column not found; keeping all transcripts.")
    }
  } else if (transcript_policy == "pick_longest") {
    warning("transcript_policy = 'pick_longest' not yet implemented; keeping all transcripts.")
  }

  # Keep only essential columns (plus optional helpful ones)
  keep <- intersect(
    c("rsid","gene_id","transcript_id","protein_id","protein_pos","aa_ref","aa_alt","consequence",
      "impact","is_canonical","is_protein_coding"),
    names(x)
  )
  x <- x[, keep, drop = FALSE]

  class(x) <- c("SNPmineHumanConsequence", class(x))
  x
}
