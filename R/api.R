#' Map human rsIDs to orthologous nonsynonymous variants
#'
#' One-liner API for SNP-to-ortholog mapping with explicit taxonomic constraints.
#'
#' @param rsids Character vector of rsIDs.
#' @param targets Optional character species codes.
#' @param taxa Optional `TaxonSet` object.
#' @param variant_db Human variant annotation table.
#' @param ortholog_db Orthology mapping table.
#' @param ref_db_targets Reference database for target species.
#' @param transcript_policy One of `"canonical"`, `"all"`, `"pick_longest"`.
#' @param one_to_many One of `"keep_all"`, `"best_hit"`, `"drop_ambiguous"`.
#' @param format One of `"long"` or `"nested"`.
#' @param min_alignment_confidence Minimum orthology confidence threshold.
#' @param drop_ambiguous Logical; drop ambiguous ortholog rows.
#' @param return_unmapped Logical; include unmapped rsIDs in output.
#'
#' @return A tibble in long or nested format.
#' @export
#'
#' @examples
#' variant_db <- tibble::tibble(
#'   rsid = c("rs1", "rs2"),
#'   human_gene_id = c("GENE1", "GENE2"),
#'   transcript_id = c("tx1", "tx2"),
#'   consequence = c("missense_variant", "missense_variant"),
#'   aa_change = c("A10V", "G20D"),
#'   is_canonical = c(TRUE, TRUE),
#'   transcript_length = c(1000, 900)
#' )
#' ortholog_db <- tibble::tibble(
#'   human_gene_id = c("GENE1", "GENE2"),
#'   target_gene_id = c("M1", "M2"),
#'   target_species = c("mmusculus", "mmusculus"),
#'   target_taxid = c(10090L, 10090L),
#'   orthology_confidence = c(0.95, 0.9)
#' )
#' ref_db_targets <- tibble::tibble(
#'   target_species = "mmusculus",
#'   target_taxid = 10090L,
#'   ref_source = "toy",
#'   build = "v1"
#' )
#' snp_map_ns(c("rs1"), targets = "mmusculus", variant_db = variant_db,
#'            ortholog_db = ortholog_db, ref_db_targets = ref_db_targets)
snp_map_ns <- function(
    rsids,
    targets = NULL,
    taxa = NULL,
    variant_db,
    ortholog_db,
    ref_db_targets,
    transcript_policy = c("canonical", "all", "pick_longest"),
    one_to_many = c("keep_all", "best_hit", "drop_ambiguous"),
    format = c("long", "nested"),
    min_alignment_confidence = 0.8,
    drop_ambiguous = TRUE,
    return_unmapped = FALSE
) {
  transcript_policy <- match.arg(transcript_policy)
  one_to_many <- match.arg(one_to_many)
  format <- match.arg(format)

  assert_targets_or_taxa(targets = targets, taxa = taxa, context = "`snp_map_ns()`")

  human <- snp_rsids_to_consequence(
    rsids = rsids,
    variant_db = variant_db,
    transcript_policy = transcript_policy
  )

  ortho <- ortholog_resolve(
    human_gene_ids = unique(human$human_gene_id),
    targets = targets,
    taxa = taxa,
    ortholog_db = ortholog_db,
    one_to_many = one_to_many
  )

  ref_scoped <- refdb_filter(ref_db_targets, taxa = taxa, targets = targets)

  out <- dplyr::inner_join(human, ortho, by = "human_gene_id")
  out <- dplyr::filter(out, .data$orthology_confidence >= min_alignment_confidence)
  out <- dplyr::inner_join(out, ref_scoped, by = c("target_species", "target_taxid"))

  if (drop_ambiguous || one_to_many == "drop_ambiguous") {
    out <- out |>
      dplyr::add_count(.data$rsid, .data$target_species, name = "n_hits") |>
      dplyr::filter(.data$n_hits == 1L) |>
      dplyr::select(-.data$n_hits)
  }

  if (return_unmapped) {
    missing_rsids <- setdiff(rsids, unique(out$rsid))
    if (length(missing_rsids) > 0) {
      unmapped <- tibble::tibble(
        rsid = missing_rsids,
        human_gene_id = NA_character_,
        transcript_id = NA_character_,
        consequence = NA_character_,
        aa_change = NA_character_,
        target_gene_id = NA_character_,
        target_species = NA_character_,
        target_taxid = NA_integer_,
        orthology_confidence = NA_real_,
        ref_source = NA_character_,
        build = NA_character_
      )
      out <- dplyr::bind_rows(out, unmapped)
    }
  }

  out <- tibble::as_tibble(out)

  if (format == "nested") {
    out <- out |>
      dplyr::group_by(.data$rsid) |>
      dplyr::summarise(mappings = list(dplyr::cur_data_all()), .groups = "drop")
  }

  tibble::as_tibble(out)
}

#' Resolve rsIDs to human coding consequences
#'
#' @param rsids Character vector of rsIDs.
#' @param variant_db Variant annotation table.
#' @param transcript_policy One of `"canonical"`, `"all"`, `"pick_longest"`.
#'
#' @return A tibble with class `SNPmineHumanConsequence`.
#' @export
snp_rsids_to_consequence <- function(rsids, variant_db, transcript_policy = c("canonical", "all", "pick_longest")) {
  transcript_policy <- match.arg(transcript_policy)
  validate_rsids(rsids)
  validate_variant_db(variant_db)

  out <- variant_db[variant_db$rsid %in% rsids, , drop = FALSE]

  if (transcript_policy == "canonical" && "is_canonical" %in% colnames(out)) {
    out <- out[out$is_canonical %in% TRUE, , drop = FALSE]
  }

  if (transcript_policy == "pick_longest" && "transcript_length" %in% colnames(out)) {
    out <- out |>
      dplyr::group_by(.data$rsid) |>
      dplyr::slice_max(.data$transcript_length, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  out <- tibble::as_tibble(out)
  class(out) <- c("SNPmineHumanConsequence", class(out))
  out
}

#' Resolve orthologs with explicit taxonomic constraints
#'
#' Primary place where `targets`/`taxa` scope is enforced.
#'
#' @param human_gene_ids Character vector of human gene IDs.
#' @param targets Optional character species codes.
#' @param taxa Optional `TaxonSet` object.
#' @param ortholog_db Ortholog mapping data.
#' @param one_to_many One of `"keep_all"`, `"best_hit"`, `"drop_ambiguous"`.
#'
#' @return A tibble of scoped orthologs.
#' @export
ortholog_resolve <- function(human_gene_ids,
                             targets = NULL,
                             taxa = NULL,
                             ortholog_db,
                             one_to_many = c("keep_all", "best_hit", "drop_ambiguous")) {
  one_to_many <- match.arg(one_to_many)
  assert_targets_or_taxa(targets = targets, taxa = taxa, context = "`ortholog_resolve()`")
  validate_ortholog_db(ortholog_db)

  out <- ortholog_db[ortholog_db$human_gene_id %in% human_gene_ids, , drop = FALSE]

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

  if (one_to_many == "best_hit") {
    out <- out |>
      dplyr::group_by(.data$human_gene_id, .data$target_species) |>
      dplyr::slice_max(.data$orthology_confidence, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  if (one_to_many == "drop_ambiguous") {
    out <- out |>
      dplyr::add_count(.data$human_gene_id, .data$target_species, name = "n_hits") |>
      dplyr::filter(.data$n_hits == 1L) |>
      dplyr::select(-.data$n_hits)
  }

  tibble::as_tibble(out)
}
