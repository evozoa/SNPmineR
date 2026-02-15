required_cols <- function(df, cols, object_name) {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "`%s` is missing required columns: %s",
        object_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate variant database schema
#'
#' @param variant_db A data frame/tibble containing human variant annotations.
#'
#' @return Invisibly returns `TRUE` if valid.
#' @export
validate_variant_db <- function(variant_db) {
  if (!is.data.frame(variant_db)) {
    stop("`variant_db` must be a data.frame or tibble.", call. = FALSE)
  }

  required_cols(
    variant_db,
    c("rsid", "human_gene_id", "transcript_id", "consequence", "aa_change"),
    "variant_db"
  )
}

#' Validate ortholog database schema
#'
#' @param ortholog_db A data frame/tibble mapping human genes to target orthologs.
#'
#' @return Invisibly returns `TRUE` if valid.
#' @export
validate_ortholog_db <- function(ortholog_db) {
  if (!is.data.frame(ortholog_db)) {
    stop("`ortholog_db` must be a data.frame or tibble.", call. = FALSE)
  }

  required_cols(
    ortholog_db,
    c("human_gene_id", "target_gene_id", "target_species", "target_taxid", "orthology_confidence"),
    "ortholog_db"
  )
}

#' Validate reference database schema
#'
#' @param refdb A data frame/tibble of target species reference records.
#'
#' @return Invisibly returns `TRUE` if valid.
#' @export
validate_ref_db <- function(refdb) {
  if (!is.data.frame(refdb)) {
    stop("`refdb` must be a data.frame or tibble.", call. = FALSE)
  }

  required_cols(
    refdb,
    c("target_species", "target_taxid", "ref_source", "build"),
    "refdb"
  )
}
