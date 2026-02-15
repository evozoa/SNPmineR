toy_variant_db <- function() {
  tibble::tibble(
    rsid = c("rs1", "rs2", "rs2"),
    human_gene_id = c("GENE1", "GENE2", "GENE2"),
    transcript_id = c("tx1", "tx2", "tx2b"),
    consequence = c("missense_variant", "missense_variant", "missense_variant"),
    aa_change = c("A10V", "G20D", "G20E"),
    is_canonical = c(TRUE, TRUE, FALSE),
    transcript_length = c(1000, 900, 1100)
  )
}

toy_ortholog_db <- function() {
  tibble::tibble(
    human_gene_id = c("GENE1", "GENE2", "GENE2"),
    target_gene_id = c("M1", "M2", "M2_alt"),
    target_species = c("mmusculus", "mmusculus", "drerio"),
    target_taxid = c(10090L, 10090L, 7955L),
    orthology_confidence = c(0.95, 0.90, 0.85)
  )
}

toy_ref_db <- function() {
  tibble::tibble(
    target_species = c("mmusculus", "drerio"),
    target_taxid = c(10090L, 7955L),
    ref_source = c("toy", "toy"),
    build = c("v1", "v1")
  )
}
