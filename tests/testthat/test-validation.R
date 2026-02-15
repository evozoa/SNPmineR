test_that("rsid validation enforces expected pattern", {
  expect_error(snp_rsids_to_consequence(c("abc"), variant_db = toy_variant_db()), "\^rs")
  expect_s3_class(snp_rsids_to_consequence(c("rs1"), variant_db = toy_variant_db()), "SNPmineHumanConsequence")
})

test_that("schema validation fails fast", {
  bad_variant <- tibble::tibble(rsid = "rs1")
  expect_error(validate_variant_db(bad_variant), "missing required columns")

  bad_ortholog <- tibble::tibble(human_gene_id = "G1")
  expect_error(validate_ortholog_db(bad_ortholog), "missing required columns")

  bad_ref <- tibble::tibble(target_species = "mmusculus")
  expect_error(validate_ref_db(bad_ref), "missing required columns")
})
