test_that("snp_map_ns returns tibble long format by default", {
  out <- snp_map_ns(
    rsids = c("rs1", "rs2"),
    targets = "mmusculus",
    variant_db = toy_variant_db(),
    ortholog_db = toy_ortholog_db(),
    ref_db_targets = toy_ref_db()
  )

  expect_s3_class(out, "tbl_df")
  expect_true(all(c("rsid", "human_gene_id", "target_species") %in% names(out)))
})

test_that("snp_map_ns can return nested format", {
  out <- snp_map_ns(
    rsids = c("rs1"),
    targets = "mmusculus",
    variant_db = toy_variant_db(),
    ortholog_db = toy_ortholog_db(),
    ref_db_targets = toy_ref_db(),
    format = "nested"
  )

  expect_s3_class(out, "tbl_df")
  expect_true("mappings" %in% names(out))
})
