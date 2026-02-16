test_that("TaxonSet creation works and stores fields", {
  tx <- taxa_set_from_taxids(c(40674, 10090), name = "Mammalia")
  expect_s3_class(tx, "TaxonSet")
  expect_equal(tx$name, "Mammalia")
  expect_true(all(c(40674L, 10090L) %in% tx$taxids))
})

test_that("refdb filtering honors taxa and targets", {
  tx <- taxa_set_from_species("mmusculus", name = "Mammalia")
  out <- refdb_filter(toy_ref_db(), taxa = tx)
  expect_equal(unique(out$target_species), "mmusculus")

  out2 <- refdb_filter(toy_ref_db(), targets = "drerio")
  expect_equal(unique(out2$target_species), "drerio")
})

test_that("refdb related operations error when neither targets nor taxa provided", {
  expect_error(refdb_load(), "requires taxonomic scoping")
  expect_error(refdb_filter(toy_ref_db()), "requires taxonomic scoping")
  expect_error(ortholog_resolve("GENE1", ortholog_db = toy_ortholog_db()), "requires taxonomic scoping")
})
