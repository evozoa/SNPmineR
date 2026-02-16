test_that("snp_map_ns enforces taxonomic scope", {
  expect_error(snp_map_ns("rs123"), "Provide either")
})

test_that("TaxonSet prints and is classed", {
  x <- taxa_set("Mammalia", species_codes = "mmusculus")
  expect_s3_class(x, "TaxonSet")
})

test_that("ortholog_resolve enforces scoping and filters targets", {
  ortholog_db <- tibble::tibble(
    human_gene_id  = c("ENSG1", "ENSG1", "ENSG2"),
    target_species = c("mmusculus", "rnorvegicus", "mmusculus"),
    target_gene_id = c("MGI1", "RGD1", "MGI2"),
    score          = c(0.9, 0.8, 0.95)
  )

  expect_error(
    ortholog_resolve("ENSG1", ortholog_db = ortholog_db),
    "Provide either `targets`"
  )

  res <- ortholog_resolve(
    human_gene_ids = c("ENSG1", "ENSG2"),
    targets = "mmusculus",
    ortholog_db = ortholog_db
  )

  expect_true(all(res$target_species == "mmusculus"))
  expect_true(all(res$human_gene_id %in% c("ENSG1", "ENSG2")))
})

test_that("ortholog_resolve one_to_many policies behave", {
  ortholog_db <- tibble::tibble(
    human_gene_id  = c("ENSG1", "ENSG1"),
    target_species = c("mmusculus", "mmusculus"),
    target_gene_id = c("MGI_A", "MGI_B"),
    score          = c(0.2, 0.9)
  )

  keep <- ortholog_resolve("ENSG1", targets = "mmusculus", ortholog_db = ortholog_db, one_to_many = "keep_all")
  expect_equal(nrow(keep), 2)

  drop <- ortholog_resolve("ENSG1", targets = "mmusculus", ortholog_db = ortholog_db, one_to_many = "drop_ambiguous")
  expect_equal(nrow(drop), 0)

  best <- ortholog_resolve("ENSG1", targets = "mmusculus", ortholog_db = ortholog_db, one_to_many = "best_hit")
  expect_equal(nrow(best), 1)
  expect_equal(best$target_gene_id, "MGI_B")
})

test_that("validate_variant_db schema is enforced via snp_rsids_to_consequence", {
  variant_db <- tibble::tibble(
    rsid = c("rs1","rs2"),
    gene_id = c("ENSG1","ENSG2"),
    transcript_id = c("ENST1","ENST2"),
    protein_id = c("ENSP1","ENSP2"),
    protein_pos = c(10, 20),
    aa_ref = c("A","G"),
    aa_alt = c("T","D"),
    consequence = c("missense_variant","missense_variant"),
    is_canonical = c(TRUE, FALSE),
    is_protein_coding = c(TRUE, TRUE)
  )

  res <- snp_rsids_to_consequence(c("rs1","rs2"), variant_db, transcript_policy = "canonical")
  expect_true(all(res$rsid == "rs1"))  # only canonical retained
  expect_s3_class(res, "SNPmineHumanConsequence")
})

test_that("snp_rsids_to_consequence returns empty tibble when no matches", {
  variant_db <- tibble::tibble(
    rsid = c("rs1"),
    gene_id = c("ENSG1"),
    transcript_id = c("ENST1"),
    protein_id = c("ENSP1"),
    protein_pos = c(10),
    aa_ref = c("A"),
    aa_alt = c("T"),
    consequence = c("missense_variant")
  )

  res <- snp_rsids_to_consequence("rs999", variant_db)
  expect_equal(nrow(res), 0)
  expect_s3_class(res, "SNPmineHumanConsequence")
})

test_that("snp_map_ns integrates consequence + ortholog resolution", {
  variant_db <- tibble::tibble(
    rsid = c("rs1","rs2"),
    gene_id = c("ENSG1","ENSG2"),
    transcript_id = c("ENST1","ENST2"),
    protein_id = c("ENSP1","ENSP2"),
    protein_pos = c(10, 20),
    aa_ref = c("A","G"),
    aa_alt = c("T","D"),
    consequence = c("missense_variant","missense_variant"),
    is_canonical = c(TRUE, TRUE),
    is_protein_coding = c(TRUE, TRUE)
  )

  ortholog_db <- tibble::tibble(
    human_gene_id  = c("ENSG1","ENSG2"),
    target_species = c("mmusculus","mmusculus"),
    target_gene_id = c("MGI1","MGI2"),
    score          = c(0.9, 0.95)
  )

  res <- snp_map_ns(
    rsids = c("rs1","rs2"),
    targets = "mmusculus",
    variant_db = variant_db,
    ortholog_db = ortholog_db,
    transcript_policy = "canonical",
    one_to_many = "keep_all"
  )

  expect_s3_class(res, "SNPmineResult")
  expect_true(all(res$target_species == "mmusculus"))
  expect_true(all(res$target_gene_id %in% c("MGI1","MGI2")))
})
test_that("coord_map_protein stub returns cached tibble", {
  x1 <- coord_map_protein("ENSP1", "TP1", c(10, 20), cache = TRUE)
  expect_true(is.data.frame(x1))
  expect_true(all(is.na(x1$target_protein_pos)))
  expect_true(all(x1$status == "stub_unimplemented"))

  # second call should hit cache and return same structure
  x2 <- coord_map_protein("ENSP1", "TP1", c(10, 20), cache = TRUE)
  expect_equal(x1, x2)
})

test_that("snp_map_ns map_coords adds coord columns even without target_protein_id", {
  variant_db <- tibble::tibble(
    rsid = c("rs1"),
    gene_id = c("ENSG1"),
    transcript_id = c("ENST1"),
    protein_id = c("ENSP1"),
    protein_pos = c(10),
    aa_ref = c("A"),
    aa_alt = c("T"),
    consequence = c("missense_variant"),
    is_canonical = c(TRUE),
    is_protein_coding = c(TRUE)
  )

  ortholog_db <- tibble::tibble(
    human_gene_id  = c("ENSG1"),
    target_species = c("mmusculus"),
    target_gene_id = c("MGI1"),
    score          = c(0.9)
    # note: no target_protein_id on purpose
  )

  res <- snp_map_ns(
    rsids = "rs1",
    targets = "mmusculus",
    variant_db = variant_db,
    ortholog_db = ortholog_db,
    map_coords = TRUE
  )

  expect_true("target_protein_pos" %in% names(res))
  expect_true("coord_status" %in% names(res))
  expect_true(all(res$coord_status == "stub_no_target_protein_id"))
})
