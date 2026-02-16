# SNPmineR

SNPmineR is an R package scaffold for mapping human SNP rsIDs to orthologous nonsynonymous variants in target species.

## One-liner workflow

Use `snp_map_ns()` for a tidyverse-style pipeline (vector-in, tibble-out):

```r
library(SNPmineR)

results <- snp_map_ns(
  rsids = c("rs1", "rs2"),
  targets = c("mmusculus"),
  variant_db = variant_db,
  ortholog_db = ortholog_db,
  ref_db_targets = ref_db_targets
)
```

## Constrain searches with targets or taxa

SNPmineR enforces explicit taxonomic scope for ortholog/refdb access. You must supply either:

- `targets = c("species_code", ...)`, or
- a `TaxonSet`, e.g. clade bundle style:

```r
mammalia <- taxa_set(list(
  name = "Mammalia",
  taxids = c(40674L),
  species_codes = c("mmusculus", "rnorvegicus"),
  source = "user"
))

results <- snp_map_ns(
  rsids = c("rs1"),
  taxa = mammalia,
  variant_db = variant_db,
  ortholog_db = ortholog_db,
  ref_db_targets = ref_db_targets
)
```

If neither `targets` nor `taxa` is provided, SNPmineR errors by design to avoid implicit global taxonomy scans.

## Development note

NEXT STEPS: implement alignment-based coordinate mapping to support high-confidence residue-level transfer between human and target transcripts/proteins.
