#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

biom_file <- args[1]
out_rds   <- args[2]

# ----------------------------
# RENAME TAXONOMIC RANKS ONLY
# ----------------------------
rename_taxonomic_ranks <- function(physeq) {

  expected_ranks <- c(
    "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species"
  )

  tax <- as.data.frame(tax_table(physeq))

  # Only rename if dimensions match expectation
  if (ncol(tax) == length(expected_ranks)) {
    colnames(tax) <- expected_ranks
  } else {
    message("Warning: taxonomy column count mismatch. Skipping strict rename.")
  }

  tax_table(physeq) <- as.matrix(tax)

  physeq
}

# ----------------------------
# MAIN
# ----------------------------
message("Importing BIOM: ", biom_file)

physeq <- import_biom(biom_file)

message("Renaming taxonomy ranks...")
physeq <- rename_taxonomic_ranks(physeq)

message("Saving RDS: ", out_rds)
saveRDS(physeq, out_rds)

message("Done.")