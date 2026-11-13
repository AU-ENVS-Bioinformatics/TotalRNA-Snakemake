library(tidyverse)
library(phyloseq)

tax_levels <- c(
  'Root',           # 0 (This is life itself)
  'Genome',         # 1 (This is for instance 'mitochondria')
  'Domain',         # 2 (This is Bacteria, Archaea or Eucarya)
  'Superkingdom',   # 3
  'Kingdom',        # 4 (This is also called Superphylum)
  'Phylum',         # 5
  'Class',          # 6
  'Order',          # 7
  'Family',         # 8
  'Genus',          # 9
  'Species',        # 10
  'Strain1',        # 11
  'Strain2',        # 12
  'Strain3',        # 13 (There is a problem with too many ranks)
  'Strain4'
)

run <- function(infile_assignments, infile_counts, outfile_otu, outfile_gg, outphyseq) {
  rRNA <- read_tsv(infile_counts) |>
    right_join(
      infile_assignments |>
        read_tsv(col_names = c("ContigID", "taxonomy")) |>
        separate_wider_delim(
          cols = taxonomy, delim = "; ",
          names = tax_levels, too_few = "align_start"
        ),
      by = "ContigID"
    )

  empty_ranks <- map(rRNA, \(x) sum(is.na(x)) == length(x)) |> keep(\(x) x) |> names()
  rRNA |>
    select(-all_of(empty_ranks)) |>
    write_tsv(file = outfile_otu)

  green_genes_ranks <- c(
    'Superkingdom',
    'Phylum',       # 5
    'Class',        # 6
    'Order',        # 7
    'Family',       # 8
    'Genus',        # 9
    'Species'       # 10
  )

  green_genes_rRNA <- rRNA |>
    filter(Genome == "Main genome") |>
    mutate_at(green_genes_ranks, str_replace_na, replacement = "") |>
    mutate(
      taxonomy = paste0(
        "k__", Superkingdom, "; p__", Phylum,
        "; c__", Class, "; o__", Order,
        "; f__", Family, "; g__", Genus,
        "; s__", Species
      )
    ) |>
    select(-any_of(tax_levels))

  green_genes_rRNA |> write_tsv(outfile_gg)

  otu <- rRNA |>
    filter(Genome == "Main genome") |>
    select(-any_of(tax_levels)) |>
    column_to_rownames("ContigID") |>
    phyloseq::otu_table(taxa_are_rows = TRUE)

  tax <- rRNA |>
    filter(Genome == "Main genome") |>
    select(ContigID, any_of(green_genes_ranks)) |>
    rename("Kingdom" = "Superkingdom") |>
    column_to_rownames("ContigID") |>
    as.matrix() |>
    phyloseq::tax_table()

  phyloseq::phyloseq(otu, tax) |>
    write_rds(file = outphyseq, compress = "bz2")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript process_crest4_phyloseq.R assignments counts out_otu out_gg out_rds")
}

run(args[[1]], args[[2]], args[[3]], args[[4]], args[[5]])