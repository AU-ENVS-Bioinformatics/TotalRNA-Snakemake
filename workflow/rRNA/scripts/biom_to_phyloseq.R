#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(biomformat)
  library(phyloseq)
  library(dplyr)
})

########################
# ARGS
########################
args <- commandArgs(trailingOnly = TRUE)

biom_file <- args[1]
out_rds   <- args[2]

########################
# FUNCTIONS
########################

clean_taxonomy <- function(physeq) {
  
  message("\n--- Cleaning taxonomy ---")
  
  tax <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)
  
  colnames(tax) <- c(
    "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species"
  )
  
  tax <- tax %>%
    mutate(across(everything(), ~ sub("^[a-z]__+", "", .))) %>%
    mutate(across(everything(), ~ na_if(., ""))) %>%
    mutate(across(everything(), ~ ifelse(. %in% c("uncultured", "unclassified"), NA, .)))
  
  tax_table(physeq) <- as.matrix(tax)
  
  return(physeq)
}

fix_incertae_sedis <- function(physeq) {
  
  message("\n--- Fixing Incertae Sedis ---")
  
  tax <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)
  
  # ✅ FIX: use t(apply()), NOT do.call()
  tax_fixed <- t(apply(tax, 1, function(row) {
    
    row <- as.character(row)
    
    for (i in seq_along(row)) {
      
      if (is.na(row[i]) || row[i] == "Incertae Sedis") {
        
        replacement <- NA
        
        if (i > 1) {
          for (j in seq(i - 1, 1)) {
            parent <- row[j]
            
            if (!is.na(parent) && parent != "Incertae Sedis") {
              replacement <- parent
              break
            }
          }
        }
        
        if (!is.na(replacement)) {
          row[i] <- paste0("Unclassified (", replacement, ")")
        } else {
          row[i] <- "Unclassified"
        }
      }
    }
    
    return(row)
  }))
  
  colnames(tax_fixed) <- colnames(tax)
  
  tax_table(physeq) <- as.matrix(tax_fixed)
  
  return(physeq)
}

add_sample_data <- function(physeq, otu_mat) {
  
  samdat <- data.frame(
    SampleID = colnames(otu_mat),
    row.names = colnames(otu_mat)
  )
  
  physeq <- merge_phyloseq(physeq, sample_data(samdat))
  
  return(physeq)
}

########################
# SAFE FILTER FUNCTION
########################
remove_multicellular <- function(physeq) {
  
  message("\n--- Removing non-microbial taxa ---")
  
  tax <- as.data.frame(tax_table(physeq))
  
  # 1. remove fully unclassified
  keep <- !(tax$Phylum %in% c(
    "Unclassified (Bacteria)",
    "Unclassified (Eukaryota)"
  ))
  
  # 2. unwanted multicellular taxa
  unwanted <- c(
    "Metazoa", "Animalia",
    "Chordata", "Arthropoda", "Mollusca",
    "Nematozoa", "Nematoda", "Annelida",
    "Rotifera", "Platyhelminthes",
    "Vertebrata", "Mammalia", "Aves",
    "Reptilia", "Amphibia", "Insecta",
    "Embryophyta", "Tracheophyta", "Magnoliophyta"
  )
  
  keep <- keep &
    !(tax$Kingdom %in% unwanted |
      tax$Phylum  %in% unwanted |
      tax$Class   %in% unwanted |
      tax$Order   %in% unwanted |
      tax$Family  %in% unwanted |
      tax$Genus   %in% unwanted)
  
  # 3. remove organelles
  keep <- keep &
    !(tax$Order %in% c("Chloroplast", "Mitochondria") |
      tax$Family %in% c("Mitochondria") |
      tax$Class  %in% c("Chloroplast"))
  
  # ensure no NA sneaks in
  keep[is.na(keep)] <- FALSE
  
  physeq <- prune_taxa(keep, physeq)
  
  return(physeq)
}

########################
# MAIN
########################

message("Loading BIOM: ", biom_file)

biom_obj <- read_biom(biom_file)

########################
# EXTRACT DATA
########################
otu_mat  <- as.matrix(biom_data(biom_obj))
obs_meta <- observation_metadata(biom_obj)

########################
# BUILD TAX TABLE
########################
tax_mat <- do.call(cbind, obs_meta)
tax_mat <- as.matrix(tax_mat)

rownames(tax_mat) <- rownames(otu_mat)

colnames(tax_mat) <- c(
  "Kingdom", "Phylum", "Class",
  "Order", "Family", "Genus", "Species"
)

########################
# FIX OTU ORIENTATION
########################
if (nrow(otu_mat) != nrow(tax_mat)) {
  otu_mat <- t(otu_mat)
}

########################
# BUILD PHYLOSEQ
########################
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)

physeq <- phyloseq(OTU, TAX)

########################
# ADD SAMPLE DATA
########################
physeq <- add_sample_data(physeq, otu_mat)

########################
# CLEAN + FIX TAXONOMY
########################
physeq <- clean_taxonomy(physeq)
physeq <- fix_incertae_sedis(physeq)

########################
# FILTER TAXA
########################
physeq <- remove_multicellular(physeq)

########################
# SAVE OUTPUT
########################
message("Saving phyloseq object: ", out_rds)
saveRDS(physeq, out_rds)

message("Done.")