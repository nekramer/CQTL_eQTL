#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
conditionalFile <- args[1]
permFile <- args[2]
outFile <- args[3]

# Read in conditional data and select columns
conditional_data <- read_delim(conditionalFile, 
                               col_names = FALSE,
                               col_select = c(1:5, 8:12, 18:21))
# Add column names
colnames(conditional_data) <- c("gene_id", "gene_chr",
                                "gene_start", "gene_end",
                                "gene_strand", "variantID",
                                "variant_chr", "variant_start",
                                "variant_end", "signal",
                                "backward_nom_pvalue",
                                "backward_r_squared", "backward_slope",
                                "backward_best_hit")

# Read in geneIDs, gene symbols from perm results to join back
perm_data <- read_csv(permFile, col_select = c("gene_id", "gene_symbol"))

# Filter conditional results for phenos in the perm data
conditional_data_perm_phenos <- conditional_data |> 
  filter(gene_id %in% perm_data$gene_id)

# Join gene symbol info from perm results to conditional results
conditional_data_final <- left_join(conditional_data_perm_phenos,
                                    perm_data,
                                    by = "gene_id") |> 
  relocate(gene_symbol, .after = gene_id)

write_csv(conditional_data_final, file = outFile)