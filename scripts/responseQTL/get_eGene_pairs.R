library(tidyverse)
source('scripts/utils.R')

args <- commandArgs(trailingOnly = TRUE)

eGenes <- read_csv(args[1])
nomData <- readQTLtools_nom(args[2]) %>%
  filter(gene_id %in% eGenes$gene_id)

eGene_snp_pairs <- right_join(eGenes, nomData)

write_csv(eGene_snp_pairs, file = paste0("output/reQTL/",
                                         args[3], "_eGene_snppairs.csv"))