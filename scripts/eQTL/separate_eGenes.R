#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
permData <- args[1]
threshold <- as.numeric(args[2]) # qval correction threshold
outFile <- args[3]


sig_eGenes <- read_csv(permData) |> 
  filter(qval < threshold) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, gene_symbol, gene_chr, gene_start, gene_end, gene_strand,
                variantID, variant_chr, variant_start, variant_end, beta, qval)


write_csv(sig_eGenes, outFile)
