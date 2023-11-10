#!/usr/bin/R
library(tidyverse)
library(qvalue)
library(org.Hs.eg.db)
source("scripts/eQTL/utils.R")

args <- commandArgs(trailingOnly = TRUE)
qtlData <- readQTLtools_perm(args[1])
outputFile <- args[2]

# Need to get gene names

# Filter for non-NA and add column with corrected q-values
qtlData <- qtlData[!is.na(qtlData$adj_beta_pval),] |> 
  mutate("qval" = qvalue(adj_beta_pval)$qvalue)


# Get gene names from orgDb
gene_names <- AnnotationDbi::select(org.Hs.eg.db, keys = qtlData$gene_id,
                      columns = c("ENSEMBL", "SYMBOL"),
                      keytype = "ENSEMBL") |> 
  dplyr::rename(gene_id = ENSEMBL) |> 
  # Collapse gene_id with many symbol mappings
  group_by(gene_id) |> 
  summarise(gene_symbol = paste(SYMBOL, collapse = ",")) |> 
  ungroup()

# Join and write to output
qtlData |> 
  left_join(gene_names, by = "gene_id") |> 
  relocate(gene_symbol, .after = gene_id) |> 
  # Write to file
  write_csv(file = outputFile)