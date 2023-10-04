#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
correction <- args[3] # Which correction to use (FDR or qval both have columns)
threshold <- as.numeric(args[4])



CTL_sig_eGenes <- read_csv(args[1]) |> 
  filter(.[[correction]] <= threshold) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand,
                variantID, variant_chr, variant_start, variant_end, beta, rsID,
                contains(correction))

FNF_sig_eGenes <- read_csv(args[2]) |> 
  filter(.[[correction]] <= threshold) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand,
                variantID, variant_chr, variant_start, variant_end, beta, rsID,
                contains(correction))


FNF_only <- FNF_sig_eGenes[!FNF_sig_eGenes$gene_id |> CTL_sig_eGenes$gene_id,]
write_csv(FNF_only, file = "output/reQTL/FNFonly_sig_eGenes.csv")