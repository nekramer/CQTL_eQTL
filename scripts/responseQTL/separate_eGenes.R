#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
correction <- args[3] # Which correction to use (FDR or qval both have columns)
threshold <- as.numeric(args[4])
fnf_sigOut <- args[5]
ctl_sigOut <- args[6]


CTL_sig_eGenes <- read_csv(args[1]) |> 
  filter(qval < threshold) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand,
                variantID, variant_chr, variant_start, variant_end, beta, rsID,
                contains(correction))

FNF_sig_eGenes <- read_csv(args[2]) |> 
  filter(qval < threshold) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand,
                variantID, variant_chr, variant_start, variant_end, beta, rsID,
                contains(correction))

write_csv(FNF_sig_eGenes, fnf_sigOut)
write_csv(CTL_sig_eGenes, ctl_sigOut)