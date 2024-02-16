#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

condFile <- args[1]
signal_top_out <- args[2]
n_signals_out <- args[3]

# Read in conditional data that has been reformatted and joined with gene symbol
# and rsID info from perm files
conditionalData <- read_csv(condFile)


## Get top variants for each independent signal by filtering if 
## backward_best_hit == 1 and write to file
signal_top_variants <- conditionalData |> 
  filter(backward_best_hit == 1)
  
write_csv(signal_top_variants, file = signal_top_out)
  
## Get number of signals per gene by counting the number of top variants and
## write to file
signal_top_variants |> 
  # Group by gene
  group_by(gene_id, gene_symbol) |> 
  # Count
  summarise(n_signal = dplyr::n()) |> 
  # Output is 2 columns: gene_id and n_signal
  write_csv(file = n_signals_out)