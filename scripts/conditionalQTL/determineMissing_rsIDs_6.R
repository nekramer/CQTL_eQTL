#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

cond_topVars_file <- args[1]
permFile <- args[2]
missing_rsID_outFile <- args[3]

# Read in variantIDs and rsIDs from perm results to join back
perm_data <- read_csv(permFile, col_select = c("variantID", "rsID"))

# Read in conditional data with top variants per signals
read_csv(cond_topVars_file) |> 
  # Join with perm_data to get variant rsIDs
  left_join(perm_data, by = "variantID", relationship = "many-to-many") |> 
  # Filter for variants where rsID is NA
  filter(is.na(rsID)) |> 
  # Only select that variantID column for rsID parsing in the next workflow step
  dplyr::select(variantID) |> 
  # Write to file
  write_csv(file = missing_rsID_outFile)
  