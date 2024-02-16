#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

cond_topVars_file <- args[1]
permFile <- args[2]
missing_rsID_resultFile <- args[3]
outFile <- args[4]

# Read in variantIDs and rsIDs from perm results to join back
perm_data <- read_csv(permFile, col_select = c("variantID", "rsID"))

# Read in missing rsID results
missing_rsID_results <- read_csv(missing_rsID_resultFile, col_select = c("variantID", "rsID"))

# Join all rsID results together, making sure there are distinct results
rsIDs <- bind_rows(perm_data, missing_rsID_results) |> distinct()

# Read in top signal variants and join rsIDs
cond_topVars_rsIDs <- read_csv(cond_topVars_file) |> 
  left_join(rsIDs, by = "variantID") |> 
  relocate(rsID, .after = "variantID")

# Write to file
write_csv(cond_topVars_rsIDs, file = outFile)

