library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

cond_nSignal_file <- args[1]
cond_topVar_file <- args[1]
perm_LD_file <- args[2]
cond_signal_LD_file <- args[3]
outFile <- args[4]

# Get eGenes that only had 1 signal and get those 
# Read in number of signals per eGene and filter for ones 
# that have more than one signal
cond_1signal_geneIDs <- read_csv(cond_nSignal_file) |> 
  filter(n_signal == 1) |> 
  pull(gene_id)

# Read in conditional top variants per signals and isolate ones that weren't for 
# secondary or more signals
cond_variant_top_signal0 <- read_csv(cond_topVar_file) |> 
  filter(signal == 0)

# Read in original perm LD files and filter for signal 0 variants
signal0_perm_LD <- fread(perm_LD_file, data.table = FALSE) |> 
  filter(variantID %in% cond_variant_top_signal0$variantID) |> 
  dplyr::select(variantID, ld_variantID, R2, ld_rsID)

# Join cond_variant_top_signal0 with LD information
cond_variant_top_signal0_LD <- cond_variant_top_signal0 |> 
  left_join(signal0_perm_LD, by = "variantID", relationship = "many-to-many") |> 
  # Rearrange columns to match original LD file structure
  relocate(rsID, variantID)

# Read in the LD file for variants with secondary or more signals
cond_signal_LD <- fread(cond_signal_LD_file, data.table = FALSE)

# Join together
cond_allsignals_LD <- bind_rows(cond_variant_top_signal0_LD,
                                cond_signal_LD)

# Write to file
fwrite(cond_allsignals_LD, file = outFile, sep = ",")