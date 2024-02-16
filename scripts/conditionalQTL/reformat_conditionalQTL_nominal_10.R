library(tidyverse)
source("scripts/utils.R")
args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
gene <- args[2]
var <- args[3]
cond_topVar_file <- args[4]
nomThresholdFile <- args[5]
outFile <- args[6]

# Read in top signal variants and subset for gene and var to get signal number
cond_topVar_gene_var <- read_csv(cond_topVar_file) |> 
  filter(gene_id == gene & variantID == var) |> 
  dplyr::select(variantID, signal)

# Read in nominal threshold to get gene nominal threshold
gene_nomThreshold <- read_delim(nomThresholdFile) |> 
  filter(gene_id == gene) |> 
  pull(pval_nominal_threshold)

# Read in data from name_gene_var nominal conditional pass, add column names,
# add columns for the gene's nominal threshold and 0/1 indicator of whether it 
# passes that threshold, and add columns for indicating the top variant and 
# the associated signal it is for
name_gene_var_nomData <- readQTLtools_nom(paste0('output/cond/', name, "_",
                                                 gene, "_", var,
                                                 '_conditionalQTL_nominal.txt')) |> 
  mutate(pval_nominal_threshold = gene_nomThreshold) |> 
  mutate(nom_sig = ifelse(nom_pval < gene_nomThreshold, 1, 0)) |> 
  mutate(signal_top_variantID = var) |> 
  left_join(cond_topVar_gene_var, by = join_by(signal_top_variantID == variantID))

# Write to file
write_csv(name_gene_var_nomData, outFile)

