library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
cond_topVar_file <- args[1]
chrom <- args[2]
nomDataDir <- args[3]
nomDataPrefix <- args[4]
outFile <- args[5]

# Read in conditional top signals and subset for genes with only one signal on the input
# chromosome
singleSignal_phenos_chrom <- read_csv(cond_topVar_file) |> 
  group_by(gene_id) |> 
  mutate(n_signal = dplyr::n()) |> 
  filter(n_signal == 1) |> 
  filter(gene_chr == chrom)

# Read in nominal data for that chromosome
nom_results_chrom_singleSignals <- fread(paste0(nomDataDir, "/", 
                                  nomDataPrefix, "_", 
                                  chrom, ".csv"), data.table = FALSE) |> 
  # filter for single signal genes
  filter(gene_id %in% singleSignal_phenos_chrom$gene_id) |> 
  # Add signal_top_variant as the variantID and set the signal as 0
  mutate(signal_top_variantID = variantID,
         signal = 0)

# Write to file
fwrite(nom_results_chrom_singleSignals, file = outFile, quote = FALSE,
       row.names = FALSE)