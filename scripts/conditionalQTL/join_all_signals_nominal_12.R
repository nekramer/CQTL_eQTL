library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
cond_topVar_file <- args[1]
name <- args[2] 
chrom <- args[3]
chrom_signal1_nomFile <- args[4]
outFile <- args[5]


read_eGene_variant_nomData <- function(multiSignal_geneRow, name){
  
  gene_id <- multiSignal_geneRow[["gene_id"]]
  variant <- multiSignal_geneRow[["variantID"]]
  
  reformatted_nomData <- read_csv(paste0("output/cond/", name, "_",
                                         gene_id, "_",
                                         variant, 
                                         "_conditionalQTL_nominal_reformat.csv"),
                                  col_types = "ccddcddccddddddddcd")
  return(reformatted_nomData)
}

# Read in conditional top signals and subset for genes with more than 1 signal on
# the chromosome
multiSignal_phenos_chrom <- read_csv(cond_topVar_file) |> 
  group_by(gene_id) |> 
  mutate(n_signal = dplyr::n()) |> 
  filter(n_signal > 1) |> 
  filter(gene_chr == chrom)

# Iterate through reading in conditional nominal results for the eGenes
# on this chromosome
conditional_nomData_chrom <- apply(multiSignal_phenos_chrom,
                                   1,
                                   read_eGene_variant_nomData,
                                   name = name) |> 
  bind_rows()

# Read in single signal nominal data and join conditional_nomData_chrom
nomData_allSignals_chrom <- fread(chrom_signal1_nomFile, data.table = FALSE) |> 
  bind_rows(conditional_nomData_chrom)

# Write to file
fwrite(nomData_allSignals_chrom, file = outFile, quote = FALSE,
       row.names = FALSE)
