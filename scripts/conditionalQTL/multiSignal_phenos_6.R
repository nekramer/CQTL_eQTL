library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

cond_nSignal_file <- args[1]

# Function to write the phenotype name (gene_id) to a file named itself
write_phenoName_toFile <- function(gene_id){
  # Create file connection with file named after gene_id
  fileConn <- file(paste0("output/cond/", gene_id, ".txt"))
  # Write gene_id
  writeLines(gene_id, fileConn)
  # Close file connection
  close(fileConn)
}


# Read in number of signals per eGene and filter for ones 
# that have more than one signal
cond_nSignals_geneIDs <- read_csv(cond_nSignal_file) |> 
  filter(n_signal > 1) |> 
  pull(gene_id)

# Write each gene_id to it's own text file for input into QTLtools
invisible(lapply(cond_nSignals_geneIDs, write_phenoName_toFile))

