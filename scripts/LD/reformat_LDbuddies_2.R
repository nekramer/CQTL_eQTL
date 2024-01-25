library(tidyverse)
library(data.table)

## This script takes in a file containing eQTL results and a file output from PLINK's --r2 command
## containing LD buddies for a given SNP and reformats them and joins them together. The output
## is a csv file.

args <- commandArgs(trailingOnly = TRUE)

# The variant ID of the SNP 
snp <- args[1]
# Logical indicating whether or not the variant ID in the LD buddy file has a chr prefix (this will have determined what the variant ID in the plink file is)
ld_chr_prefix <- args[2]
# The original eQTL results file containing the SNP
eQTL_file <- args[3]
# The column name of the variant column in the eQTL file
eQTL_col <- args[4]
# The output file to write to
outFile <- args[5]

# Read in original eQTL results
eQTL_data <- read_csv(eQTL_file)

# Add chr prefix back to snp if necessary for subsetting original eQTL data 
# (i.e. if ld_chr_prefix is false, then we removed it to match the LD reference)
if (ld_chr_prefix == "FALSE"){
  snp_og <- paste0("chr", snp)
} else {
  snp_og <- snp
}
# Subsetting the eQTL results for the specific SNP by its original ID
eQTL_data_snp <- eQTL_data[which(eQTL_data[, eQTL_col] == snp_og),]

## Check for .ld file, if it exists, read it in and join it to the eQTL results
if (file.exists(paste0("output/ld/", snp, ".ld"))){
  # Grab only snp IDs and R2 values from the ld file and rename to match eQTL results
  ldbuddyData <- fread(paste0("output/ld/", snp, ".ld"), data.table = FALSE) |> 
    dplyr::select(SNP_A, SNP_B, R2)

  # Only join if data is not empty
  if (nrow(ldbuddyData) > 0){
    # Add chr prefix back to SNP_A and SNP_B ids in ld buddy file if necessary
    if (ld_chr_prefix == "FALSE"){
      ldbuddyData$SNP_A <- paste0("chr", ldbuddyData$SNP_A)
      ldbuddyData$SNP_B <- paste0("chr", ldbuddyData$SNP_B)
    }
    
    # Rename columns to match
    ldbuddyData <- ldbuddyData |>  
      dplyr::rename(!!eQTL_col := SNP_A,
                    ld_variantID = SNP_B)
    # Join
    snp_eQTL_ld <- left_join(eQTL_data_snp, ldbuddyData, by = eQTL_col)
  } else {
    # Otherwise, do same thing if ld file didn't exist: add snp as an ld buddy to itself with R2 = 1
    snp_eQTL_ld <- eQTL_data_snp |> 
      mutate(ld_variantID = snp_og,
             R2 = 1)
  }
  
} else {
  ## If .ld file doesn't exist, then just add the snp as an ld buddy to itself with R2 = 1
  snp_eQTL_ld <- eQTL_data_snp |> 
    mutate(ld_variantID = snp_og,
           R2 = 1)
}

# Write to file
write_csv(snp_eQTL_ld, outFile)
