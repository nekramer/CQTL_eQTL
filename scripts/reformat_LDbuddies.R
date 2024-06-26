library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
leadFile <- args[1]
ldbuddyFile <- args[2]
snp <- args[3]
outFile <- args[4]


# Read in lead data, subsetting for that lead
leadData <- read_csv(leadFile) |> 
  filter(variantID == snp)

# Read in ld buddy output from PLINK, grabbing relevant columns
if (file.exists(ldbuddyFile)){
  ldbuddyData <- fread(ldbuddyFile, data.table = FALSE) |> 
    dplyr::select(SNP_A, SNP_B, R2) |> 
    dplyr::rename(variantID = SNP_A,
                ld_variantID = SNP_B)
                # Join ld buddies to lead snp
lead_ld <- left_join(leadData, ldbuddyData, by = "variantID")
} else {
  # Add LD buddy as itself if no LD buddies found
  lead_ld <- leadData |> 
    mutate(ld_variantID = snp,
    R2 = 1)
}

# Write to file
write_csv(lead_ld, outFile)