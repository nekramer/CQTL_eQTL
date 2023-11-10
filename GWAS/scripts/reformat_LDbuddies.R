library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

lead <- args[1]
OA <- args[2]
subset <- args[3] # ALL or EUR

# Get alleles for LD buddies
alleleInfo <- fread(paste0(OA, "/leads/", subset, "_1000G_snps/", lead, "_1000G.bim"),
                    data.table = FALSE,
                    select = c(2, 5, 6),
                    col.names = c("ldbuddy_CHR:hg19POS", "ldbuddy_A1", "ldbuddy_A2")) |> 
  distinct(`ldbuddy_CHR:hg19POS`, .keep_all = TRUE)

buddies <- fread(args[4], data.table = FALSE, 
                 select = c(3,6,7),
                 col.names = c("CHR:hg19POS", "ldbuddy_CHR:hg19POS", "ldbuddy_R2")) |> 
  left_join(alleleInfo, by = "ldbuddy_CHR:hg19POS")


# Read in set of leads and subset for input lead
lead_data <- read_csv(paste0(OA, "/leads/", OA, "_leads_liftOver_final.csv"),
                      col_types = "cdddccccddddddddddddddd") |>
            filter(`CHR:hg19POS` == lead)


lead_ldbuddies <- right_join(lead_data, buddies, by = "CHR:hg19POS")


write_csv(lead_ldbuddies, file = paste0(OA, "/leads/", subset, "_ld/",
                                        lead, "_ld.csv"))




