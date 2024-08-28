library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

lead <- args[1]
RA <- args[2]
ancestry <- args[3]

buddies <- fread(args[4], data.table = FALSE, 
                 select = c(3, 6, 7),
                 col.names = c("variantID", "ldbuddy_variantID", "ldbuddy_R2"))

# Read in set of leads and subset for input lead
lead_data <- read_csv(paste0(RA, "/", ancestry, "/leads/", RA, "_", ancestry, 
                             "_leads.csv")) |>
  filter(variantID == lead)

lead_ldbuddies <- right_join(lead_data, buddies, by = "variantID")


write_csv(lead_ldbuddies, file = paste0(RA, "/", ancestry, "/leads/", RA, "_", ancestry,
                                        "_", lead, "_ld.csv"))