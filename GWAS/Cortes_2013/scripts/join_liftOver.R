library(tidyverse)
library(data.table)

# Read in original lead data
as_leads <- read_csv("IGAS_2013_Table1.csv")

# Read in liftOver
as_leads_liftover <- fread("IGAS_leads_liftOver.bed", 
                           select = c(2, 4),
                           data.table = FALSE,
                           col.names = c("hg38pos", "SNP"))
# Join
as_leads_all <- left_join(as_leads, as_leads_liftover, by = "SNP") |> 
  relocate(hg38pos, .after = "ncbi36_pos") |> 
  # Add variant ID column for matching with 1000G 
  mutate(variantID_v1 = paste0(gsub("chr", "", chr), ":", 
                            hg38pos, ":",
                            risk_allele, ":",
                            non_risk_allele)) |> 
  relocate(variantID_v1, .after = "SNP") |> 
  mutate(variantID_v2 = paste0(gsub("chr", "", chr), ":",
                               hg38pos, ":",
                               non_risk_allele, ":",
                               risk_allele)) |> 
  relocate(variantID_v2, .after = "variantID_v1")

write_csv(as_leads_all, file = "IGAS_2013_leads_hg38.csv")
