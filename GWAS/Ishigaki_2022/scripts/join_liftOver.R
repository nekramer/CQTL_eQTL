library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

RA <- args[1]
ancestry <- args[2]

leads_ld <- read_csv(paste0(RA, "/", ancestry, 
                         "/leads/", RA, "_", ancestry, "_leads_ld.csv"))

# Read in lead liftOver data 
leads_liftover <- fread(paste0(RA, "/", ancestry, "/leads/",
                               RA, "_", ancestry, "_lead_liftOver.bed"), 
                        select = c(2, 4), 
                        data.table = FALSE,
                        col.names = c("hg38pos", "rsID")) |> 
  distinct()


# Read in ld buddy liftOver data
buddies_liftover <- fread(paste0(RA, "/", ancestry, "/leads/",
                                 RA, "_", ancestry, "_buddies_liftOver.bed"),
                          select = c(2, 4),
                          data.table = FALSE,
                          col.names = c("ldbuddy_hg38pos", "ldbuddy_variantID"),
                          quote = "") |> 
  # Remove any leading quotes
  mutate(ldbuddy_variantID = gsub('\"', "", ldbuddy_variantID))


# Join together
leads_ld_liftover <- left_join(leads_ld, leads_liftover, by = "rsID") |>
  relocate(hg38pos, .after = hg19pos) |> 
  mutate(variantID_hg38 = paste0(gsub("chr", "", chrom),
                                 ":",
                                 hg38pos,
                                 ":",
                                 NEA,
                                 ":",
                                 EA)) |> 
  relocate(variantID_hg38, .after = variantID) |> 
  left_join(buddies_liftover, by = "ldbuddy_variantID") |> 
  relocate(ldbuddy_hg38pos, .after = ldbuddy_pos)

write_csv(leads_ld_liftover, 
          file = paste0(RA, "/", ancestry, 
                        "/leads/", RA, "_", 
                        ancestry, "_leads_ld_hg38.csv"))
