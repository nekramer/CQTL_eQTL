library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

OA <- args[1]
subset <- args[2]
liftOver_data <- args[3]

# Read in lead file with ld buddies
leads_ld <- read_csv(paste0(OA, "/leads/", 
                         subset, "_", OA, 
                         "_leads_ld_rsid.csv"),
                  col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddcdccc")

# Read in liftOver data (ld buddies ID's by chrom:hg19pos)
ld_liftover <- fread(liftOver_data, select = c(2, 4), 
                        data.table = FALSE,
                        col.names = c("ldbuddy_hg38pos", "ldbuddy_rsID")) |> 
  distinct(.keep_all = TRUE)


# Join together
leads_ld_liftover <- left_join(leads_ld, ld_liftover, by = "ldbuddy_rsID") |>
  mutate(`ldbuddy_CHR:hg38POS` = paste0(chrom, ":", ldbuddy_hg38pos)) |>
  relocate(`ldbuddy_CHR:hg38POS`, .after = `ldbuddy_CHR:hg19POS`) |> 
  dplyr::select(-ldbuddy_hg38pos)

write_csv(leads_ld_liftover, 
          file = paste0(OA, "/leads/", subset, "_", OA, "_leads_ld_liftOver.csv"))