library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

OA <- args[1]
liftOver_data <- args[2]


# Read in lead file
leads <- read_csv(paste0(OA, "/leads/", OA, "_leads_rsid.csv"),
                  col_types = "cddcccdddddddddddddddddddddddddddddddddddddddddd")

# Read in liftOver data 
leads_liftover <- fread(liftOver_data, select = c(2, 4), 
                        data.table = FALSE,
                        col.names = c("hg38pos", "rsID"))

# Join together
leads_liftover <- left_join(leads, leads_liftover, by = "rsID") |>
    mutate(`CHR:hg38POS` = paste0(chrom, ":", hg38pos)) |>
    relocate(`CHR:hg38POS`, .after = `CHR:hg19POS`) |>
    relocate(hg38pos, .after = hg19pos)


write_csv(leads_liftover, 
        file = paste0(OA, "/leads/", OA, "_leads_liftOver.csv"))
