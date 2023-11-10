library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
OA <- args[1]
chr <- args[2]

# Read in original summary stats
summstats <- read_csv(paste0(OA, "/summary_stats/", OA, "_chr", chr, ".csv"),
                      col_types = "cddccddddd")

# Read in liftOver data 
summ_liftover <- fread(paste0(OA, "/summary_stats/", OA, "_chr", chr, "_liftOver.bed"), 
                        select = c(2, 4), 
                        data.table = FALSE,
                        col.names = c("hg38pos", "CHR:hg19POS"))

# Join together
summstats_liftover <- left_join(summstats, summ_liftover, by = "CHR:hg19POS") |>
  mutate(`CHR:hg38POS` = paste0(chrom, ":", hg38pos)) |>
  relocate(`CHR:hg38POS`, .after = `CHR:hg19POS`) |>
  relocate(hg38pos, .after = hg19pos)

write_csv(summstats_liftover, 
          file = paste0(OA, "/summary_stats/", OA, "_chr", chr, "_liftOver.csv"))