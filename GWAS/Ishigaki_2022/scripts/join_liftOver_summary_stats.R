library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

RA <- args[1]
ancestry <- args[2]

summary_stats <- read_csv(paste0(RA, "/",
                                 ancestry, 
                                 "/summary_stats/",
                                 RA, "_", ancestry, ".csv"))

# Read in liftover data
liftover <- fread(paste0(RA, "/", ancestry, "/summary_stats/",
                               RA, "_", ancestry, "_liftOver.bed"), 
                        select = c(2, 4), 
                        data.table = FALSE,
                        col.names = c("hg38pos", "SNP"))

# Join
summary_liftover <- left_join(summary_stats, liftover, by = "SNP") |>
  relocate(SNP) |> 
  relocate(hg38pos, .after = hg19pos) |> 
  mutate(SNP_hg38 = paste0(gsub("chr", "", chrom),
                           ":",
                           hg38pos,
                           ":",
                           NEA,
                           ":",
                           EA))

write_csv(summary_liftover, file = paste0(RA, "/",
                                          ancestry, "/summary_stats/",
                                          RA, "_", ancestry, "_hg38.csv"))