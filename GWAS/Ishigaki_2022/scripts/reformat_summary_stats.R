library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)


RA <- args[1]
ancestry <- args[2]

path <- paste0(RA, "/", ancestry, "/summary_stats")

## Read in
summary_stats <- fread(list.files(path, full.name = TRUE), 
                       data.table = FALSE) |> 
  separate_wider_delim(cols = "SNP",
                       delim = "_",
                       names = c("chrom", "hg19pos", "NEA", "EA"),
                       cols_remove = FALSE,
                       too_many = "merge", too_few = "align_start") |> 
  mutate(SNP = gsub("_", ":", SNP),
         chrom = paste0("chr", chrom))

write_csv(summary_stats, file = paste0(RA, "/",
                                       ancestry, "/summary_stats/",
                                       RA, "_", ancestry, ".csv"))