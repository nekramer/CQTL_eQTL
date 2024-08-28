library(data.table)

args <- commandArgs(trailingOnly = TRUE)


RA <- args[1]
ancestry <- args[2]

bed <- fread(paste0(RA, "/", ancestry, "/summary_stats/",
                    RA, "_", ancestry, "_temp.bed"), data.table = FALSE)


fwrite(bed, file = paste0(RA, "/", ancestry, "/summary_stats/",
                          RA, "_", ancestry, "_temp_noheader.bed"),
       quote = FALSE, row.names = FALSE, col.names = FALSE,
       sep = "\t")