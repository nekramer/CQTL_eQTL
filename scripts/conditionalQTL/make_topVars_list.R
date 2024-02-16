library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)


cond_topVar_file <- args[1]
outFile <- args[2]


read_csv(cond_topVar_file, col_select = "variantID") |> 
  distinct() |> 
  write_delim(file = outFile, col_names = FALSE)
  