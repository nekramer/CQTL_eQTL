library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

reQTL_data <- read_csv(args[1])
LD_data <- fread(args[2], data.table = FALSE) |> 
  dplyr::select(c("rsID", "ld_variantID", "R2", "ld_rsID"))


joined_data <- left_join(reQTL_data, LD_data, by = c("rsID"))

write_csv(joined_data, file = args[3])



