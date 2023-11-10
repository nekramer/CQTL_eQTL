library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

OA <- args[1]
subset <- args[2]


path <- paste0(OA, "/leads/", subset, "_ld")
ld_files <- list.files(path,
                       pattern = "*_ld.csv")

ld_data <- list()

for (file in ld_files){
  
  data <- read_csv(paste0(path, "/", file), 
                   col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddcdcc")
  ld_data[[file]] <- data
}

ld_data <- ld_data |> bind_rows()

original_leads <- read_csv(paste0(OA, "/leads/", OA, 
                                  "_leads_liftOver_final.csv"),
                           col_types = "cdddccccddddddddddddddd")


full_data <- full_join(ld_data, original_leads, by = colnames(original_leads))

# Write to file
write_csv(full_data, file = paste0(OA, "/leads/", subset, "_", OA, "_leads_ld.csv"))
