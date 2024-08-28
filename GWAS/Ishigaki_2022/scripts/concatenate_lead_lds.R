library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

RA <- args[1]
ancestry <- args[2]

path <- paste0(RA, "/", ancestry, "/leads")
ld_files <- list.files(path,
                       pattern = "*_ld.csv")

ld_data <- list()

for (file in ld_files){
  
  data <- read_csv(paste0(path, "/", file), 
                   col_types = "cccdccddddcd")
  ld_data[[file]] <- data
}

ld_data <- ld_data |> bind_rows()

original_leads <- read_csv(paste0(RA, "/", ancestry, "/leads/", RA, "_",
                                  ancestry, 
                                  "_leads.csv"),
                            col_types = "cccdccdddd")

if (nrow(ld_data) > 0){
  full_data <- full_join(ld_data, original_leads, by = colnames(original_leads)) |> 
    separate_wider_delim(cols = "ldbuddy_variantID", delim = ":",
                         names = c(NA, "ldbuddy_pos", NA, NA), 
                         cols_remove = FALSE, too_many = "merge")
} else {
  full_data <- original_leads |> 
    mutate(ldbuddy_variantID = NA,
           ldbuddy_R2 = NA)
}



# Write to file
write_csv(full_data, file = paste0(RA, "/", ancestry, "/leads/",
                                   RA, "_", ancestry, "_leads_ld.csv"))
