library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

all_fnf_leads <- read_csv(args[1])
fnf_response <- read_csv(args[2])

fnf_noresponse <- all_fnf_leads |> 
  filter(!gene_id %in% fnf_response$gene_id)


write_csv(fnf_noresponse, file = args[3])