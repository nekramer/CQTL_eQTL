library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

leadFile <- read_csv(args[1],
                     col_types = "ccccccccdddddddcdddddddddccccddddddddddddddddddddddd") |> 
  # Filter just for leads
  filter(`100Signals` == "Lead")

#OA_phenotypes = c('AllOA', 'KneeHipOA', 'KneeOA', 'HipOA', 'TJR', 'THR', 
                  #'TKR', 'SpineOA', 'FingerOA', 'ThumbOA', 'HandOA')
OA_phenotypes <- c('HipOA')

# Go through each OA phenotype to subset its leads
for (OA in OA_phenotypes){
  
  leads_subtype <- leadFile |> 
    filter(PHENO == OA) |> 
    # Filter out extra columns
    dplyr::select(c(-Signal, -`Novel?`, -`100Signals`, -PHENO, -Direction,
                    -`R2 with lead variant`, -EA_summstats, -NEA_summstats)) |> 
    dplyr::rename(rsID = SNV) |> 
    dplyr::rename(metaP = P) |> 
    separate_wider_delim(cols = "CHR:POS", delim = ":", 
                         names = c("chrom", "hg19pos"),
                         cols_remove = FALSE) |>
    dplyr::rename(`CHR:hg19POS` = `CHR:POS`) |> 
    mutate(hg19pos = as.numeric(hg19pos))

  write_csv(leads_subtype, paste0(OA, "/leads/", OA, "_leads.csv"))
  
}