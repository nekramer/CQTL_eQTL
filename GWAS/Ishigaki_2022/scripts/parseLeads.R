library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

leadFile <- read_csv(args[1],
                     col_types = "dccdcccddddddcdddddcdddddcdddddcdddddddddddddcdddddcdddddcdddddcd") |> 
  # Remove commas from position and convert to numeric
  mutate(Position = gsub(",", "", Position)) |> 
  mutate(Position = as.numeric(Position))


## Possible lead pop/serostatus combos
lead_types <- leadFile |> 
  dplyr::select(Signal_Population, Signal_Serostatus) |> 
  distinct()

## Go through each lead type combo and subset its leads
for (i in 1:nrow(lead_types)){
  pop <- lead_types[[i, "Signal_Population"]]
  sero <- lead_types[[i, "Signal_Serostatus"]]
  
  leads_pop_sero <- leadFile |> 
    filter(Signal_Population == pop & 
             Signal_Serostatus == sero) |> 
    # Filter out extra columns
    dplyr::select(-`Locus ID`, -Signal_Population,
                  -Signal_Serostatus) |> 
    # Rename columns
    dplyr::rename(rsID = `Rs ID`,
                  variantID = `Variant ID`,
                  chrom = `Chr.`,
                  hg19pos = Position) |> 
    # Add chr prefix to chrom
    mutate(chrom = paste0("chr", chrom)) |> 
    # Split variantID to get alleles; second allele listed is alternate, 
    # which is the effect allele in this study
    separate_wider_delim(cols = "variantID",
                         delim = ":",
                         names = c(NA, NA, "NEA", "EA"),
                         cols_remove = FALSE)
  
  
  if (pop == "Multi"){
    leads_pop_sero <- leads_pop_sero |> 
      dplyr::select(rsID,
                    variantID,
                    chrom,
                    hg19pos, 
                    EA,
                    NEA,
                    paste0("OR_", tolower(pop)),
                    paste0("L95_", tolower(pop)),
                    paste0("U95_", tolower(pop)),
                    paste0("Pvalue_", tolower(pop)))
  } else {
    leads_pop_sero <- leads_pop_sero |> 
      dplyr::select(rsID,
                    variantID,
                    chrom,
                    hg19pos, 
                    EA,
                    NEA,
                    paste0("OR_", pop),
                    paste0("L95_", pop),
                    paste0("U95_", pop),
                    paste0("Pvalue_", pop))
  }  
  
    
  if (pop == "Multi"){
    pop <- "multi_ancestry"
  }
  
  if (sero == "seroposi"){
    write_csv(leads_pop_sero, paste0("SeropositiveRA/", pop,
                                     "/leads/SeropositiveRA_", pop, "_leads.csv"))
  } else {
    write_csv(leads_pop_sero, paste0("AllRA/", pop,
                                     "/leads/AllRA_", pop, "_leads.csv"))
  }
  
  
}


