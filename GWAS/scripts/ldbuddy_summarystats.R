library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

OA <- args[1]
subset <- args[2]

# Read in leads with ld buddies
leads_ld <- read_csv(paste0(OA, "/leads/", 
                            subset, "_", OA, "_leads_ld_liftOver.csv"),
                     col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdccc")

all_leads_ld <- list()

for (chr in unique(leads_ld$chrom)){
  
  # Subset for chrom
  leads_ld_chr <- leads_ld |> 
    filter(chrom == chr)
  
  # Read in summary stats for the chrom
  chr_gwas <- fread(paste0(OA, "/summary_stats/", 
                           OA, "_chr", chr, ".csv"), data.table = FALSE,
                    select = c("CHR:hg19POS", "p"))
  
  # Join with ld_buddies to get pvals
  leads_ld_gwas <- left_join(leads_ld_chr, chr_gwas, 
                             by = join_by(`ldbuddy_CHR:hg19POS` == `CHR:hg19POS`))
  
  all_leads_ld[[chr]] <- leads_ld_gwas
  
}

final_leads_ld <- bind_rows(all_leads_ld) |> 
  filter((!is.na(p)) | (is.na(`ldbuddy_CHR:hg19POS`)))

write_csv(final_leads_ld, 
          file = paste0(OA, "/leads/", subset, "_", OA, "_leads_ld_final.csv"))