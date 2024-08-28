library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)


RA <- args[1]
ancestry <- args[2]

# Read in leads/ld buddies
leads_ld <- fread(paste0(RA, "/", ancestry, "/leads/",
                         RA, "_", ancestry, "_leads_ld_hg38.csv"),
                  data.table = FALSE)

# Read in summary stats
summary_stats <- fread(paste0(RA, "/", ancestry, "/summary_stats/",
                              RA, "_", ancestry, "_hg38.csv"), 
                       data.table = FALSE, select = c("SNP",
                                                      "NEA", 
                                                      "EA", 
                                                      "Beta", 
                                                      "SE", 
                                                      "Pval")) |> 
  dplyr::rename(ldbuddy_NEA = NEA,
                ldbuddy_EA = EA,
                ldbuddy_beta = Beta,
                ldbuddy_beta_se = SE,
                ldbuddy_pval = Pval)

allele_freqs <- list()
# Read in allele frequencies
if (ancestry == "multi_ancestry"){
  # If multi-ancestry, read in PLINK .frq files with minor allele freqs
  for (chr in 1:22){
    chr_allele_freq <- fread(paste0("AFs/chr", chr, "_1000G_GRCh37_ALL_AF.frq"), 
                             data.table = FALSE, select = c("SNP", "A1", "MAF")) |>
      dplyr::rename("MA" = "A1") |> 
      separate_wider_delim(cols = "SNP", delim = ":",
                           names = c(NA, NA, "A1", "A2"), cols_remove = FALSE,
                           too_many = "merge") |> 
      relocate(SNP) |> 
      mutate(A2_AF = case_when(MA == A2 ~ MAF,
                               MA == A1 ~ 1 - MAF)) |> 
      dplyr::select(-MA, -MAF)
      
    allele_freqs[[paste0("chr", chr)]] <- chr_allele_freq
  }
} else {
  for (chr in 1:22){
    chr_allele_freq <- fread(paste0("AFs/chr", chr, "_1000G_GRCh37_",
                                    ancestry, "_AF.txt"), 
                             data.table = FALSE,
                             col.names = c("SNP", "A2_AF")) |> 
      separate_wider_delim(cols = "SNP", delim = ":", 
                           names = c(NA, NA, "A1", "A2"), cols_remove = FALSE,
                           too_many = "merge") |> 
      relocate(SNP)
    allele_freqs[[paste0("chr", chr)]] <- chr_allele_freq
  }
  
}
allele_freqs <- bind_rows(allele_freqs)

# Join together
leads_ld_summary_af <- left_join(leads_ld, summary_stats, 
                                 by = join_by(ldbuddy_variantID == SNP)) |> 
  left_join(allele_freqs, by = join_by(ldbuddy_variantID == SNP)) |> 
  mutate(ldbuddy_EAF = case_when(ldbuddy_EA == A2 ~ as.numeric(A2_AF),
                                 ldbuddy_EA == A1 & is.na(A2_AF) ~ 1 - as.numeric(A2_AF))) |> 
  dplyr::select(-A1,-A2,-A2_AF)


fwrite(leads_ld_summary_af, file = paste0(RA, "/", 
                                          ancestry, "/leads/",
                                          RA, "_", ancestry,
                                          "_leads_ld_hg38_final.csv"),
       sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

