library(tidyverse)
library(data.table)
source("scripts/utils.R")

ctl_betas <- list()
fnf_betas <- list()

for (chr in 1:22){
  print(chr)
  CTL_nom <- readQTLtools_nom(paste0("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final_chr",
                                           chr, ".txt")) %>%
    mutate("eGene_eSNP" = paste0(gene_id, ":", variantID))
  FNF_nom <- readQTLtools_nom(paste0("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final_chr",
                                     chr, ".txt")) %>%
    mutate("eGene_eSNP" = paste0(gene_id, ":", variantID))
  
  
  # get whichever eGene-eSNP pairs were only nominally significant in one condition
  CTL_nom_only <- CTL_nom %>% filter(!eGene_eSNP %in% FNF_nom$eGene_eSNP) %>%
    dplyr::rename(beta_CTL = beta)
  FNF_nom_only <- FNF_nom %>% filter(!eGene_eSNP %in% CTL_nom$eGene_eSNP) %>%
    dplyr::rename(beta_FNF = beta)
  
  
  # Get corresponding beta for opposite condition
  ctl_opposite_betas <- fread(paste0("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_chr", chr, ".txt"), data.table = FALSE) %>%
    mutate("eGene_eSNP" = paste0(gene_id, ":", variantID)) %>%
    filter(eGene_eSNP %in% CTL_nom_only$eGene_eSNP) %>%
    dplyr::select(eGene_eSNP, beta) %>%
    dplyr::rename(beta_FNF = beta)
  
  fnf_opposite_betas <- fread(paste0("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_chr", chr, ".txt"), data.table = FALSE) %>%
    mutate("eGene_eSNP" = paste0(gene_id, ":", variantID)) %>%
    filter(eGene_eSNP %in% FNF_nom_only$eGene_eSNP) %>%
    dplyr::select(eGene_eSNP, beta) %>%
    dplyr::rename(beta_CTL = beta)
  
  # Join back together
  ctl_betas_all <- left_join(CTL_nom_only, ctl_opposite_betas, by = "eGene_eSNP")
  fnf_betas_all <- left_join(FNF_nom_only, fnf_opposite_betas, by = "eGene_eSNP")
  
  # Add to list to join all chromosomes
  ctl_betas[[chr]] <- ctl_betas_all
  fnf_betas[[chr]] <- fnf_betas_all
}

# Join chromosome subsets
ctl_nomsig_allbetas <- bind_rows(ctl_betas)
fnf_nomsig_allbetas <- bind_rows(fnf_betas)


write_csv(ctl_nomsig_allbetas, 
          file = "/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/CTL_nomsigonly_allbetas.csv")
write_csv(fnf_nomsig_allbetas, 
          file = "/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/FNF_nomsigonly_allbetas.csv")