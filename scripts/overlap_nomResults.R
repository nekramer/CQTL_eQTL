library(tidyverse)
library(data.table)
source("scripts/utils.R")


CTL_nom <- readQTLtools_nom("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt")
FNF_nom <- readQTLtools_nom("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt")

both_nom <- inner_join(CTL_nom, FNF_nom, by = c("gene_id", "variantID"), suffix = c("_CTL", "_FNF"))

write_csv(both_nom, 
          file = "/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/nom_overlap.csv")