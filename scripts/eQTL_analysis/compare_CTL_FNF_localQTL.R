library(tidyverse)
library(dplyr)
library(tidyr)
library(VennDiagram)
source("scripts/utils.R")
library(scales)
library(qvalue)


sig_CTL_eGenes <- read_csv("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05) %>%
  mutate(condition = "Control")

sig_FNF_eGenes <- read_csv("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05) %>%
  mutate(condition = "FN-f")

venn.diagram(list(sig_CTL_eGenes$gene_id, sig_FNF_eGenes$gene_id), filename = "eGene_venndiagram.tiff",
             category.names = c("Control", "FN-f"),
             output = TRUE,
             lwd = 0, 
             fill = c("#7bc5ee", "#fc7971"),
             label.col = "grey15",
             cat.col = c("#7bc5ee", "#fc7971"),
             cat.pos = c(360, 360), 
             cex = 3,
             cat.cex = 2.1,
             rotation.degree = 180)



# Plotting betas ----------------------------------------------------------

overlap <- inner_join(sig_CTL_eGenes, sig_FNF_eGenes, by = c("gene_id", "variantID"),
                      suffix = c("_CTL", "_FNF"))
# 
# overlap_egene_info <- overlap %>% dplyr::select(gene_id, gene_chr_CTL, 
#                                                 gene_start_CTL, gene_end_CTL, 
#                                                 gene_name_CTL, variantID_CTL)
# colnames(overlap_egene_info) <- c("gene_id", "gene_chr", "gene_start", "gene_end", "gene_name")
# 
# write_csv(overlap_egene_info, 
#           file = "/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/both_sigeGene.csv")
# 
# 
# variants_opposite_effect <- overlap %>% filter((beta_CTL < 0 & beta_FNF > 0) | (beta_CTL > 0 & beta_FNF < 0))
# write_csv(variants_opposite_effects[, c("variantID_CTL", "variant_chr_CTL", "variant_start_CTL")],
#           file = "/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/variants_opposite_effects.csv")
# 
# 
ggplot(data = overlap, mapping = aes(x = beta_CTL, y = beta_FNF)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "lightgrey") +
  geom_point(color = "grey25") +
  theme_minimal() +
  xlab("Control effect size") +
  ylab("FN-f effect size")

#  ----------------------------------------------------------
CTL_nom <- readQTLtools_nom("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt") %>%
  mutate("eGene_eSNP" = paste0(gene_id, ":", variantID))
FNF_nom <- readQTLtools_nom("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt") %>%
  mutate("eGene_eSNP" = paste0(gene_id, ":", variantID))



CTL_nom_only <- CTL_nom %>% filter(!eGene_eSNP %in% FNF_nom$eGene_eSNP) %>%
  dplyr::rename(beta_CTL = beta)
FNF_nom_only <- FNF_nom %>% filter(!eGene_eSNP %in% CTL_nom$eGene_eSNP) %>%
  dplyr::rename(beta_FNF = beta)


both_nom <- read_csv("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/nom_overlap.csv")


ggplot(data = both_nom, mapping = aes(x = beta_CTL, y = beta_FNF)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "lightgrey") +
  geom_point(color = "grey25") +
  geom_point(data = CTL_nom_only, mapping = aes(x = bet)) + 
  theme_minimal() +
  xlab("Control effect size") +
  ylab("FN-f effect size")


# TSS distances -----------------------------------------------------------

# Read in nominally significant variants and plot histograms of distance from TSS site
CTL_nom <- readQTLtools_nom("CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt")
FNF_nom <- readQTLtools_nom("FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt")

ggplot(CTL_nom, mapping = aes(x = dist)) +
  geom_histogram(fill = "#7bc5ee", color = "grey25") +
  theme_minimal() +
  ylab("eQTL number") +
  scale_x_continuous(labels = label_number(suffix = "k", scale = 1e-3, big.mark = "")) +
  theme(axis.title.x = element_blank())

ggsave("CTL_nom_TSS.pdf", units = "in", width = 6, height = 5)

ggplot(FNF_nom, mapping = aes(x = dist)) +
  geom_histogram(fill = "#fc7971", color = "grey25") +
  theme_minimal() +
  ylab("eQTL number") +
  scale_x_continuous(labels = label_number(suffix = "k", scale = 1e-3, big.mark = "")) +
  theme(axis.title.x = element_blank())
ggsave("FNF_nom_TSS.pdf", units = "in", width = 6, height = 5)


# pi1 sharing of CTL and FNF ----------------------------------------------


sig_CTL_eGenes <- read_csv("CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05)

sig_FNF_eGenes <- read_csv("FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05)

ourdata_pi1s <- list()

for (n in 1:100){
  
  # CTL tested in FNF
  fnf_ctl_overlap <- sig_FNF_eGenes[which(sig_FNF_eGenes$variantID %in% sig_CTL_eGenes$variantID &
                                            sig_FNF_eGenes$gene_id %in% sig_CTL_eGenes$gene_id),]
  
  # For any not found, sample from uniform distribution
  fnf_ctl_pi1 <- 1- pi0est(c(fnf_ctl_overlap$FDR, 
                             runif(n = length(which(!sig_CTL_eGenes$variantID %in% sig_FNF_eGenes$variantID)))))$pi0
  
  # FNF tested in CTL
  ctl_fnf_overlap <- sig_CTL_eGenes[which(sig_CTL_eGenes$variantID %in% sig_FNF_eGenes$variantID &
                                            sig_CTL_eGenes$gene_id %in% sig_FNF_eGenes$gene_id),]
  
  ctl_fnf_pi1 <- 1- pi0est(c(ctl_fnf_overlap$FDR, 
                             runif(n = length(which(!sig_FNF_eGenes$variantID %in% sig_CTL_eGenes$variantID)))))$pi0
  
  pi1 <- data.frame("combo" = c("Control tested in FN-f", "FN-f tested in Control"),
                    "pi1" = c(fnf_ctl_pi1, ctl_fnf_pi1))
  ourdata_pi1s[[n]] <- pi1
  
}


ourdata_pi1s <- bind_rows(ourdata_pi1s)

ggplot(ourdata_pi1s, mapping = aes(x = combo, y = pi1, fill = combo)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  theme_minimal() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4),
                     labels = c(0.1, 0.2, 0.3, 0.4)) +
  scale_fill_manual(values = c("#B5D9E9", "#FAC3A9")) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none")
ggsave("ctlfnf_pi1comparisons.pdf", units = "in",
       width = 6, height = 5)

# gtex pi1 ----------------------------------------------------------------

gtex_egene_path <- "/proj/phanstiel_lab/External/consortium/gtex/eqtls/egenes/"
gtex_signif_path <- "/proj/phanstiel_lab/External/consortium/gtex/eqtls/signif/"


sig_CTL_eGenes <- read_csv("CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05)
sig_CTL_eGenes$eGene_eSNP <- paste0(sig_CTL_eGenes$gene_id, ":", sig_CTL_eGenes$variantID)

sig_FNF_eGenes <- read_csv("FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05) 
sig_FNF_eGenes$eGene_eSNP <- paste0(sig_FNF_eGenes$gene_id, ":", sig_FNF_eGenes$variantID)

gtex_tissues <- gsub(".v8.egenes.txt.gz", "", list.files(gtex_egene_path))

gtex_eGene_pi1_CTL <- list()

for (tissue in gtex_tissues){

  # Read in GTEx dataset and grab eGene overlaps
  gtex <- read_delim(paste0(gtex_egene_path, tissue, ".v8.egenes.txt.gz"),
                     delim = "\t") %>%
    # Remove decimals from ENS IDs
    mutate(across(gene_id, gsub, pattern = "\\..*", replacement = "")) %>%
    semi_join(sig_CTL_eGenes, by = "gene_id")
  
  
  # For eGene-eSNP pairs not found in gtex, sample from uniform distribution
  # Estimate pi1 from 1 - pi0
  pi1 <- 1 -pi0est(gtex$pval_beta, lambda = seq(0.2, 0.8, 0.01))$pi0
  # pi1 <- 1- pi0est(c(gtex$pval_nominal, 
  #                    runif(n = length(which(!sig_CTL_eGenes$gene_id %in% gtex$gene_id)))),
  #                  lambda = seq(0.2, 0.8, 0.01))$pi0
  
  gtex_eGene_pi1_CTL[[tissue]] <- data.frame("tissue" = tissue,
                                         "pi1" = pi1)
}

gtex_eGene_pi1_CTL %>% bind_rows() %>%
  write_csv(file = paste0("CTL_eGene_GTEx_pi1_padj.csv"))

gtex_eGene_pi1_FNF <- list()

for (tissue in gtex_tissues){
  # Read in GTEx dataset and grab eGene overlaps
  gtex <- read_delim(paste0(gtex_egene_path, tissue, ".v8.egenes.txt.gz"),
                     delim = "\t") %>%
    # Remove decimals from ENS IDs
    mutate(across(gene_id, gsub, pattern = "\\..*", replacement = "")) %>%
    semi_join(sig_FNF_eGenes, by = "gene_id")
  
  
  # For eGene-eSNP pairs not found in gtex, sample from uniform distribution
  # Estimate pi1 from 1 - pi0
  pi1 <- 1 -pi0est(gtex$pval_beta, lambda = seq(0.2, 0.8, 0.01))$pi0
  # pi1 <- 1- pi0est(c(gtex$pval_nominal, 
  #                    runif(n = length(which(!sig_FNF_eGenes$gene_id %in% gtex$gene_id)))),
  #                  lambda = seq(0.2, 0.8, 0.01))$pi0
  
  gtex_eGene_pi1_FNF[[tissue]] <- data.frame("tissue" = tissue,
                                             "pi1" = pi1)
}

gtex_eGene_pi1_FNF %>% bind_rows() %>%
  write_csv(file = paste0("FNF_eGene_GTEx_pi1_updated_adj.csv"))


# Read in FNF, add Condition column, convert tissue order to factor with CTL levels
FNF_eGene_pi1 <- read_csv("FNF_eGene_GTEx_pi1_updated_adj.csv") %>% mutate(Condition = "FNF") %>%
  arrange(pi1) %>%
  # Clean tissue labels
  mutate(across(tissue, gsub, pattern = "_", replacement = " ")) %>%
  mutate(across(tissue, factor, levels = .$tissue))



CTL_eGene_pi1 <- read_csv("CTL_eGene_GTEx_pi1_padj.csv") %>% mutate(Condition = "CTL") %>%
  arrange(pi1) %>%
  # Clean tissue labels
  mutate(across(tissue, gsub, pattern = "_", replacement = " ")) %>%
  mutate(across(tissue, factor, levels = FNF_eGene_pi1$tissue))

# Combine CTL and FNF for plotting
eGene_pi1 <- bind_rows(CTL_eGene_pi1, FNF_eGene_pi1) %>%
  # Convert Condition column to factor
  mutate(across(Condition, as.factor))


# Plot pi1 estimates
ggplot(eGene_pi1, mapping = aes(x = pi1, y = tissue)) +
  geom_point(aes(color = Condition)) +
  theme_minimal() +
  xlim(0, 1) +
  xlab("pi1") +
  #xlab("\u03C0\u2081") +
  scale_color_manual(values = c(brewer.pal(n = 6, "YlGnBu")[3], brewer.pal(n = 6, "YlGnBu")[5]),
                     labels = c("Control", "FN-f")) +
  theme(axis.title.y=element_blank(),
        legend.title = element_blank())

ggsave(filename = "GTEx_eGene_pi1_updated2.pdf", units = "in",
       width = 11, height = 13)


# Steinberg comparison ----------------------------------------------------


low_grade <- read_delim("/proj/phanstiel_lab/External/public/Steinberg_2021/eQTL_LowGrade_FastQTL.txt", 
                        delim = "\t", 
                        col_types = "ccddddddddcdddddddddcc") %>% 
  distinct(Gene, .keep_all = TRUE) %>%
  separate(Gene, into = c("gene_name", "gene_id"), sep = "_")

high_grade <- read_delim("/proj/phanstiel_lab/External/public/Steinberg_2021/eQTL_HighGrade_FastQTL.txt", 
                         delim = "\t", 
                         col_types = "ccddddddddcdddddddddcc") %>%
  distinct(Gene, .keep_all = TRUE) %>%
  separate(Gene, into = c("gene_name", "gene_id"), sep = "_")

sig_CTL_eGenes <- read_csv("CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05)

sig_FNF_eGenes <- read_csv("FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_perm1Mb_FDR.txt") %>%
  filter(FDR < 0.05) 

steinberg_pi1s <- list()

for (n in 1:100){
  
  ctl_low_overlap <- sig_CTL_eGenes[which(sig_CTL_eGenes$gene_id %in% low_grade$gene_id),]
  ctl_lowpi1 <- 1-pi0est(c(runif(length(which(!sig_CTL_eGenes$gene_id %in% low_grade$gene_id))),
                           ctl_low_overlap$nom_pval))$pi0
  
  ctl_high_overlap <- sig_CTL_eGenes[which(sig_CTL_eGenes$gene_id %in% high_grade$gene_id),]
  ctl_highpi1 <- 1-pi0est(c(runif(length(which(!sig_CTL_eGenes$gene_id %in% high_grade$gene_id))),
                            ctl_high_overlap$nom_pval))$pi0
  
  fnf_low_overlap <- sig_FNF_eGenes[which(sig_FNF_eGenes$gene_id %in% low_grade$gene_id),]
  fnf_lowpi1 <- 1-pi0est(c(runif(length(which(!sig_FNF_eGenes$gene_id %in% low_grade$gene_id))),
                           fnf_low_overlap$nom_pval))$pi0
  
  fnf_high_overlap <- sig_FNF_eGenes[which(sig_FNF_eGenes$gene_id %in% high_grade$gene_id),]
  fnf_highpi1 <- 1-pi0est(c(runif(length(which(!sig_FNF_eGenes$gene_id %in% high_grade$gene_id))),
                            fnf_high_overlap$nom_pval))$pi0
  
  pi1 <- data.frame("combo" = c("ctl_low", "ctl_high", "fnf_low", "fnf_high"),
                    "pi1" = c(ctl_lowpi1, ctl_highpi1, fnf_lowpi1, fnf_highpi1))
  steinberg_pi1s[[n]] <- pi1
  
}

steinberg_pi1s <- bind_rows(steinberg_pi1s)
steinberg_pi1s$combo <- factor(steinberg_pi1s$combo,
                               levels = c("ctl_low", "ctl_high",
                                          "fnf_low", "fnf_high"))

ggplot(steinberg_pi1s, mapping = aes(x = combo, y = pi1, fill = combo)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  theme_minimal() +
  scale_y_continuous(breaks = c(0.2, 0.3, 0.4, 0.5),
                     labels = c(0.2, 0.3, 0.4, 0.5)) +
  scale_fill_manual(values = c("#B5D9E9", "#B5D9E9", "#FAC3A9", "#FAC3A9")) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none")
ggsave("steinberg_pi1comparisons.pdf", units = "in",
       width = 6, height = 5)
