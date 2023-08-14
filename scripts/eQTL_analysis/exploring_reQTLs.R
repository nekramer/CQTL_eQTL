library(tidyverse)
library(vcfR)
library(grid)
library(cowplot)
source("scripts/utils.R")
library(coloc)
library(data.table)
library(plotgardener)
#library(gtx)

reformatRSIDdata <- function(dataPath){
  data <- as.data.frame(read_csv(dataPath))
  data[,20] <- gsub("\\[|\\]|'", "", data[,20])
  data <- data %>%
    # Split multiple rsids
    mutate(`0` = str_split(`0`, " ")) %>%
    unnest(`0`) %>%
    # Rename column
    dplyr::rename(rsid = `0`)
  return(data)
}

CTL_sig_reQTLs <- reformatRSIDdata("output/reQTL/CTL_reQTLs_rsids.csv")
FNF_sig_reQTLs <- reformatRSIDdata("output/reQTL/FNF_reQTLs_rsids.csv")

# vcfFile <- "/proj/phanstiel_lab/Data/processed/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_ALL_qc.vcf.gz"
vcfFile <- "output/reQTL/ALL_leadVars.vcf.gz"
#vcfFile <- "/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/vcf/testvariantfilter.recode.vcf.gz"

# Plot the reQTLs ---------------------------------------------------------


CTL_sig_reQTLs <- reformatRSIDdata("output/reQTL/CTL_reQTLs_rsids_old.csv")
FNF_sig_reQTLs <- reformatRSIDdata("output/reQTL/FNF_reQTLs_rsids_old.csv")
# Read in normalized expresion data for CTL and FNF

CTL_normQuant <- read_delim("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch_DNAKitBatch/output/normquant/CTL_CPMadjTMM_invNorm.bed.gz")
FNF_normQuant <- read_delim("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch_DNAKitBatch/output/normquant/FNF_CPMadjTMM_invNorm.bed.gz")


vcf <- vcfR2tidy(read.vcfR(vcfFile, verbose = FALSE))
# Extract genotype matrix and join with info

# Fix sample swaps
geno <- vcf$gt %>% mutate(Indiv = case_when(Indiv == "AM7278" ~ "AM7280", Indiv == "AM7280" ~ "AM7278", TRUE ~ Indiv))
geno_data <- left_join(geno, vcf$fix, by = c("ChromKey", "POS"))


# CTL sig reQTLs
for (var in CTL_sig_reQTLs$variantID){
  var_subset <- geno_data %>% filter(ID == var)
  ref <- unique(var_subset$REF)
  alt <- unique(var_subset$ALT)
  
  eGene_id <- CTL_sig_reQTLs %>% filter(variantID == var) %>% pull(gene_id)
  eGene_name <- CTL_sig_reQTLs %>% filter(variantID == var) %>% pull(gene_name)
  rsid <- CTL_sig_reQTLs %>% filter(variantID == var) %>% pull(rsid)
  pval <- CTL_sig_reQTLs %>% filter(variantID == var) %>% pull(interaction_FDR)
  
  CTL_expression <- CTL_normQuant %>%
    filter(gene_id == eGene_id) %>%
    dplyr::select(-`#chr`, -start, -end, -gene_id, -gene_name, -strand) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Indiv")
  
  CTL_alldata <- left_join(var_subset, CTL_expression) %>%
    # Put genotypes in factor order of ref to alt
    mutate(across(gt_GT_alleles, factor, levels = c(paste0(ref, "/", ref),
                                                    paste0(ref, "/", alt),
                                                    paste0(alt, "/", alt)))) %>%
    group_by(gt_GT_alleles) %>%
    mutate(numGeno = dplyr::n())
  numGenos_CTL <- CTL_alldata %>% dplyr::select(gt_GT_alleles, numGeno) %>% distinct()
  
  # Get corresponding data in FNF
  FNF_expression <- FNF_normQuant %>%
    filter(gene_id == eGene_id) %>%
    dplyr::select(-`#chr`, -start, -end, -gene_id, -gene_name, -strand) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Indiv")
  
  FNF_alldata <- left_join(var_subset, FNF_expression) %>%
    # Put genotypes in factor order of ref to alt
    mutate(across(gt_GT_alleles, factor, levels = c(paste0(ref, "/", ref),
                                                    paste0(ref, "/", alt),
                                                    paste0(alt, "/", alt)))) %>%
    group_by(gt_GT_alleles) %>%
    mutate(numGeno = dplyr::n())
  
  text_geno1 <- textGrob(paste0("n = ", 
                                numGenos_CTL %>% 
                                  filter(gt_GT_alleles == levels(numGenos_CTL$gt_GT_alleles)[1]) %>% 
                                  pull(numGeno)), 
                         gp=gpar(fontsize=12))
  text_geno2 <- textGrob(paste0("n = ", 
                                numGenos_CTL %>% 
                                  filter(gt_GT_alleles == levels(numGenos_CTL$gt_GT_alleles)[2]) %>% 
                                  pull(numGeno)), 
                         gp=gpar(fontsize=12))
  text_geno3 <- textGrob(paste0("n = ", 
                                numGenos_CTL %>% 
                                  filter(gt_GT_alleles == levels(numGenos_CTL$gt_GT_alleles)[3]) %>% 
                                  pull(numGeno)), 
                         gp=gpar(fontsize=12))
  
  # Plot
  ctl <- ggplot(CTL_alldata, aes(x = gt_GT_alleles, y = V1)) + 
    geom_boxplot(color = "#7bc5ee") +
    geom_jitter(position = position_jitter(width = .1),
                color = "#7bc5ee")  +
    coord_cartesian(clip = "off") +
    labs(title = paste0("Control, p = ", pval)) +
    theme_light() +
    ylab(label = paste0(eGene_name, " normalized expression")) +
    theme(axis.title.x = element_text(color = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size =12),
          plot.title = element_text(size = 14, hjust = 0.5)) +
    annotation_custom(text_geno1, 
                      xmin = 1, xmax = 1, ymin=min(CTL_alldata$V1) - 0.625,ymax=min(CTL_alldata$V1) - 0.625) +
    annotation_custom(text_geno2, xmin = 2, 
                      xmax = 2, ymin=min(CTL_alldata$V1) -0.625,ymax=min(CTL_alldata$V1)-0.625) +
    annotation_custom(text_geno3, xmin = 3, 
                      xmax = 3, ymin=min(CTL_alldata$V1)-0.625,ymax=min(CTL_alldata$V1)-0.625)
  
  fnf <- ggplot(FNF_alldata, aes(x = gt_GT_alleles, y = V1)) + 
    geom_boxplot(color = "#fc7971") +
    geom_jitter(position = position_jitter(width = .1),
                color = "#fc7971")  +
    coord_cartesian(clip = "off") +
    labs(title = "FN-f") +
    theme_light() +
    theme(axis.title.x = element_text(color = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title.y = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 12),
          plot.title = element_text(size = 14, hjust = 0.5)) +
    annotation_custom(text_geno1, xmin = 1, 
                      xmax = 1, ymin=min(FNF_alldata$V1)-0.625,ymax=min(FNF_alldata$V1)-0.625) +
    annotation_custom(text_geno2, xmin = 2, 
                      xmax = 2, ymin=min(FNF_alldata$V1)-0.625,ymax=min(FNF_alldata$V1)-0.625) +
    annotation_custom(text_geno3, xmin = 3, 
                      xmax = 3, ymin=min(FNF_alldata$V1)-0.625,ymax=min(FNF_alldata$V1)-0.625)
  
  p <- plot_grid(ctl, fnf)
  title <- ggdraw() + draw_label(rsid, fontface='bold', size = 16)
  
  plot_grid(title, p, ncol = 1, rel_heights=c(0.1, 1))
  ggsave(paste0("output/reQTL/CTL_reQTLs_cleaned/", eGene_name, "_reQTL.pdf"), units = "in",
         width = 8, height = 6)
}

# FNF sig reQTLs
for (var in FNF_sig_reQTLs$variantID){
  var_subset <- geno_data %>% filter(ID == var)
  ref <- unique(var_subset$REF)
  alt <- unique(var_subset$ALT)
  
  eGene_id <- FNF_sig_reQTLs %>% filter(variantID == var) %>% pull(gene_id)
  eGene_name <- FNF_sig_reQTLs %>% filter(variantID == var) %>% pull(gene_name)
  rsid <- FNF_sig_reQTLs %>% filter(variantID == var) %>% pull(rsid)
  pval <- FNF_sig_reQTLs %>% filter(variantID == var) %>% pull(interaction_FDR)
  
  FNF_expression <- FNF_normQuant %>%
    filter(gene_id == eGene_id) %>%
    dplyr::select(-`#chr`, -start, -end, -gene_id, -gene_name, -strand) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Indiv")
  
  FNF_alldata <- left_join(var_subset, FNF_expression) %>%
    # Put genotypes in factor order of ref to alt
    mutate(across(gt_GT_alleles, factor, levels = c(paste0(ref, "/", ref),
                                                    paste0(ref, "/", alt),
                                                    paste0(alt, "/", alt)))) %>%
    group_by(gt_GT_alleles) %>%
    mutate(numGeno = dplyr::n())
  numGenos_FNF <- FNF_alldata %>% dplyr::select(gt_GT_alleles, numGeno) %>% distinct()
  
  # Get corresponding data in CTL
  CTL_expression <- CTL_normQuant %>%
    filter(gene_id == eGene_id) %>%
    dplyr::select(-`#chr`, -start, -end, -gene_id, -gene_name, -strand) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Indiv")
  
  CTL_alldata <- left_join(var_subset, CTL_expression) %>%
    # Put genotypes in factor order of ref to alt
    mutate(across(gt_GT_alleles, factor, levels = c(paste0(ref, "/", ref),
                                                    paste0(ref, "/", alt),
                                                    paste0(alt, "/", alt)))) %>%
    group_by(gt_GT_alleles) %>%
    mutate(numGeno = dplyr::n())
  
  text_geno1 <- textGrob(paste0("n = ", 
                                numGenos_FNF %>% 
                                  filter(gt_GT_alleles == levels(numGenos_FNF$gt_GT_alleles)[1]) %>% 
                                  pull(numGeno)), 
                         gp=gpar(fontsize=12))
  text_geno2 <- textGrob(paste0("n = ", 
                                numGenos_FNF %>% 
                                  filter(gt_GT_alleles == levels(numGenos_FNF$gt_GT_alleles)[2]) %>% 
                                  pull(numGeno)), 
                         gp=gpar(fontsize=12))
  text_geno3 <- textGrob(paste0("n = ", 
                                numGenos_FNF %>% 
                                  filter(gt_GT_alleles == levels(numGenos_FNF$gt_GT_alleles)[3]) %>% 
                                  pull(numGeno)), 
                         gp=gpar(fontsize=12))
  
  # Plot
  fnf <- ggplot(FNF_alldata, aes(x = gt_GT_alleles, y = V1)) + 
    geom_boxplot(color = "#fc7971") +
    geom_jitter(position = position_jitter(width = .1),
                color = "#fc7971")  +
    coord_cartesian(clip = "off") +
    theme_light() +
    labs(title = paste0("FN-f, p = ", pval)) +
    ylab(label = paste0(eGene_name, " normalized expression")) +
    theme(axis.title.x = element_text(color = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size =12),
          plot.title = element_text(size = 14, hjust = 0.5)) +
    annotation_custom(text_geno1, 
                      xmin = 1, xmax = 1, ymin=min(FNF_alldata$V1)-0.625,ymax=min(FNF_alldata$V1)-0.625) +
    annotation_custom(text_geno2, xmin = 2, 
                      xmax = 2, ymin=min(FNF_alldata$V1)-0.625,ymax=min(FNF_alldata$V1)-0.625) +
    annotation_custom(text_geno3, xmin = 3, 
                      xmax = 3, ymin=min(FNF_alldata$V1)-0.625,ymax=min(FNF_alldata$V1)-0.625)
  
  ctl <- ggplot(CTL_alldata, aes(x = gt_GT_alleles, y = V1)) + 
    geom_boxplot(color = "#7bc5ee") +
    geom_jitter(position = position_jitter(width = .1),
                color = "#7bc5ee")  +
    coord_cartesian(clip = "off") +
    theme_light() +
    labs(title = "Control") +
    theme(axis.title.x = element_text(color = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title.y = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 12),
          plot.title = element_text(size = 14, hjust = 0.5)) +
    annotation_custom(text_geno1, xmin = 1, 
                      xmax = 1, ymin=min(CTL_alldata$V1)-0.625,ymax=min(CTL_alldata$V1)-0.625) +
    annotation_custom(text_geno2, xmin = 2, 
                      xmax = 2, ymin=min(CTL_alldata$V1)-0.625,ymax=min(CTL_alldata$V1)-0.625) +
    annotation_custom(text_geno3, xmin = 3, 
                      xmax = 3, ymin=min(CTL_alldata$V1)-0.625,ymax=min(CTL_alldata$V1)-0.625)
  
  p <- plot_grid(fnf, ctl)
  title <- ggdraw() + draw_label(rsid, fontface='bold', size = 16)
  plot_grid(title, p, ncol = 1, rel_heights=c(0.1, 1))

  ggsave(paste0("output/reQTL/FNF_reQTLs_cleaned/", eGene_name, "_reQTL.pdf"), units = "in",
         width = 8, height = 6)
}

# Number of eQTLs ---------------------------------------------------------

# CTL_nom <- readQTLtools_nom("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch_DNAKitBatch/output/qtl/CTL_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt")
# 
# FNF_nom <- readQTLtools_nom("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch_DNAKitBatch/output/qtl/FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb_final.txt")
# 
# CTL_reGenes_qtls <- CTL_nom %>% filter(gene_id %in% CTL_sig_reQTLs$gene_id) %>%
#   group_by(gene_id) %>%
#   summarise(n = dplyr::n()) %>%
#   arrange(desc(n))


CTL_sig_reQTLs <- read_csv("output/reQTL/CTL_sig_reQTLs.csv") %>%
  group_by(gene_id) %>%
  summarise(n = dplyr::n()) %>%
  arrange(desc(n))


ggplot(CTL_sig_reQTLs, mapping = aes(x = n)) +
  geom_histogram(fill = "#7bc5ee", color = "grey25") +
  theme_minimal() +
  xlab("Number of reQTLs") +
  ylab("Number of response eGenes")
ggsave("output/reQTL/CTL_numreQTLs.pdf", units = "in", width = 6, height = 5)

FNF_sig_reQTLs <- read_csv("output/reQTL/FNF_sig_reQTLs.csv") %>%
  group_by(gene_id) %>%
  summarise(n = dplyr::n()) %>%
  arrange(desc(n))

ggplot(FNF_sig_reQTLs, mapping = aes(x = n)) +
  geom_histogram(fill = "#fc7971", color = "grey25") +
  theme_minimal() +
  xlab("Number of reQTLs") +
  ylab("Number of response eGenes")

ggsave("output/reQTL/FNF_numreQTLs.pdf", units = "in", width = 6, height = 5)
# Overlap with differential gene expression -------------------------------

CTL_sig_eGenes <- read_csv("output/reQTL/CTL_sig_reQTLs_old.csv")
FNF_sig_eGenes <- read_csv("output/reQTL/FNF_sig_reQTLs_old.csv")

de_genes <- read_csv("../AllelicImbalance/AIanalysis/data/sig_differential_genes.csv")

CTL_deoverlap <- CTL_sig_eGenes[which(CTL_sig_eGenes$gene_name %in% de_genes$symbol),]
FNF_deoverlap <- FNF_sig_eGenes[which(FNF_sig_eGenes$gene_name %in% de_genes$symbol),]



# Overlap with ASE sites --------------------------------------------------
CTL_sig_reQTLs <- reformatRSIDdata("output/reQTL/CTL_reQTLs_rsids_old.csv")
FNF_sig_reQTLs <- reformatRSIDdata("output/reQTL/FNF_reQTLs_rsids_old.csv")

CTL_qtl_afc <- read_delim("output/reQTL/CTLqtl_afc.csv")
load("../AllelicImbalance/AIanalysis/data/2023-01-09_AIresCTL.rda")
CTL_ase <- as.data.frame(notNA_resCTL) %>%
  rownames_to_column(var = "sid") %>%
  # Mark if significant
  mutate(AIsig = ifelse(padj < 0.05, "TRUE", "FALSE"))

CTL_commonSNPs <- intersect(CTL_ase$sid, CTL_qtl_afc$sid)


CTL_qtl_afc_subset <- CTL_qtl_afc %>% filter(sid %in% CTL_commonSNPs)
CTL_ase_subset <- CTL_ase %>% filter(sid %in% CTL_commonSNPs)

CTL_afc_ase <- left_join(CTL_ase_subset, CTL_qtl_afc_subset)

ggplot(CTL_afc_ase, mapping = aes(x = log2FoldChange, y = log2_aFC)) +
  geom_point(aes(color = AIsig)) +
  theme_light() +
  ylim(c(-4, 4)) +
  xlim(c(-4, 4)) +
  geom_abline(slope = 1, intercept = 0)

FNF_qtl_afc <- read_delim("output/reQTL/FNFqtl_afc.csv")
load("../AllelicImbalance/AIanalysis/data/2023-01-09_AIresFNF.rda")
FNF_ase <- as.data.frame(notNA_resCTL) %>%
  rownames_to_column(var = "sid") %>%
  # Mark if significant
  mutate(AIsig = ifelse(padj < 0.05, "TRUE", "FALSE"))

FNF_commonSNPs <- intersect(FNF_ase$sid, FNF_qtl_afc$sid)


FNF_qtl_afc_subset <- FNF_qtl_afc %>% filter(sid %in% FNF_commonSNPs)
FNF_ase_subset <- FNF_ase %>% filter(sid %in% FNF_commonSNPs)

FNF_afc_ase <- left_join(FNF_ase_subset, FNF_qtl_afc_subset)

ggplot(FNF_afc_ase, mapping = aes(x = log2FoldChange, y = log2_aFC)) +
  geom_point(aes(color = AIsig)) +
  theme_light() +
  ylim(c(-4, 4)) +
  xlim(c(-4, 4)) +
  geom_abline(slope = 1, intercept = 0)


# Colocalization ----------------------------------------------------------

FNF_sig_reQTLs <- read_csv("output/reQTL/FNF_sigreGenes_GRCh37pos.csv") %>%
  mutate(snp = paste0(variant_chr, ":", GRCh37pos))
FNF_sig_reQTLs$snp <- gsub("chr", "", FNF_sig_reQTLs$snp)  


load("output/reQTL/KneeHipOA_FNF_colocalizations.rda")
FNF_colocalizations[sapply(FNF_colocalizations, is.null)] <- NULL
FNF_colocalizations_filtered <- FNF_colocalizations %>% keep(~.$summary["PP.H4.abf"] >= 0.2)

coloc1 <- FNF_colocalizations_filtered[[6]]

eGene <- FNF_sig_reQTLs %>% filter(snp == coloc1$results$snp[1]) %>% pull(gene_id)

FNF_sig_reQTLs_coloc1 <- FNF_sig_reQTLs %>% filter(gene_id == eGene)

min_region <- min(FNF_sig_reQTLs_coloc1$GRCh37pos) - 100000
max_region <- max(FNF_sig_reQTLs_coloc1$GRCh37pos) + 100000

gwas_chr16 <- fread("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/KneeHipOA/summary_statistics/KneeHipOA_chr16.txt.gz", data.table = FALSE)

gwas_subset <- gwas_chr16 %>% filter(pos >= min_region & pos <= max_region) %>% filter(!EffectAlleleFrequency >= 1) %>% distinct(snp, .keep_all = TRUE)


besthit <- paste0("chr", FNF_sig_reQTLs_coloc1 %>% filter(best_hit == 1) %>% pull(snp))

# Got LD results for best hit against 1000G
coloc1_ld <- fread("output/reQTL/FNF_kneehipccoloc_chr16_snp69964374.csv", data.table = FALSE) 
coloc1_ld <- coloc1_ld[,c("SNP_B", "R2")]
colnames(coloc1_ld)[1] <- "snp"

coloc1_gwas_subset_ld <- left_join(gwas_subset, coloc1_ld)



coloc1_gwas_subset_ld <- as.data.frame(dplyr::group_by(coloc1_gwas_subset_ld,
                                                   LDgrp = cut(
                                                     coloc1_gwas_subset_ld$R2,
                                                     c(0, 0.2, 0.4, 0.6, 0.8, 1))))

coloc1_gwas_subset_ld$LDgrp <- addNA(coloc1_gwas_subset_ld$LDgrp)



# FNF_sig_reQTLs_coloc1 <- left_join(FNF_sig_reQTLs_coloc1, coloc1_ld)
# 
# FNF_sig_reQTLs_coloc1 <- as.data.frame(dplyr::group_by(FNF_sig_reQTLs_coloc1,
#                                                        LDgrp = cut(
#                                                          FNF_sig_reQTLs_coloc1$R2,
#                                                          c(0, 0.2, 0.4, 0.6, 0.8, 1))))
# FNF_sig_reQTLs_coloc1$LDgrp <- addNA(FNF_sig_reQTLs_coloc1$LDgrp)
# FNF_sig_reQTLs_coloc1$p <- FNF_sig_reQTLs_coloc1$nom_pval
# FNF_sig_reQTLs_coloc1 <- FNF_sig_reQTLs_coloc1 %>%
#   dplyr::rename(chrom = variant_chr) %>%
#   dplyr::rename(pos = GRCh37pos)



#################### Knee
############ 
load("output/reQTL/KneeOA_FNF_colocalizations.rda")
FNF_colocalizations[sapply(FNF_colocalizations, is.null)] <- NULL
FNF_colocalizations_filtered <- FNF_colocalizations %>% keep(~.$summary["PP.H4.abf"] >= 0.2)

coloc2 <- FNF_colocalizations_filtered[[2]]
eGene2 <- FNF_sig_reQTLs %>% filter(snp == coloc2$results$snp[1]) %>% pull(gene_id)

FNF_sig_reQTLs_coloc2 <- FNF_sig_reQTLs %>% filter(gene_id == eGene2)

min_region <- min(FNF_sig_reQTLs_coloc2$GRCh37pos) - 100000
max_region <- max(FNF_sig_reQTLs_coloc2$GRCh37pos) + 100000

kneegwas_chr16 <- fread("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/KneeOA/summary_statistics/KneeOA_chr16.txt.gz", data.table = FALSE)

kneegwas_subset <- kneegwas_chr16 %>% filter(pos >= min_region & pos <= max_region) %>% filter(!EffectAlleleFrequency >= 1) %>% distinct(snp, .keep_all = TRUE)


besthit2 <- paste0("chr", FNF_sig_reQTLs_coloc2 %>% filter(best_hit == 1) %>% pull(snp))

# Got LD results for best hit against 1000G
coloc2_ld <- fread("output/reQTL/FNF_kneehipccoloc_chr16_snp69964374.csv", data.table = FALSE) 
coloc2_ld <- coloc2_ld[,c("SNP_B", "R2")]
colnames(coloc2_ld)[1] <- "snp"

coloc2_gwas_subset_ld <- left_join(kneegwas_subset, coloc2_ld)
coloc2_gwas_subset_ld <- as.data.frame(dplyr::group_by(coloc2_gwas_subset_ld,
                                                       LDgrp = cut(
                                                         coloc2_gwas_subset_ld$R2,
                                                         c(0, 0.2, 0.4, 0.6, 0.8, 1))))

coloc2_gwas_subset_ld$LDgrp <- addNA(coloc2_gwas_subset_ld$LDgrp)


FNF_sig_reQTLs_coloc2 <- left_join(FNF_sig_reQTLs_coloc2, coloc2_ld)

# Read in nominal results 
FNF_region_nominal_results <- read_csv("output/reQTL/FNF_region_nominal_results_rsids_GRCh37.csv") %>%
 mutate(snp = paste0(variant_chr, ":", GRCh37pos))
FNF_region_nominal_results$snp <- gsub("chr", "", FNF_region_nominal_results$snp)

FNF_region_nominal_results <- left_join(FNF_region_nominal_results, coloc2_ld)
FNF_region_nominal_results$interaction_pval <- NA
FNF_region_nominal_results$interaction_FDR <- NA
FNF_region_nominal_results$gene_name <- NA
FNF_region_nominal_results$FDR <- NA
FNF_region_nominal_results$pval_nominal_threshold <- NULL

# bind it all together
FNF_colocRegion_all <- bind_rows(FNF_sig_reQTLs_coloc2, FNF_region_nominal_results)

# Reformat for plotgardener

FNF_colocRegion_all <- FNF_colocRegion_all %>%
  dplyr::rename(chrom = variant_chr) %>%
  dplyr::rename(pos = GRCh37pos) %>%
  dplyr::rename(p = nom_pval)


FNF_colocRegion_all <- as.data.frame(dplyr::group_by(FNF_colocRegion_all,
                                                       LDgrp = cut(
                                                         FNF_colocRegion_all$R2,
                                                         c(0, 0.2, 0.4, 0.6, 0.8, 1))))

FNF_colocRegion_all$LDgrp <- addNA(FNF_colocRegion_all$LDgrp)

FNF_colocRegion_all <- FNF_colocRegion_all %>% distinct(snp, .keep_all = TRUE)

##### Plot
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

pageCreate(width = 6, height = 5, default.units = "inches")

man1 <- plotManhattan(data = coloc1_gwas_subset_ld, assembly = "hg19",
                      chrom = "chr16", chromstart = min_region, chromend = max_region,
                      fill = colorby("LDgrp",
                                     palette = colorRampPalette(c(
                                       "#262C74",
                                       "#98CDED", "#499A53",
                                       "#EEA741", "#DD3931", "#262C74"
                                     ))), range = c(0, 18),
                      x = 0.5, y = 0, width = 4, height = 1.25, sigLine = TRUE, 
                      sigVal = 1.3e-8, lty = 2,
                      baseline = TRUE, lwd = 0.5,
                      leadSNP = list(snp = "16:69964374",
                                     fill = "#DD3931",
                                     pch = 18, 
                                     cex = 0.75,
                                     fontsize = 8,
                                     fontcolor =NA))
annoSegments(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = -1*log10(2.27e-5), 
             y1 = -1*log10(2.27e-5), plot = man1, default.units = "native")
plotText("KneeHipOA", x = 4.5, y = 0.2,just = c("right", "top"), fontsize = 10,
         fontcolor = "#414141", fontface = "bold")
plotText("PP4 = 0.794", x = 4.5, y = 0.4,just = "right", fontsize = 8,fontcolor = "#414141" )
man2 <- plotManhattan(data = coloc2_gwas_subset_ld, assembly = "hg19",
                      chrom = "chr16", chromstart = min_region, chromend = max_region,
                      fill = colorby("LDgrp",
                                     palette = colorRampPalette(c(
                                       "#262C74",
                                       "#98CDED", "#499A53",
                                       "#EEA741", "#DD3931", "#262C74"
                                     ))), range = c(0, 18),
                      x = 0.5, y = 1.4, width = 4, height = 1.25, sigLine = TRUE,
                      sigVal = 1.3e-8,
                      lty = 2, baseline = TRUE, lwd = 0.5,
                      leadSNP = list(snp = "16:69964374",
                                     fill = "#DD3931",
                                     pch = 18,
                                     cex = 0.75,
                                     fontsize = 8,
                                     fontcolor =NA))
annoSegments(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = -1*log10(2.27e-5), 
             y1 = -1*log10(2.27e-5), plot = man2, default.units = "native")
plotText("KneeOA", x = 4.5, y = 1.6,just = c("right", "top"), fontsize = 10,
         fontcolor = "#414141", fontface = "bold")
plotText("PP4 = 0.943", x = 4.5, y = 1.79,just = "right", fontsize = 8,fontcolor = "#414141" )

man3 <- plotManhattan(data = FNF_colocRegion_all, assembly = "hg19",
                      chrom = "chr16", chromstart = min_region, chromend = max_region,
                      fill = colorby("LDgrp",
                                     palette = colorRampPalette(c(
                                       "#262C74",
                                       "#98CDED", "#499A53",
                                       "#EEA741", "#DD3931", "#262C74"
                                     ))), range = c(0, 18),
                      x = 0.5, y = 2.8, width = 4, height = 1.25, baseline = TRUE,
                      leadSNP = list(snp = "16:69964374",
                                     fill = "#DD3931",
                                     pch = 18,
                                     cex = 0.75,
                                     fontsize = 8,
                                     fontcolor ="#DD3931"))

annoSegments(plot = man3, x0 = unit(0, "npc"), x1 = unit(1, "npc"),
             y0 = -1*log10(0.0000292), y1 = -1*log10(0.0000292))


plotText("reQTL", x = 4.5, y = 3,just = c("right", "top"), fontsize = 10,
         fontcolor = "#414141", fontface = "bold")
annoYaxis(plot = man1, at = seq(0, 18, 5), axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.15, y = 0.7, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches"
)

annoYaxis(plot = man2, at = seq(0, 18, 5), axisLine = TRUE, fontsize = 8)

plotText(
  label = "-log10(p-value)", x = 0.15, y = 2.1, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches"
)
annoYaxis(plot = man3, at = seq(0, 18, 5), axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.15, y = 3.5, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches"
)


gene_track <- plotGenes(chrom = "chr16", chromstart = min_region, chromend = max_region,
          assembly = "hg19", x = 0.5, y = 4.2, width = 4, height = 0.5,
          geneHighlights = data.frame("gene" = "PDXDC2P-NPIPB14P",
                                      "color" = "#669fd9"),
          geneOrder = c("PDXDC2P-NPIPB14P", "WWP2"))


plotGenomeLabel(chrom = "chr16", chromstart = min_region, 
                chromend = max_region, assembly = "hg19",
                x = 0.5, y = 4.75, length = 4, fontsize = 8)

plotLegend(legend = c("0.8 - 0.1",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "No LD information"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 4.45, y = 0.2, width = 0.1, height = 0.5, border = FALSE, 
           fontsize = 8)

pageGuideHide()
