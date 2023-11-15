library(tidyverse)
library(data.table)
library(coloc)

args <- commandArgs(trailingOnly = TRUE)
eGenes <- args[1]
eGenes_ld <- args[2]
gwasPath <- args[3]
gwas_background <- args[4] # ALL or EUR LD buddies from 1000G
nominalPrefix <- args[5]
eqtlN <- as.numeric(args[6])
outFile <- args[7]

# Functions ---------------------------------------------------------------


performColocalization <- function(QTLs, eGene_ld, GWAS_leads, OAtype, 
                                  nominalPrefix, gwasPath){
  
  # Look up variant in GWAS
  qtl_rsID <- QTLs[["rsID"]]
  print(qtl_rsID)
  geneid <- QTLs[["gene_id"]]
  
  # This is looking up in GWAS ld buddies, which includes the leads themselves with R2 = 1
  gwas_variant <- GWAS_leads %>% 
    filter(ldbuddy_rsID == qtl_rsID)
  gwas_R2 <- gwas_variant$ldbuddy_R2
  
  # Look up gwas variant ld in eGene_ld dataset
  eGene_R2 <- eGene_ld |> 
    filter(rsID == qtl_rsID) |> 
    filter(ld_rsID == gwas_variant$rsID) |> 
    pull(R2)
  
  
  # Make sure LD is high enough (trying either GWAS or QTL)
  if (any(gwas_R2 > 0.5, eGene_R2 > 0.5)){
  
    chrom <- gwas_variant$chrom
    print(chrom)
    
    # GWAS --------------------------------------------------------------------
    # Grab gwas sumstats with 100 kb on either side of variant
    gwas_sum_chrom <- fread(paste0(gwasPath,
                                   OAtype,
                                   "/summary_stats/",
                                   OAtype, "_chr",
                                   chrom, ".csv"), data.table = FALSE)
    
    # Get signal region of + or - 250 kb 
    min_region <- gwas_variant$hg38pos - 250000
    max_region <- gwas_variant$hg38pos + 250000
    
    gwas_sum_chrom_subset <- gwas_sum_chrom %>%
      filter(hg38pos >= min_region & hg38pos <= max_region) %>%
      # GWAS gives EAF, so convert to MAF
      mutate(MAF = ifelse(EAF < 0.5, 
                          EAF, 1 - EAF)) |> 
      # Create variantID in same format as QTL variantID
      mutate(variantID = paste0("chr", `CHR:hg38POS`, ":", EA, ":", NEA))
    
    # Get Case/Control number for gwas
    case_control_sizes <- read_csv(paste0(gwasPath, "Case_Control_sampleSizes.csv")) |> 
      filter(OAsubtype == OAtype)
    fraction_cases <- case_control_sizes$Max_Cases/(case_control_sizes$Max_Cases + case_control_sizes$Max_Controls)
    
    gwasN <- case_control_sizes$Max_Cases + case_control_sizes$Max_Controls
    
    # eQTL data ---------------------------------------------------------------
    
    # Pull input nominal QTL data for that region and eGene signal
    qtlData <- read_csv(paste0(nominalPrefix, "_chr", chrom, ".csv")) %>%
      filter(gene_id == geneid) %>%
      filter(variant_start >= min_region & variant_end <= max_region) 
      
    # run coloc --------------------------------------------------------------
     coloc_result <- coloc.abf(dataset1 = list(pvalues = qtlData$nom_pval,
                                               N = eqtlN,
                                               MAF = qtlData$maf,
                                               type = "quant",
                                               beta = qtlData$beta,
                                               snp = qtlData$variantID),
                               dataset2 = list(beta = gwas_sum_chrom_subset$BETA,
                                               s = fraction_cases,
                                               N = gwasN,
                                               type = "cc", 
                                               snp = gwas_sum_chrom_subset$variantID,
                                               MAF = gwas_sum_chrom_subset$MAF,
                                               pvalues = gwas_sum_chrom_subset$p))
  } else {
    coloc_result <- NULL
  }
  
  return(coloc_result)
}



# get significant eGene-snp pairs to  colocalize ----------------------------

eGene_snps <- read_csv(eGenes)

# Read in eGenes with LD information
eGene_ld <- fread(eGenes_ld, data.table = FALSE)

# Iterate through GWAS and perform colocalization -------------------------

OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA",
                "THR", "ThumbOA", "TJR", "TKR")
all_colocs <- list()

for (subtype in OAsubtypes){
  print(paste0("Processing ", subtype))
  
  # Read in corresponding OA subtype leads and their LD buddies
  subtype_leads <- read_csv(paste0(gwasPath, subtype, "/leads/", gwas_background, "_", subtype, 
                                   "_leads_ld_final.csv"),
                            col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd")
  
  # Subset QTLs for ones found in GWAS
  eGene_snps_gwasoverlaps <- eGene_snps %>%
    filter(rsID %in% subtype_leads$ldbuddy_rsID)
  
  if (nrow(eGene_snps_gwasoverlaps) > 0){
    colocalizations <- apply(eGene_snps_gwasoverlaps, 1, performColocalization,
                             eGene_ld = eGene_ld,
                             GWAS_leads = subtype_leads,
                             OAtype = subtype,
                             nominalPrefix = nominalPrefix,
                             gwasPath = gwasPath)
    if (!is.null(colocalizations)){
      names(colocalizations) <- eGene_snps_gwasoverlaps$rsID
    }
    
  } else {
    colocalizations <- NULL
  }
  
  all_colocs[[subtype]] <- colocalizations
}

save(all_colocs, file = outFile)
