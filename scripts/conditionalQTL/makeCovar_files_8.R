library(tidyverse)
library(vcfR)


args <- commandArgs(trailingOnly = TRUE)

cond_topVar_file <- args[1]
covar_file <- args[2]
vcf <- args[3]
file_prefix <- args[4]

# Function for iterating through all eGenes and all combinations of 
# top variants to make covariate files
make_covars <- function(gene, multipleSignal_phenos, covariates, geno_data,
                        file_prefix){
 
  # Filter multipleSignal_phenos to get all the top variants for the gene
  gene_top_variants <- multipleSignal_phenos |> 
    filter(gene_id == gene) |> 
    pull(variantID)
  
  # Go through each variant and create covariate files conditioning on all other
  # lead
  for (var in gene_top_variants){
    conditioning_variants <- gene_top_variants[gene_top_variants != var]
    
    # Pull geno data for conditioning variants
    conditioning_variants_geno <- geno_data |> 
      filter(ID %in% conditioning_variants) |> 
      pivot_wider(names_from = ID, values_from = gt_GT)
    
    
    # Add geno data to covariates
    covariates_conditioning_variants <- covariates |> 
      left_join(conditioning_variants_geno, by = "Indiv") |> 
      column_to_rownames(var = "Indiv") |> 
      t() |> 
      as.data.frame() |> 
      rownames_to_column(var = "covariate")
    
    # Write to file indicating which eGene it's for and which 
    # variant was not conditioned on (i.e. what signal it is for)
    write_delim(covariates_conditioning_variants,
                file = paste0("output/covar/",
                              file_prefix, 
                              "_", gene, "_",
                              var, ".txt"),
                delim = "\t")
    
  }
}



# Read in top signal variants and subset for eGenes with more than 1 signal
multiSignal_phenos <- read_csv(cond_topVar_file) |> 
  group_by(gene_id) |> 
  mutate(n_signal = dplyr::n()) |> 
  filter(n_signal > 1)

# Read in original covariate file and reformat for joining more covariates
covariates <- read_delim(covar_file) |> 
  column_to_rownames(var = "covariate") |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column(var = "Indiv")

# Read in vcf 
read_vcf <- vcfR2tidy(read.vcfR(vcf))
geno_data <- read_vcf$fix |> 
  dplyr::select(POS, ID) |> 
  left_join(read_vcf$gt |> 
              dplyr::select(POS, Indiv, gt_GT),
            by = "POS", relationship = "many-to-many") |> 
  dplyr::select(-POS) 

# Make covariate files for every eGene and combination of conditioned on top variants
invisible(lapply(unique(multiSignal_phenos$gene_id),
       make_covars,
       multipleSignal_phenos = multiSignal_phenos,
       covariates = covariates,
       geno_data = geno_data,
       file_prefix = file_prefix))

