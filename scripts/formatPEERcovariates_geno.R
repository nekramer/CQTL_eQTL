#!/usr/bin/R
library(tidyverse)
library(janitor)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

# Read in PEER factors and subset for numPEER
if (file.exists(args[2])){
  print("PEER factor number detected. Reading specific number of PEER factors.")
  numPEER <- read_csv(args[2], show_col_types = FALSE, col_names = FALSE) %>% pull(X1)
  peer <- read_csv(args[1]) %>%
    dplyr::select(paste0("PEER", 1:numPEER), Donor)
} else {
  peer <- read_csv(args[1], show_col_types = FALSE)
  npeer <- ncol(peer)
  print(glue('Reading in {npeer - 1} PEER factors.'))
}

# Read in donorInfo to get sex (age?)
print("Reading in donor samplesheet...")
donorInfo <- read_csv(args[3], show_col_types = FALSE) %>% 
  dplyr::select(Donor, Sex)



# Read in DNA samplesheet
dnaBatches <- c("GenotypingBatch", "DNAReagentBatch")[c(as.logical(args[5]), 
                                              as.logical(args[6]))]
print("Reading in DNA samplesheet...")
dnaInfo <- read_csv(args[4], show_col_types = FALSE) %>% 
  dplyr::select(Donor, all_of(dnaBatches))

# Read in RNA samplesheet
rnaBatches <- c("RNAextractionKitBatch", "SequencingBatch")[c(as.logical(args[8]),
                                                              as.logical(args[9]))]
condition <- args[10]
print("Reading in RNA samplesheet...")
rnaInfo <- read_csv(args[7], show_col_types = FALSE) %>% 
  distinct(Sample, .keep_all = TRUE) %>%
  filter(Condition == condition) %>%
  dplyr::select(Donor, all_of(rnaBatches))

# Join covariates thus far
print("Joining PEER factors with DNA and RNA covariates")
covariates <- peer %>%
  full_join(donorInfo) %>%
  full_join(dnaInfo) %>%
  full_join(rnaInfo) %>%
  # Transpose
  t() %>%
  as.data.frame() %>%
  # Rename rows and columns
  row_to_names(row_number = which(rownames(.) == "Donor"), 
               remove_rows_above = FALSE) %>%
  rownames_to_column(var = "covariate")

# Read in genoPCs
numgenoPC <- as.numeric(args[12])
print(glue('Reading in {numgenoPC} genotyping PCs'))

genoPCs <- read_delim(args[11], delim = " ", show_col_types = FALSE) %>%
  dplyr::select(-SampleID) %>%
  mutate("covariate" = paste0("geno_PC", 1:nrow(.))) %>%
  # Put donor columns in same order of covariates
  relocate(all_of(colnames(covariates))) %>%
  # Filter for numgenoPC
  filter(covariate %in% paste0("geno_PC", 1:numgenoPC))

# Join covariates with genoPCs
print("Adding genotyping PCs to covariates")
covariates <- rbind(covariates, genoPCs)

# Write to tab-delimited file
print("Writing to output")
write_delim(covariates, file = args[13], delim = "\t")
print("Success!")
  