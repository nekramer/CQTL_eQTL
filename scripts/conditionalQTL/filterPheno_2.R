#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
phenoFile <- args[1]
nomThresholdFile <- args[2]
outFile <- args[3]

# Read in nominal threshold file, where phenotype ID's are in the first column
nomThresholds <- read_delim(nomThresholdFile, col_names = FALSE)

# Read in phenotype data, where phenotype ID's are in the fourth column
phenoData <- read_delim(phenoFile) |> 
  # Filter for phenotype IDs in nomThresholds
  filter(if_all(4, ~ . %in% nomThresholds[[1]])) |> 
  # Write to file
  write_delim(file = outFile, delim = "\t")


