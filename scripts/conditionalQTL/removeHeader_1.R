library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# This input file is the nominal p-value thresholds file
inputFile <- args[1]
outFile <- args[2]

# Read in agnostic of delimiter
nomThresholds <- read_delim(inputFile)

# Write to tab-delimited file without header
write_delim(nomThresholds, outFile, delim = "\t", col_names = FALSE)