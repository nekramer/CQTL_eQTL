library(tidyverse)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

selfSM <- read_delim(args[1])
bestSM <- read_delim(args[2])
outfile <- args[3]

seq_sample <- selfSM$`#SEQ_ID`
sample_condition <- str_extract(fileName, "(CTL|FNF)")
best_sample <- bestSM$`#SEQ_ID`

# Check selfSM FREEMIX and CHIPMIX
if (selfSM$FREEMIX > 0.02 & selfSM$CHIPMIX > 0.02){
  
  # Further check for potential swaps that can be corrected within bestSM
  if (bestSM$FREEMIX < 0.02 & bestSM$CHIPMIX < 0.02){
    res <- "swap"
    
  } else {
    
    res <- "check"
  }
  
  
  
} else {
  
  res <- "OK"
  
}

fwrite(data.frame("self" = seq_sample,
                  "condition" = sample_condition,
                  "best" = best_sample,
                  "res" = res),
       file = outfile,
       quote = FALSE,
       row.names = FALSE,
       col.names = FALSE)
