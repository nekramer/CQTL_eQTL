library("SNPlocs.Hsapiens.dbSNP155.GRCh37")
library(BSgenome)
library(readr)

reformatRSIDdata <- function(dataPath){
  data <- as.data.frame(read_csv(dataPath))
  data[,18] <- gsub("\\[|\\]|'", "", data[,18])
  data <- data %>%
    # Split multiple rsids
    mutate(`0` = str_split(`0`, " ")) %>%
    unnest(`0`) %>%
    # Rename column
    dplyr::rename(rsid = `0`) %>%
    # Keep first RSID for simplicity
    distinct(variantID, .keep_all = TRUE)
  return(data)
}

CTL_eGene_snppairs <- reformatRSIDdata("output/reQTL/CTL_eGene_snppairs_rsids.csv")

FNF_eGene_snppairs <- reformatRSIDdata("output/reQTL/FNF_eGene_snppairs_rsids.csv")

