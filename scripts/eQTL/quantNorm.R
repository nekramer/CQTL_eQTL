#!/usr/bin/R
library(tximeta)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(DESeq2)
library(edgeR)
library(plyranges)

# FUNCTIONS -------------------------------------------------------------------
# Function to get most upstream transcript TSS
gene_tss <- function(geneRow, txdb_transcripts){
  
  gene_transcripts <- geneRow[[7]]
  strand <- geneRow[[4]]

  gene_transcript_info <- txdb_transcripts |> 
    filter(tx_name_trunc %in% gene_transcripts) 
  
  # Get most upstream transcript
  if (strand == "-"){
    upstream_transcript <- gene_transcript_info |> 
      filter(start == max(start)) |> 
      # Just pick first one if multiple
      distinct(tx_name, .keep_all = TRUE)
  } else if (strand == "+"){
    upstream_transcript <- gene_transcript_info |> 
      filter(start == min(start)) |> 
      # Just pick first one if multiple
      distinct(tx_name, .keep_all = TRUE)
  }

  tss <- unique(upstream_transcript$start)
  print(tss)
  if (length(tss) == 0) {
    tss <- NA
  }

  return(tss)
}

# Function to inverse normalize a row of gene counts
inverseNormGene <- function(geneRow){
  normValues <- qnorm((rank(as.numeric(geneRow),
                            na.last = "keep") - 0.5)/sum(!is.na(as.numeric(geneRow))))
  return(normValues)
}

# READ IN ---------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
condition <- args[2]
sampleQuants <- args[3:length(args)]
samples <- basename(dirname(sampleQuants))

coldata <- read_csv(samplesheet) %>%
  # Get distinct samples
  distinct(Sample, .keep_all = TRUE) %>%
  # Put in order of input quants
  arrange(match(Sample, samples)) %>%
  # add files
  mutate(files = sampleQuants) %>%
  # names column
  mutate(names = Sample) %>%
  # Condition groups
  mutate(Condition = as.factor(Condition)) %>% 
  mutate(group = as.numeric(Condition))

# Import salmon transcript quantification-------------------------------------
se <- tximeta(coldata)

# Convert to gene-level scaled transcripts -------------------------------------
gse <- summarizeToGene(se)

# Filter out lowly expressed genes ---------------------------------------------

# At least 10 counts in more than 5% of samples from all
keep <- rowSums(assay(gse) >= 10) >= ceiling(nrow(colData(gse))*0.05)

gse_filtered <- gse[keep,]

# TMM normalization ---------------------------------------------------------
gse_quant <- calcNormFactors(gse_filtered, method = "TMM")

# Grab CPM counts -----------------------------------------------------------
CQTL_CPMadjTMM <- as.data.frame(cpm(gse_quant))

# Gene info -----------------------------------------------------------------
gene_info <- as.data.frame(rowRanges(gse_filtered)) |> 
  dplyr::select(seqnames, start, end, strand, gene_id, gene_name, tx_ids, gene_biotype) |> 
  # Filter out pseudogenes
  filter(!str_detect(gene_biotype, "pseudogene")) |> 
  dplyr::select(-gene_biotype) |> 
  mutate(seqnames = paste0("chr", seqnames)) |> 
  mutate(length = end - start)


# Get TSS for each gene ---------------------------------------------------

# Grab all transcripts from TxDb
txdb_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene) |> 
  as.data.frame() |> 
  mutate(tx_name_trunc = gsub("\\..*", "", tx_name))
  
# Apply fxn to each gene row
tss <- apply(gene_info, 1, gene_tss, txdb_transcripts = txdb_transcripts)

# Replace gene start with TSS
gene_info$start <- as.numeric(tss)

gene_info <- gene_info |> 
  # Filter any NA
  filter(!is.na(start)) |>
  # Update end
  mutate(end = start + 1) |>
  #mutate(end = start + length) |> 
  # Remove transcript column 
  dplyr::select(-tx_ids)

# Inverse normal transformation -----------------------------------------------

if (condition == "CTL"){
  
  # Grab CTL
  CTL_CPMadjTMM <- dplyr::select(CQTL_CPMadjTMM, contains("CTL"))  %>%
    rename_with(.fn = ~ unlist(lapply(str_split(.x, "_"), `[[`, 2)))
  
  # Inverse normalize
  CTL_CPMadjTMM_invNorm <- as.data.frame(t(apply(CTL_CPMadjTMM, 1, inverseNormGene))) %>%
    rownames_to_column("gene_id")
  colnames(CTL_CPMadjTMM_invNorm) <- c("gene_id", colnames(CTL_CPMadjTMM))
  
  # Join with gene info
  CTL_CPMadjTMM_invNorm <- CTL_CPMadjTMM_invNorm |>  
    left_join(gene_info, by = "gene_id") |> 
    dplyr::select(-gene_name) |> 
    relocate(seqnames, start, end, gene_id, length, strand) |> 
    dplyr::rename(`#chr` = seqnames) |> 
    arrange(`#chr`, start) |>
    write_delim(file = "output/normquant/CTL_CPMadjTMM_invNorm.bed", 
                delim = "\t")
  
} else if (condition == "FNF"){
  # Grab FNF
  FNF_CPMadjTMM <- dplyr::select(CQTL_CPMadjTMM, contains("FNF")) %>%
    rename_with(.fn = ~ unlist(lapply(str_split(.x, "_"), `[[`, 2)))
  
  # Inverse normalize
  FNF_CPMadjTMM_invNorm <- as.data.frame(t(apply(FNF_CPMadjTMM, 1, inverseNormGene))) %>%
    rownames_to_column("gene_id")
  colnames(FNF_CPMadjTMM_invNorm) <- c("gene_id", colnames(FNF_CPMadjTMM))
  
  # Join with gene info
  FNF_CPMadjTMM_invNorm <- FNF_CPMadjTMM_invNorm |>  
    left_join(gene_info, by = "gene_id") |> 
    dplyr::select(-gene_name) |> 
    relocate(seqnames, start, end, gene_id, length, strand) |> 
    dplyr::rename(`#chr` = seqnames) |> 
    arrange(`#chr`, start) |>
    write_delim(file = "output/normquant/FNF_CPMadjTMM_invNorm.bed", 
              delim = "\t")
  
} else if (condition == "ALL"){
  
  # Inverse normalize across all data
  CPMadjTMM_invNorm <- as.data.frame(t(apply(CQTL_CPMadjTMM, 1, inverseNormGene))) %>%
    rownames_to_column("gene_id")
  colnames(CPMadjTMM_invNorm) <- c("gene_id", colnames(CQTL_CPMadjTMM))
  
  # Join with gene info
  CPMadjTMM_invNorm <- CPMadjTMM_invNorm |> 
    left_join(gene_info, by = "gene_id") |> 
    dplyr::select(-gene_name) |> 
    relocate(seqnames, start, end, gene_id, length, strand) |> 
    dplyr::rename(`#chr` = seqnames) |> 
    arrange(`#chr`, start) |>
    write_delim(file = "output/normquant/ALL_CPMadjTMM_invNorm.bed", 
              delim = "\t")
}