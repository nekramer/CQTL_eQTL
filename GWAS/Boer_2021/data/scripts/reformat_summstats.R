library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

OA <- args[1]

OA_summarystats <- fread(paste0(OA, "/summary_stats/KP.Format.GO.FILTER.GW.", OA, ".FULL.09052019.txt.gz"),
                         data.table = FALSE)

# Remove alleles from SNPID
OA_summarystats$SNPID <- unlist(lapply(str_split(OA_summarystats$SNPID, 
                                                 '_', n = 2), `[[`, 1))

OA_summarystats <- OA_summarystats |> 
  # Make alleles uppercase
  mutate(EffectAllele = toupper(EffectAllele),
         AlternateAllele = toupper(AlternateAllele)) |> 
  # Rename columns for consistency with lead data and plotgardener
  dplyr::rename(`CHR:hg19POS` = SNPID,
                EA = EffectAllele,
                NEA = AlternateAllele,
                EAF = EffectAlleleFrequency,
                BETA = EffectSize.Beta,
                OR = EffectSize.OR,
                p = Pvalue,
                chrom = Chromosome, 
                hg19pos = Position) |> 
  relocate(chrom, .after = `CHR:hg19POS`) |> 
  relocate(hg19pos, .after = chrom)


# Split by chromosome and write to file
for (chr in 1:22){
  chrom_data <- OA_summarystats |> 
    filter(chrom == chr)
  
  write_csv(chrom_data, file = paste0(OA, "/summary_stats/", OA, "_chr", chr, ".csv"))
  
}




