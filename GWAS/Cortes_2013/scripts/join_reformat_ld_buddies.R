library(tidyverse)
library(data.table)


as_leads <- read_csv("IGAS_2013_leads_hg38.csv")
as_summary <- fread("23749187-GCST005529-EFO_0003898.h.tsv.gz", 
                    data.table = FALSE) |> 
  dplyr::select(hm_variant_id, hm_rsid, 
                hm_chrom, hm_pos, hm_other_allele, hm_effect_allele,
                beta, standard_error, p_value) |> 
  dplyr::rename(ldbuddy_variantID_v1 = hm_variant_id,
                ldbuddy_rsID = hm_rsid,
                ldbuddy_other_allele = hm_other_allele,
                ldbuddy_effect_allele = hm_effect_allele,
                ldbuddy_beta = beta,
                ldbuddy_beta_se = standard_error,
                ldbuddy_p = p_value) |> 
  mutate(ldbuddy_variantID_v1 = gsub("_", ":", ldbuddy_variantID_v1),
         ldbuddy_variantID_v2 = paste0(hm_chrom, ":",
                                       hm_pos, ":",
                                       ldbuddy_effect_allele, ":",
                                       ldbuddy_other_allele)) |> 
  dplyr::select(-hm_chrom, -hm_pos)


ld_files <- list.files(pattern = "_ld.ld")

ld_file_data <- list()
for (file in ld_files){
  
  ld_data <- fread(file, data.table = FALSE,
                   select = c(3, 6, 7),
                   col.names = c("variantID", "ldbuddy_variantID", "ldbuddy_R2"))
  ld_file_data[[file]] <- ld_data
  
}

ld_file_data <- bind_rows(ld_file_data)
  

as_leads_ld_v1 <- as_leads |> 
  left_join(ld_file_data, by = join_by(variantID_v1 == variantID)) |> 
  mutate(actual_variantID = variantID_v1) |> 
  filter(!is.na(ldbuddy_variantID))

as_leads_ld_v2 <- as_leads |> 
  left_join(ld_file_data, by = join_by(variantID_v2 == variantID)) |> 
  mutate(actual_variantID = variantID_v2) |> 
  filter(!is.na(ldbuddy_variantID))


as_leads_ld <- bind_rows(as_leads_ld_v1, 
                         as_leads_ld_v2)

# Join back in leads that weren't found in either variantID version
missing_as_leads <- as_leads |> 
  filter(!SNP %in% as_leads_ld$SNP) |> 
  # Add in self as LD buddy
  mutate(ldbuddy_variantID = variantID_v1,
         ldbuddy_R2 = 1,
         actual_variantID = variantID_v1)

as_leads_ld_all <- bind_rows(as_leads_ld, 
                             missing_as_leads) |>
  dplyr::select(-variantID_v1, -variantID_v2) |> 
  dplyr::rename(variantID = actual_variantID) |> 
  relocate(variantID, .after = "SNP") 


as_summary_subset_v1 <- as_summary |> 
  filter(ldbuddy_variantID_v1 %in% as_leads_ld_all$ldbuddy_variantID) |> 
  dplyr::select(-ldbuddy_variantID_v2) |> 
  dplyr::rename(ldbuddy_variantID = ldbuddy_variantID_v1)

as_summary_subset_v2 <- as_summary |> 
  filter(ldbuddy_variantID_v2 %in% as_leads_ld_all$ldbuddy_variantID) |> 
  dplyr::select(-ldbuddy_variantID_v1) |> 
  dplyr::rename(ldbuddy_variantID = ldbuddy_variantID_v2)


as_summary_total <- bind_rows(as_summary_subset_v1,
                              as_summary_subset_v2)

as_leads_ld_all_final <- as_leads_ld_all |> 
  left_join(as_summary_total, by = "ldbuddy_variantID")


write_csv(as_leads_ld_all_final,
          file = "IGAS_2013_leads_hg38_ld.csv")
