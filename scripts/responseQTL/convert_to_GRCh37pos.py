import pandas as pd
# import numpy as np
import allel
# #from pandarallel import pandarallel

# def getPositions(variantRSIDs, chrom, dbSNP = '/proj/phanstiel_lab/References/genomes/GENCODE.GRCh37.p13/dbSNP/dbSNP155.GRCh37.p13.vcf.gz'):
    
#     def perVariantPosition(variant, vcf):
        
#         # Get index of rsid
#         rsid = variant[17]
#         rsidIndex = np.where(vcf['variants/ID'] == rsid)[0][0]
        
#         # Grab GRCh37 position
#         pos = vcf['variants/POS'][rsidIndex]
        
#         return(pos)
        
#     # Make sure input is DataFrame
#     variantRSIDs = pd.DataFrame(variantRSIDs)
    
#     # Grab chrom subset of vcf
#     vcf = allel.read_vcf(dbSNP, region = chrom, fields = ['POS', 'ID'])
    
#     # Initialize pandarallel
#     #pandarallel.initialize(progress_bar = True)

#     # Apply perVariantRSID to every row, using pandarallel apply
#     #variantPositions = variantRSIDs.parallel_apply(perVariantPosition, axis = 1, vcf = vcf)
#     variantPositions = variantRSIDs.apply(perVariantPosition, axis = 1, vcf = vcf)
#     return(variantPositions)
  
  
# variantRSIDs_CTL = pd.read_csv("/pine/scr/n/e/nekramer/CQTL_GIT/CQTL/eQTL/output/reQTL/CTL_eGene_snppairs_rsids.csv")
# variantRSIDs_FNF = pd.read_csv("/pine/scr/n/e/nekramer/CQTL_GIT/CQTL/eQTL/output/reQTL/FNF_eGene_snppairs_rsids.csv")


# variantsAll_CTL = []
# variantsAll_FNF = []

# for chr in range(1, 23):
#   print(chr)
#   variantRSIDs_CTL_subset = variantRSIDs_CTL[variantRSIDs_CTL["variant_chr"] == ("chr" + str(chr))]
#   variantPositions_CTL_subset = getPositions(variantRSIDs_CTL_subset, ("chr" + str(chr)))
#   variantsAll_CTL_subset = pd.concat([variantRSIDs_CTL_subset, variantPositions_CTL_subset], axis = 1)
#   variantsAll_CTL.append(variantsAll_CTL_subset)
  
  
#   variantRSIDs_FNF_subset = variantRSIDs_FNF[variantRSIDs_FNF["variant_chr"] == ("chr" + str(chr))]
#   variantPositions_FNF_subset = getPositions(variantRSIDs_FNF_subset, ("chr" + str(chr)))
#   variantsAll_FNF_subset = pd.concat([variantRSIDs_FNF_subset, variantPositions_FNF_subset], axis = 1)
#   variantsAll_FNF.append(variantsAll_FNF_subset)


# variantsAll_CTL = pd.DataFrame(variantsAll_CTL)
# variantsAll_FNF = pd.DataFrame(variantsAll_FNF)

# variantsAll_CTL.to_csv('output/reQTL/CTL_eGene_snppairs_rsids_GRCh37pos.csv', index = False)
# variantsAll_FNF.to_csv('output/reQTL/FNF_eGene_snppairs_rsids_GRCh37pos.csv', index = False)
vcf = allel.read_vcf('output/reQTL/FNF_sigreGenes_GRCh37dbSNP_subset.vcf.gz', fields = ['POS', 'ID'])
vcfdf = pd.DataFrame(vcf)
vcfdf.columns = ['rsid', 'GRCh37pos']
vcfdf.to_csv('output/reQTL/FNF_sigreGenes_GRCh37dbSNP_subset.csv', index = False)
