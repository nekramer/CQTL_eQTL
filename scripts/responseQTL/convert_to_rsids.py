import pandas as pd
import allel
from pandarallel import pandarallel

def getRsids(variantPositions, dbSNP = '/proj/phanstiel_lab/References/genomes/GENCODE.GRCh38.p13/dbSNP/dbSNP155.GRCh38.p13.vcf.gz'):
    
    def perVariantRSID(variant):
        #regionString = str(variant[9]) + ':' + str(variant[10]) + '-' + str(variant[10])
        regionString = str(variant[7]) + ':' + str(variant[8]) + '-' + str(variant[8])
        rsid = allel.read_vcf(dbSNP, region = regionString, fields = ['ID'])
        
        # Catch if rsid position not found
        if rsid is None:
            #return(str(variant[9]) + ':' + str(variant[10]))
            return(str(variant[7]) + ':' + str(variant[8]))
        else:
            return(rsid['variants/ID'])
        
    # Make sure input is DataFrame
    variantPositions = pd.DataFrame(variantPositions)
    
    # Initialize pandarallel
    pandarallel.initialize(progress_bar = True)

    # Apply perVariantRSID to every row, using pandarallel apply
    variantRSIDs = variantPositions.parallel_apply(perVariantRSID, axis = 1)
    return(variantRSIDs)
  
#variantPositions_CTL = pd.read_csv("/pine/scr/n/e/nekramer/CQTL_GIT/CQTL/eQTL/output/reQTL/CTL_sig_reQTLs.csv")
# variantPositions_CTL = pd.read_csv("/pine/scr/n/e/nekramer/CQTL_GIT/CQTL/eQTL/output/reQTL/CTL_eGene_snppairs.csv")
# variantRsids_CTL = getRsids(variantPositions_CTL)
# 
# variantsAll_CTL = pd.concat([variantPositions_CTL, variantRsids_CTL], axis = 1)
# #variantsAll_CTL.to_csv('output/reQTL/CTL_reQTLs_rsids.csv', index = False)
# variantsAll_CTL.to_csv('output/reQTL/CTL_eGene_snppairs_rsids.csv', index = False)

# variantPositions_nom = pd.read_csv("/pine/scr/n/e/nekramer/CQTL_GIT/CQTL/eQTL/output/reQTL/chr16nominal_results.csv")
# variantRsids_nom = getRsids(variantPositions_nom)
# 
# variantsAll_nom = pd.concat([variantPositions_nom, variantRsids_nom], axis = 1)
# #variantsAll_FNF.to_csv('output/reQTL/FNF_reQTLs_rsids.csv', index = False)
# variantsAll_nom.to_csv('output/reQTL/chr16nominal_results_rsids.csv', index = False)

vcf = allel.read_vcf('output/reQTL/chr16_GRCh38dbSNP.vcf.gz', fields = ['POS', 'ID'])
vcfdf = pd.DataFrame(vcf)
vcfdf.columns = ['rsid', 'GRCh37pos']
vcfdf.to_csv('output/reQTL/chr16_GRCh38dbSNP.csv', index = False)
