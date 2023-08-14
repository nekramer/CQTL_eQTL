# utils.py>
import allel
import pandas as pd
from pandarallel import pandarallel

# Function to convert variant positions to variant rsids based on a dbSNP reference vcf
def getRsids(variantPositions, dbSNP = '/proj/phanstiel_lab/References/genomes/GENCODE.GRCh38.p13/dbSNP/dbSNP155.GRCh38.p13.vcf.gz'):
    
    def perVariantRSID(variant):
        regionString = str(variant[9]) + ':' + str(variant[10]) + '-' + str(variant[10])
        rsid = allel.read_vcf(dbSNP, region = regionString, fields = ['ID'])
        
        # Catch if rsid position not found
        if rsid is None:
            return(str(variant[9]) + ':' + str(variant[10]))
        else:
            return(rsid['variants/ID'])
        
    # Make sure input is DataFrame
    variantPositions = pd.DataFrame(variantPositions)
    
    # Initialize pandarallel
    pandarallel.initialize(progress_bar = False)

    # Apply perVariantRSID to every row, using pandarallel apply
    variantRSIDs = variantPositions.parallel_apply(perVariantRSID, axis = 1)
    return(variantRSIDs)
