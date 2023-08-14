import allel
import pandas as pd
from pandarallel import pandarallel
import sys
import numpy as np

# Function to convert variant positions to variant rsids based on a dbSNP reference vcf
def getRsids(variantPositions, dbSNP):
    
    def perVariantRSID(variant):
      a1 = str(variant[7]).split(":")[2]
      a2 = str(variant[7]).split(":")[3]
        
      regionString = str(variant[8]) + ':' + str(variant[9]) + '-' + str(variant[9])
      rsid = allel.read_vcf(dbSNP, region = regionString, fields = ['ID', 'REF', 'ALT'])
      
      # Catch if rsid position not found
      if rsid is None:
        return(str(variant[8]) + ':' + str(variant[9]) + ":" + a1 + ":" + a2)
      else:

        # Check for matching alleles
        
        refs = []
        alts = []
        
        # Expand combinations of alleles when there are multiple
        for i in range(len(rsid['variants/REF'])):
          
          if isinstance(rsid['variants/REF'][i], np.ndarray):
            refalleles = list(filter(None, rsid['variants/REF'][i]))
          else:
            refalleles = rsid['variants/REF'][i]
            
          if isinstance(rsid['variants/ALT'][i], np.ndarray):
            altalleles = list(filter(None, rsid['variants/ALT'][i]))
          else:
            altalleles = rsid['variants/ALT'][i]  
          
          expandedalleles = np.meshgrid(refalleles, altalleles)
          refs.append(expandedalleles[0].flatten().tolist())
          alts.append(expandedalleles[1].flatten().tolist())

        ID = None  
        # Iterate through expanded groups of alleles
        for i in range(len(refs)):
          ref_group = refs[i]
          alt_group = alts[i]
          
          if (a1 in ref_group and a2 in alt_group) or (a2 in ref_group and a1 in alt_group):  
            ID = rsid['variants/ID'][i]
        
        if ID == None:
          return(str(variant[8]) + ':' + str(variant[9]) + ":" + a1 + ":" + a2)
        else:
          return(ID)
          
    # Make sure input is DataFrame
    variantPositions = pd.DataFrame(variantPositions)
    
    # Initialize pandarallel
    pandarallel.initialize(progress_bar = False)

    # Apply perVariantRSID to every row, using pandarallel apply
    variantRSIDs = variantPositions.parallel_apply(perVariantRSID, axis = 1)
    #variantRSIDs = variantPositions.apply(perVariantRSID, axis = 1)
    return(variantRSIDs)
  

chrom = sys.argv[2]

nom_maf_results = pd.read_csv(sys.argv[1])

nom_maf_rsids = getRsids(nom_maf_results, dbSNP = '/work/users/n/e/nekramer/References/GRCh38.p13/dbSNP/dbSNP155.GRCh38.p13_chr' + str(chrom) + '.vcf.gz')

nom_maf_results_rsIDs = pd.concat([nom_maf_results, nom_maf_rsids], axis = 1)
nom_maf_results_rsIDs.to_csv(sys.argv[3], index = False)

  
