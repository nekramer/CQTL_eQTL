import pandas as pd
import sys

## Inputs
lds = sys.argv[1]
OA = str(sys.argv[2])
subset = str(sys.argv[3])
dbSNP = sys.argv[4] # Path to dbSNP file directories
prefix = sys.argv[5] # dbSNP file prefix

# Read in ld buddies (joined with leads)
lds = pd.read_csv(lds)

# Separate out leads without ld buddies
ld_yes = lds[lds['ldbuddy_CHR:hg19POS'].notna()]
ld_no = lds[lds['ldbuddy_CHR:hg19POS'].isna()]
ld_no = ld_no.drop(columns = ['ldbuddy_A1', 'ldbuddy_A2'])
ld_no['ldbuddy_rsID'] = ld_no['ldbuddy_CHR:hg19POS']
ld_no['ldbuddy_ref'] = ld_no['ldbuddy_CHR:hg19POS']
ld_no['ldbuddy_alt'] = ld_no['ldbuddy_CHR:hg19POS']

ld_yes[['ld_chr', 'ld_pos']] = ld_yes['ldbuddy_CHR:hg19POS'].str.split(':', expand = True)
ld_yes = ld_yes.astype({'ld_pos': 'int'})



lds_all = [ld_no]
for chr in range(1, 23):
    ld_chr = ld_yes[ld_yes['ld_chr'] == str(chr)]
    if ld_chr.shape[0] == 0:
        continue
    else:
        dbSNP_chr = pd.read_csv(dbSNP + '/' + prefix + '_chr' + str(chr) + '.csv')
        dbSNP_chr = dbSNP_chr.rename(columns={"rsID": "ldbuddy_rsID",
                                              "ref": "ldbuddy_ref",
                                              "alt": "ldbuddy_alt"})
        print("Read in dbSNP for chr" + str(chr))


        ld_dbSNP = pd.concat([ld_chr.merge(dbSNP_chr, how = 'left', 
                        left_on = ['ld_pos', 'ldbuddy_A1', 'ldbuddy_A2'], 
                        right_on = ['pos', 'ldbuddy_ref', 'ldbuddy_alt']),
        ld_chr.merge(dbSNP_chr, how = 'left',
                        left_on = ['ld_pos', 'ldbuddy_A1', 'ldbuddy_A2'],
                        right_on = ['pos', 'ldbuddy_alt', 'ldbuddy_ref'])])
        
        # Remove duplicates that are NAs

        ld_dbSNP = ld_dbSNP.sort_values(by = 'ldbuddy_rsID')
        ld_dbSNP_dedup  = ld_dbSNP.drop_duplicates(subset = ['ldbuddy_CHR:hg19POS'])

        # Only fill in rsID NA's and drop extra columns
        ld_dbSNP_dedup.loc[ld_dbSNP_dedup['ldbuddy_rsID'].isnull(), 'ldbuddy_rsID'] = ld_dbSNP_dedup.loc[ld_dbSNP_dedup['ldbuddy_rsID'].isnull(), 'ldbuddy_CHR:hg19POS']
        ld_dbSNP_dedup = ld_dbSNP_dedup.drop(columns = ['ldbuddy_A1', 'ldbuddy_A2', 'ld_chr', 'ld_pos', 'pos'])

        lds_all.append(ld_dbSNP_dedup)


# Concatenate all chromosomes
if len(lds_all) > 0:
    lds_dbSNP_all_final = pd.concat(lds_all)
else:
    lds_dbSNP_all_final = pd.DataFrame(columns = ['rsID', 'chrom', 'hg19pos', 'CHR:hg19POS', 'EA', 'NEA', 
                                                    'EAF', 'FreqSE', 'MinFreq', 'MaxFreq', 'BETA', 'SE', 'metaP', 'HetISq', 'HetChiSq', 
                                                    'HetDf', 'HetPVal', 'NCASES', 'NCONTROLS', 'N', 'OR', 'OR_U95', 'OR_L95',
                                                    'EAF_summstats', 'OR_AllOA', 'Pvalue_AllOA', 'OR_KneeHipOA', 'Pvalue_KneeHipOA',
                                                    'OR_HipOA', 'Pvalue_HipOA', 'OR_KneeOA', 'Pvalue_KneeOA', 'OR_TJR', 'Pvalue_TJR',
                                                    'OR_THR', 'Pvalue_THR', 'OR_TKR', 'Pvalue_TKR', 'OR_HandOA', 'Pvalue_HandOA',
                                                    'OR_FingerOA', 'Pvalue_FingerOA', 'OR_ThumbOA', 'Pvalue_ThumbOA',
                                                    'OR_SpineOA', 'Pvalue_SpineOA', 'OR_EarlyAllOA', 'Pvalue_EarlyAllOA', 
                                                    'ldbuddy_CHR:hg19POS', 'ldbuddy_R2', 'ldbuddy_rsID', 'ldbuddy_ref', 'ldbuddy_alt'])

# Write to file
lds_dbSNP_all_final.to_csv(OA + '/leads/' + subset + '_'+ OA + '_leads_ld_rsid.csv', index = False)