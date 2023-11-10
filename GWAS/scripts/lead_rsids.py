import pandas as pd
import sys

## Inputs
leads = sys.argv[1]
OA = sys.argv[2]
dbSNP = sys.argv[3] # Path to dbSNP file directories
prefix = sys.argv[4] # dbSNP file prefix

# Read in leads
leads = pd.read_csv(leads)
leads = leads.drop(columns = ['rsID'])


leads_all = []
for chr in range(1, 23):
    lead_chr = leads[leads['chrom'] == chr]
    if lead_chr.shape[0] == 0:
        continue
    else:

        dbSNP_chr = pd.read_csv(dbSNP + '/' + prefix + '_chr' + str(chr) + '.csv')
        print("Read in dbSNP for chr" + str(chr))
        lead_dbSNP = lead_chr.merge(dbSNP_chr, how = 'left', left_on = ['hg19pos'], right_on = ['pos'])
        
        # Replace D and I with ref and alt
        # Find entries with EA = D or I

        lead_dbSNP.loc[(lead_dbSNP["EA"] == "D") | (lead_dbSNP["EA"] == "I"), "EA"] = lead_dbSNP.loc[(lead_dbSNP["EA"] == "D") | (lead_dbSNP["EA"] == "I"), "ref"]
        lead_dbSNP.loc[(lead_dbSNP["EA"] == "D") | (lead_dbSNP["EA"] == "I"), "NEA"] = lead_dbSNP.loc[(lead_dbSNP["EA"] == "D") | (lead_dbSNP["EA"] == "I"), "alt"]

        # fill in any variants with no rsIDs 
        if len(lead_dbSNP.loc[lead_dbSNP['rsID'].isnull(), 'rsID']) != 0:
            lead_dbSNP.loc[lead_dbSNP['rsID'].isnull(), 'rsID'] = str(lead_dbSNP.loc[lead_dbSNP['rsID'].isnull(), 'chrom'].tolist()[0]) + ':' + str(lead_dbSNP.loc[lead_dbSNP['rsID'].isnull(), 'hg19pos'].tolist()[0]) + ':' + str(lead_dbSNP.loc[lead_dbSNP['rsID'].isnull(), 'EA'].tolist()[0]) + ':' + str(lead_dbSNP.loc[lead_dbSNP['rsID'].isnull(), 'NEA'].tolist()[0])
            
        lead_dbSNP = lead_dbSNP.drop(columns = ['pos', 'ref', 'alt'])

        # Remove duplicates
        lead_dbSNP_dedup = lead_dbSNP.drop_duplicates(subset = ['rsID'])
        lead_dbSNP_dedup = lead_dbSNP_dedup.drop_duplicates(subset = ['hg19pos'])

        leads_all.append(lead_dbSNP_dedup)


# Concatenate all chromosomes
if len(leads_all) > 0:
    leads_dbSNP_all_final = pd.concat(leads_all)

    rsid_column = leads_dbSNP_all_final.pop('rsID')
    leads_dbSNP_all_final.insert(0, 'rsID', rsid_column)
else:
    leads_dbSNP_all_final = pd.DataFrame(columns = ['rsID', 'chrom', 'hg19pos', 'CHR:hg19POS', 'EA', 'NEA', 
                                                    'EAF', 'FreqSE', 'MinFreq', 'MaxFreq', 'BETA', 'SE', 'metaP', 'HetISq', 'HetChiSq', 
                                                    'HetDf', 'HetPVal', 'NCASES', 'NCONTROLS', 'N', 'OR', 'OR_U95', 'OR_L95',
                                                    'EAF_summstats', 'OR_AllOA', 'Pvalue_AllOA', 'OR_KneeHipOA', 'Pvalue_KneeHipOA',
                                                    'OR_HipOA', 'Pvalue_HipOA', 'OR_KneeOA', 'Pvalue_KneeOA', 'OR_TJR', 'Pvalue_TJR',
                                                    'OR_THR', 'Pvalue_THR', 'OR_TKR', 'Pvalue_TKR', 'OR_HandOA', 'Pvalue_HandOA',
                                                    'OR_FingerOA', 'Pvalue_FingerOA', 'OR_ThumbOA', 'Pvalue_ThumbOA',
                                                    'OR_SpineOA', 'Pvalue_SpineOA', 'OR_EarlyAllOA', 'Pvalue_EarlyAllOA'])

# Write to file
leads_dbSNP_all_final.to_csv(str(OA) + '/leads/' + str(OA) + '_leads_rsid.csv', index = False)