import pandas as pd
import sys
import numpy as np

## Inputs
qtl_results = sys.argv[1]
dbSNP = sys.argv[2] # Path to dbSNP file directories
prefix = sys.argv[3] # dbSNP file prefix

# Check length of arguments to determine if there's a suffix or not
if len(sys.argv) == 7:
    suffix = '_' + sys.argv[4] # dbSNP file suffix
    variant_type = sys.argv[5]
    outfile = sys.argv[6]
else:
    suffix = ''
    variant_type = sys.argv[4]
    outfile = sys.argv[5]

# Read in QTL results and grab snp chromosome, position, and alleles with variantID
qtl_pos = pd.read_csv(qtl_results)
if variant_type == 'lead':
    qtl_pos[['chr', 'pos', 'a1', 'a2']] = qtl_pos['variantID'].str.split(':', expand = True)
elif variant_type == 'ld':
    qtl_pos[['chr', 'pos', 'a1', 'a2']] = qtl_pos['ld_variantID'].str.split(':', expand = True)
qtl_pos = qtl_pos.astype({'pos': 'int'})

print('Read in QTL results',flush=True)
qtl_dbSNP_all = []

# Iterate through chromosomes
for chr in range(1, 23):
    print('Processing chr' + str(chr),flush=True)
    # Subset qtl_pos for chr
    qtl_pos_chr = qtl_pos[qtl_pos['chr'] == 'chr' + str(chr)]

    if qtl_pos_chr.shape[0] == 0:
        print('No variants for chr' + str(chr),flush=True)
        continue

    else:
        # Read in dbSNP file for chr
        dbSNP_chr = pd.read_csv(dbSNP + '/' + prefix + '_chr' + str(chr) + suffix + '.csv')

        print('Successfully read in dbSNP file for chr' + str(chr),flush=True)
        
        if variant_type == 'ld':
            # Rename columns to not conflict with lead variant
            dbSNP_chr = dbSNP_chr.rename(columns={"rsID": "ld_rsID",
                                                "ref": "ld_ref",
                                                "alt": "ld_alt"})
        
            # Join qtl_pos_chr and dbSNP_chr trying to match alleles (both to ref and alt)
            qtl_dbSNP = pd.concat([qtl_pos_chr.merge(dbSNP_chr, how = 'left', 
                            left_on = ['pos', 'a1', 'a2'], 
                            right_on = ['pos', 'ld_ref', 'ld_alt']),
            qtl_pos_chr.merge(dbSNP_chr, how = 'left',
                            left_on = ['pos', 'a1', 'a2'],
                            right_on = ['pos', 'ld_alt', 'ld_ref'])])

        elif variant_type == 'lead':

            # Join qtl_pos_chr and dbSNP_chr trying to match alleles (both to ref and alt)
            qtl_dbSNP = pd.concat([qtl_pos_chr.merge(dbSNP_chr, how = 'left', 
                            left_on = ['pos', 'a1', 'a2'], 
                            right_on = ['pos', 'ref', 'alt']),
            qtl_pos_chr.merge(dbSNP_chr, how = 'left',
                            left_on = ['pos', 'a1', 'a2'],
                            right_on = ['pos', 'alt', 'ref'])])
    
        print('Successfully joined qtl data and dbSNP data for chr' + str(chr),flush=True)
    
        # Remove duplicates that are NAs

        if variant_type == 'lead':
            # Sort to put NA's last
            qtl_dbSNP = qtl_dbSNP.sort_values(by = 'rsID')
            # Keep first
            qtl_dbSNP_dedup = qtl_dbSNP.drop_duplicates(subset = ['gene_id', 'variantID'])

            ## Try checking NA's just based on position 
            qtl_dbSNP_notNA = qtl_dbSNP_dedup[qtl_dbSNP_dedup['rsID'].notnull()]
            qtl_dbSNP_NAs = qtl_dbSNP_dedup[qtl_dbSNP_dedup['rsID'].isnull()]

            ## Reset rsID column
            qtl_dbSNP_NAs = qtl_dbSNP_NAs.drop(columns=['rsID'])

            print(str(qtl_dbSNP_NAs.shape[0]) + ' variantIDs not found upon initial join for chr' + str(chr) + '. Fuzzy checking just based on position..',flush=True)
            
            ## Remove alleles from dbSNP_chr 
            dbSNP_chr = dbSNP_chr.drop(columns = ['ref', 'alt'])
            # Get one rsID per position (picking first if dupicates)
            dbSNP_chr = dbSNP_chr.drop_duplicates(subset = ['pos'])
            
            qtl_dbSNP_NAs = qtl_dbSNP_NAs.merge(dbSNP_chr, how = 'left',
                            left_on = ['pos'],
                            right_on = ['pos'])
            qtl_dbSNP_NAs.loc[qtl_dbSNP_NAs['ref'].isnull(), 'ref'] = qtl_dbSNP_NAs.loc[qtl_dbSNP_NAs['ref'].isnull(), 'variantID'].str.split(':').str[2]
            qtl_dbSNP_NAs.loc[qtl_dbSNP_NAs['alt'].isnull(), 'alt'] = qtl_dbSNP_NAs.loc[qtl_dbSNP_NAs['alt'].isnull(), 'variantID'].str.split(':').str[3]
        
            print('After fuzzy checking, ' + str(qtl_dbSNP_NAs[qtl_dbSNP_NAs['rsID'].isnull()].shape[0]) + ' rsIDs not found in chr' + str(chr) + ' dbSNP. Setting rsID to variantID...', flush=True)
            # Join back together
            qtl_dbSNP_dedup_v2 = pd.concat([qtl_dbSNP_NAs, qtl_dbSNP_notNA])

            # If there are any remaining NA's, set rsID to variantID
            qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['rsID'].isnull(), 'rsID'] = qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['rsID'].isnull(), 'variantID']
            qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['ref'].isnull(), 'ref'] = qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['ref'].isnull(), 'variantID'].str.split(':').str[2]
            qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['alt'].isnull(), 'alt'] = qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['alt'].isnull(), 'variantID'].str.split(':').str[3]

            qtl_dbSNP_final = qtl_dbSNP_dedup_v2.drop(columns = ['chr', 'pos', 'a1', 'a2'])

            qtl_dbSNP_final = qtl_dbSNP_final.rename(columns={'ref': 'lead_ref', 'alt': 'lead_alt'})
            

        elif variant_type == 'ld':
            qtl_dbSNP = qtl_dbSNP.sort_values(by = 'ld_rsID')
            qtl_dbSNP_dedup = qtl_dbSNP.drop_duplicates(subset = ['gene_id', 'variantID', 'ld_variantID'])

            ## Try checking NA's based on position and if allele is contained within one of the dbSNP alleles (may be an indel with various repeats)
            qtl_dbSNP_notNA = qtl_dbSNP_dedup[qtl_dbSNP_dedup['ld_rsID'].notnull()]
            qtl_dbSNP_NAs = qtl_dbSNP_dedup[qtl_dbSNP_dedup['ld_rsID'].isnull()]

            ## Reset rsID column
            qtl_dbSNP_NAs = qtl_dbSNP_NAs.drop(columns=['ld_rsID'])

            print(str(qtl_dbSNP_NAs.shape[0]) + ' variantIDs not found upon initial join for chr' + str(chr) + '. Fuzzy checking just based on position..',flush=True)

            ## Remove alleles from dbSNP_chr
            dbSNP_chr = dbSNP_chr.drop(columns = ['ref', 'alt'])
             # Get one rsID per position (picking first if dupicates)
            dbSNP_chr = dbSNP_chr.drop_duplicates(subset = ['pos'])

            qtl_dbSNP_NAs = qtl_dbSNP_NAs.merge(dbSNP_chr, how = 'left',
                            left_on = ['pos'],
                            right_on = ['pos'])

            
            print('After fuzzy checking, ' + str(qtl_dbSNP_NAs[qtl_dbSNP_NAs['rsID'].isnull()].shape[0]) + ' rsIDs not found in chr' + str(chr) + ' dbSNP. Setting rsID to variantID...', flush=True)
            # Join back together
            qtl_dbSNP_dedup_v2 = pd.concat([qtl_dbSNP_NAs, qtl_dbSNP_notNA])

            # Only fill in rsID NA's and drop alleles
            qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['ld_rsID'].isnull(), 'ld_rsID'] = qtl_dbSNP_dedup_v2.loc[qtl_dbSNP_dedup_v2['ld_rsID'].isnull(), 'ld_variantID']
            qtl_dbSNP_dedup_v2 = qtl_dbSNP_dedup_v2.drop(columns = ['ld_ref', 'ld_alt'])

            # Remove unnecesary columns 
            qtl_dbSNP_final = qtl_dbSNP_dedup_v2.drop(columns = ['chr', 'pos', 'a1', 'a2'])
        

        # Append to list
        qtl_dbSNP_all.append(qtl_dbSNP_final)


# Concatenate all chromosomes
qtl_dbSNP_all_final = pd.concat(qtl_dbSNP_all)

# Move variant columns to front
rsid_column = qtl_dbSNP_all_final.pop('rsID')
varid_column = qtl_dbSNP_all_final.pop('variantID')
qtl_dbSNP_all_final.insert(0, 'rsID', rsid_column)
qtl_dbSNP_all_final.insert(1, 'variantID', varid_column)


# Write to file
print('Writing all chromosomes to file...',flush=True)
qtl_dbSNP_all_final.to_csv(outfile, index = False)