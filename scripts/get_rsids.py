import pandas as pd
import sys

## Inputs
qtl_results = sys.argv[1]
dbSNP = sys.argv[2] # Path to dbSNP file directories
prefix = sys.argv[3] # dbSNP file prefix
suffix = sys.argv[4] # dbSNP file prefix
outfile = sys.argv[5]

# Read in QTL results and grab snp chromosome, position, and alleles with variantID
qtl_pos = pd.read_csv(qtl_results)
qtl_pos[['chr', 'pos', 'a1', 'a2']] = qtl_pos['variantID'].str.split(':', expand = True)
qtl_pos = qtl_pos.astype({'pos': 'int'})


qtl_dbSNP_all = []

# Iterate through chromosomes
for chr in range(1, 23):
    print('Processing chr' + str(chr))
    # Subset qtl_pos for chr
    qtl_pos_chr = qtl_pos[qtl_pos['chr'] == 'chr' + str(chr)]

    # Read in dbSNP file for chr
    dbSNP_chr = pd.read_csv(dbSNP + '/' + prefix + '_chr' + str(chr) + '_' + suffix + '.csv')
    print('Successfully read in dbSNP file for chr' + str(chr))
    # Join qtl_pos_chr and dbSNP_chr trying to match alleles (both to ref and alt)
    qtl_dbSNP = pd.concat([qtl_pos_chr.merge(dbSNP_chr, how = 'left', 
                      left_on = ['pos', 'a1', 'a2'], 
                      right_on = ['pos', 'ref', 'alt']),
    qtl_pos_chr.merge(dbSNP_chr, how = 'left',
                      left_on = ['pos', 'a1', 'a2'],
                      right_on = ['pos', 'alt', 'ref'])])
    print('Successfully joined qtl data and dbSNP data for chr' + str(chr))
    # Remove duplicates
    qtl_dbSNP_dedup = qtl_dbSNP.drop_duplicates()

    # Remove a1 and a2 columns 
    qtl_dbSNP_final = qtl_dbSNP_dedup.drop(columns = ['a1', 'a2'])


    # Rename NA's to variantID and grab alleles
    qtl_dbSNP_final.loc[qtl_dbSNP_final['rsID'].isnull(), 'rsID'] = qtl_dbSNP_final.loc[qtl_dbSNP_final['rsID'].isnull(), 'variantID']
    qtl_dbSNP_final.loc[qtl_dbSNP_final['ref'].isnull(), 'ref'] = qtl_dbSNP_final.loc[qtl_dbSNP_final['rsID'].isnull(), 'variantID'].str.split(':').str[2]
    qtl_dbSNP_final.loc[qtl_dbSNP_final['alt'].isnull(), 'alt'] = qtl_dbSNP_final.loc[qtl_dbSNP_final['rsID'].isnull(), 'variantID'].str.split(':').str[3]

    # Append to list
    qtl_dbSNP_all.append(qtl_dbSNP_final)


# Concatenate all chromosomes
qtl_dbSNP_all_final = pd.concat(qtl_dbSNP_all)

# Write to file
print('Writing to file...')
qtl_dbSNP_all_final.to_csv(outfile, index = False)