import pandas as pd
import sys


# sys.argv[1]: nominal eQTL results, chromosome
# sys.argv[2]: MAF file, chromosome
# sys.argv[3]: outfile prefix
# sys.argv[4]: chromsome

nomResults_chr = pd.read_csv(sys.argv[1])

maf_chr = pd.read_csv(sys.argv[2])

# Join nomResults_chr and maf_chr by variantID
nomResults_maf_chr = nomResults_chr.merge(maf_chr, how = 'left', left_on = 'variantID', right_on = 'variantID')

# Write to file
nomResults_maf_chr.to_csv(sys.argv[3] + '_chr' + str(sys.argv[4]) + '.csv', index = False)
