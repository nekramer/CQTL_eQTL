import pandas as pd
import sys

# sys.argv[1] = donors.txt file with VCF donors

# Read in file containing sample names from vcf file
donorFile = pd.read_csv(sys.argv[1], header = None)

# String split donorFile sample names and grab 4th element, which just contains donor name
donorFile["new_donor"] = donorFile[0].str.split("_", expand = True)[3]

# Write to file for use in bcftools
donorFile.to_csv("samples.txt", sep = " ", index = False, header = False)