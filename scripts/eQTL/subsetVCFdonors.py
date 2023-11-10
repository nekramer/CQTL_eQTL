import pandas as pd
import sys
import difflib

# sys.argv[1] = donors.txt file with VCF donors
# sys.argv[2] = string of donors from RNA samplesheet

# Read in file containing sample names from vcf file and use as array
donorFile = pd.read_csv(sys.argv[1], header = None)
vcfDonors = donorFile.iloc[:,0].values

# Convert RNA donor string to array
rnaDonors = sys.argv[2].split(",")

keepDonors = []
# Check that rnaDonors have a match in vcfDonors and append to keepDonors
for rnaDonor in rnaDonors:
  if rnaDonor not in vcfDonors:
    print("Donor " + rnaDonor + " not found in VCF file")
  else:
    keepDonors.append(rnaDonor)


samples = pd.DataFrame(keepDonors)
samples.to_csv("subset.txt", sep = " ", index = False, header = False)