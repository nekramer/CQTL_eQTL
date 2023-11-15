import sys

freqFile = sys.argv[1]
filePrefix = sys.argv[2]

# Split up chromosomes
for chr in range(1, 23):
  # Append to outfile
  with open('output/coloc/' + filePrefix + '_MAFs_chr' + str(chr) + '.csv', 'a+') as output:
    # Header
    output.write('variantID,A1,A2,maf,ma\n')
    # Iterate through freqFile
    with open(freqFile, 'r') as f:
      for line in f:
        # Split columns  
        snpdata = line.rstrip().split()
        # Skip header
        if snpdata[0] == 'CHR':
           continue
    
        chrom = snpdata[0]
        variant = snpdata[1]
        a1 = snpdata[2]
        a2 = snpdata[3]
        a1_freq = float(snpdata[4])
        a2_freq = 1 - a1_freq

        # Whichever freq is smaller is the MAF
        if a1_freq < a2_freq:
            maf = a1_freq
            ma = a1
        else:
            maf = a2_freq
            ma = a2

        # Write to output      
        if chrom == str(chr):
            output.write(variant + ',' + a1 + ',' + a2 + ',' + str(maf) + ',' + ma + '\n')