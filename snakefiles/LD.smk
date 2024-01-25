#!/usr/bin/env python3
# Nicole Kramer
# 1/23/2024
# This snakemake workflow will take a file of variants from eQTL results and will get and join LD information for these variants 
# based on an LD reference of choice. The output will contain all the same information from the input file, but each variant
# will have rows for each LD buddy, with LD R2 and allele information for each variant.
import pandas as pd
import re, os

## Load config file
configfile: "config/config_LD.yaml"

## Read in samplesheet for eQTL input files
input_files = pd.read_csv(config["samplesheet"], sep = ",")

## Grab variant lists reading in files from samplesheet and grabbing the variant column and chrom columns, then put into dictionary based on file
file_dict = {}
for index, row in input_files.iterrows():
    file_dict[row["file"]] = pd.read_csv(row["file"], sep = ",", usecols = [row["variant_col"], row["chrom_col"]])
    
## Grab file column names for the variant_col
file_var_columns = {}
for index, row in input_files.iterrows():
    file_var_columns[row["file"]] = pd.read_csv(row["file"], sep = ",", usecols = [row["variant_col"]], nrows = 1).columns[0]

## Get basename of each input eQTL file for naming output files
file_basenames = []
for key in file_dict.keys():
    file_nodir = os.path.basename(key)
    file_basename = os.path.splitext(file_nodir)[0]
    file_basenames.append(file_basename)

# Convert file dictionary to variant-based dictionary that contains the chromosome and file name for each variant
variant_dict = {}
for key, value in file_dict.items():
    for index, row in value.iterrows():
        # Get the column name for the variant_col from file_var_columns
        variant_col_name = file_var_columns[key]

        # Check if the variant id needs the chr prefix removed or not
        # And just chromosome number to put in dictionary
        if config["ld_chr_prefix"] == "FALSE":
            varID = re.sub('chr', '', row[0])
            variant_dict[varID] = {'chrom': re.sub('chr', '', row[1]), 'file': key}
        else:
            variant_dict[row[0]] = {'chrom': re.sub('chr', '', row[1]), 'file': key}


rule all:
    input:
        [expand('output/ld/{snp}_getLDbuddies.done', snp = variant_dict.keys())],
        [expand('output/ld/{fileBasename}_LD_' + config['out_suffix'] + '.csv', fileBasename = file_basenames)]

# This rule will get the LD information for each variant in the input files based on the LD reference of choice.
# The variants will each be a separate wildcard and have their own output file from PLINK. Because variants
# will not always be found in the LD reference, the output file will be a dummy file indicating PLINK ran and 
# future jobs will check for the presence of the .ld file.
rule get_LDbuddies:
    input:
        lambda wildcards: variant_dict[wildcards.snp]['file']
    output:
        temp('output/ld/{snp}_getLDbuddies.done')
    params:
        ldref_prefix = config["ldref_prefix"],
        ldref_chrom = lambda wildcards: variant_dict[wildcards.snp]['chrom']
    log:
        'output/ld/get_LDbuddies_{snp}.log'
    shell:
        """
        module load plink
        # Prevent bash strict mode from exiting on error
        set +e
        plink --bfile {params.ldref_prefix}{params.ldref_chrom} --ld-snp {wildcards.snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out output/ld/{wildcards.snp}
        touch {output}
        """

# This rule will reformat the PLINK output of LD buddies for a given variant and join it back with the data from the original eQTL file.
# Each variant is still treated separately.
rule reformat_LDbuddies:
    input:
        rules.get_LDbuddies.output
    output:
        temp('output/ld/{snp}_ld_reformat.csv')
    params:
        version = config['Rversion'],
        eqtl_file = lambda wildcards: variant_dict[wildcards.snp]['file'],
        eqtl_var_col = lambda wildcards: file_var_columns[variant_dict[wildcards.snp]['file']],
        ld_chr_prefix = config["ld_chr_prefix"]
    log:
        out = 'output/ld/reformat_LDbuddies_{snp}.out',
        err = 'output/ld/reformat_LDbuddies_{snp}.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/LD/reformat_LDbuddies_2.R {wildcards.snp} {params.ld_chr_prefix} {params.eqtl_file} {params.eqtl_var_col} {output} 1> {log.out} 2> log.err
        """

# This rule will join all the reformatted LD buddy files for each variant into one file for original file.
rule join_LDbuddies:
    input:
        lambda wildcards: expand('output/ld/{snp}_ld_reformat.csv', snp = variant_dict.keys())
    output:
        'output/ld/{fileBasename}_LD_' + config['out_suffix'] + '.csv'
    log:
        out = 'output/ld/join_LDbuddies_{fileBasename}.out',
        err = 'output/ld/join_LDbuddies_{fileBasename}.err'
    run:
        file_snps = []
        for file in input:
            data = pd.read_csv(file)
            file_snps.append(data)

        final_data = pd.concat(file_snps)
        final_data.to_csv(output[0], index = False)