#!/usr/bin/env python3
import pandas as pd
import os, shutil
import re
import glob
import yaml
import numpy as np
from kneed import KneeLocator
import os.path
from os import path
import subprocess

## Load config file
configfile: "config/config_QTLtools_eQTL.yaml"

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep = ",")

## Convert samplesheet columns to strings
samples = samples.astype(str)

## Concatenate Sequencing_Directory to Read1 and Read2 for full read paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Group Seq_Reps
samples['id'] = samples[['Proj', 'Donor']].agg('_'.join, axis=1) + '_R_' + samples[['Condition', 'Time', 'Tech_Rep']].agg('_'.join, axis=1)

## Extract grouped read1 and read2s
read1 = samples.groupby(['id'])['Read1'].apply(list).to_dict()
read2 = samples.groupby(['id'])['Read2'].apply(list).to_dict()

## Get vcf file path of post-imputed, qc'd gzipped vcf file
vcf = config["vcf"]
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_ALL_qc.vcf.gz", vcf_file).span()[0]]

## Number of PEER factors
Nk = config['PEERfactors']

rule_all_inputs = ['output/qc/multiqc_report.html',
                    'config/config_reQTL_final.yaml',
                    [expand('output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz', condition = ['ALL', 'CTL', 'FNF'])],
                    [expand('output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz.tbi', condition = ['ALL', 'CTL', 'FNF'])],
                    [expand('output/covar/{condition}_PEERfactors_k{Nk}.txt', condition = ['CTL', 'FNF'], Nk = Nk)],
                    [expand('output/covar/{condition}_PEERfactors_k{Nk}_variance.txt', condition = ['CTL', 'FNF'], Nk = Nk)]]

        
include: 'genoCovariate.smk'
batches = ['RNAKitBatch', 'RNASequencingBatch', 'genoBatch', 'DNAKitBatch']

if config['iteratePEER'] == 'TRUE':
    filePrefix = '{condition}_PEER_k{Nk}_genoPC'
else:
    filePrefix = '{condition}_PEER_kneedle_genoPC'

for b in batches:
    b_include = config[b]
    if b_include == "TRUE":
        filePrefix += '_{}'.format(b)

# Final files
peerCov = 'output/covar/' + filePrefix + '.txt'
peerQTL_perm = 'output/qtl/' + filePrefix + '_perm1Mb.txt'
peerQTL_nominal = 'output/qtl/' + filePrefix + '_nom1Mb.txt'
nomThreshold = 'output/qtl/' + filePrefix + '_nom1Mb_thresholds.csv'
nomFilter = 'output/qtl/' + filePrefix + '_nom1Mb_filter.csv'
nom_rsid = 'output/qtl/' + filePrefix + '_nom1Mb_final_rsids.csv'
peerMultipleTestingFinal = 'output/qtl/' + filePrefix + '_perm1Mb_FDR.csv'
peerMultipleTesting_rsid = 'output/qtl/' + filePrefix + '_perm1Mb_FDR_rsids.csv'

# Log files
PEER_eQTL_out = 'output/logs/' + filePrefix + '_eQTL.out'
PEER_eQTL_err = 'output/logs/' + filePrefix + '_eQTL.err'
PEER_nominal_eQTL_out = 'output/logs/' + filePrefix + '_nominal_eQTL.out'
PEER_nominal_eQTL_err = 'output/logs/' + filePrefix + '_nominal_eQTL.err'
get_nomThreshold_out = 'output/logs/' + filePrefix + '_nomThreshold.out'
get_nomThreshold_err = 'output/logs/' + filePrefix + '_nomThreshold.err'
nomFilter_out = 'output/logs/' + filePrefix + '_nomFilter.out'
nomFilter_err = 'output/logs/' + filePrefix + '_nomFilter.err'
nom_rsid_out = 'output/logs/' + filePrefix + '_nom_rsid.out'
nom_rsid_err = 'output/logs/' + filePrefix + '_nom_rsid.err'
PEER_multipleTesting_perm_out_final = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm.out'
PEER_multipleTesting_perm_err_final = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm.err'
PEER_multipleTesting_perm_rsid_out = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm_rsid.out'
PEER_multipleTesting_perm_rsid_err = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm_rsid.err'


if config['iteratePEER'] == 'TRUE':
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(nomThreshold, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(nom_rsid, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerMultipleTesting_rsid, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
else:
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(nomThreshold, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(nom_rsid, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerMultipleTesting_rsid, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand('output/covar/{condition}_PEERkneedle.txt', condition = ['CTL', 'FNF'])]])

## Define rules
rule all:
    input:
        rule_all_inputs

include: "eQTL.smk"

# Permutation pass
rule PEER_eQTL:
    input:
        vcf = rules.filterVCFvariants.output.vcf,
        vcfIndex = rules.filterVCFvariants.output.index,
        bed = rules.indexQuant_Condition.output.bed,
        bedIndex = rules.indexQuant_Condition.output.index,
        cov = peerCov
    output:
        peerQTL_perm
    params:
        version = config['QTLToolsVersion']
    log:
        out = PEER_eQTL_out,
        err = PEER_eQTL_err
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """

# Nominal pass    
rule PEER_nominal_eQTL:
    input:
        vcf = rules.filterVCFvariants.output.vcf,
        vcfIndex = rules.filterVCFvariants.output.index,
        bed = rules.indexQuant_Condition.output.bed,
        bedIndex = rules.indexQuant_Condition.output.index,
        cov = peerCov
    output:
        peerQTL_nominal
    params:
        version = config['QTLToolsVersion']
    log:
        out = PEER_nominal_eQTL_out,
        err = PEER_nominal_eQTL_err
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --nominal 1 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """ 

# Determining other significant variants using threshold determined by lead variants in permutation pass
rule get_nomThreshold:
    input:
        permData = rules.PEER_eQTL.output
    output:
        nomThreshold
    params:
        version = config['Rversion'],
        FDRthreshold = config['FDRthreshold']
    log:
        out = get_nomThreshold_out,
        err = get_nomThreshold_err
    shell:
        """
        module load r/{params.version}
        Rscript scripts/get_nomThreshold.R {input.permData} {params.FDRthreshold} {output} 1> {log.out} 2> {log.err}
        """

# Filtering nominal variants based on nominal threshold
rule nominal_Filter:
    input:
        nomData = rules.PEER_nominal_eQTL.output,
        nomThreshold = rules.get_nomThreshold.output
    output:
        nomFilter
    params:
        version = config['pythonVersion']
    log:
        out = nomFilter_out,
        err = nomFilter_err
    shell:
        """
        module load python/{params.version}
        python3 scripts/correct_nomQTLs.py {input.nomData} {input.nomThreshold} {output} 1> {log.out} 2> {log.err}
        """

# Add multiple testing correction to permutation pass lead variants/eGenes 
rule PEER_multipleTesting_perm:
    input:
        qtlResult = peerQTL_perm,
        geneInfo = rules.indexQuant_Condition.output.bed
    output:
        peerMultipleTestingFinal 
    params:
        version = config['Rversion']
    log:
        out = PEER_multipleTesting_perm_out_final,
        err = PEER_multipleTesting_perm_err_final,
    shell:
        """
        module load r/{params.version}
        Rscript scripts/correctQTLs.R {input.qtlResult} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        """

# Get LD buddies for significant lead variants

# Convert perm SNP positions to rsIDs
rule convert_perm_rsids:
    input:
        rules.PEER_multipleTesting_perm.output
    output:
        peerMultipleTesting_rsid
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix']
    log:
        out = PEER_multipleTesting_perm_rsid_out,
        err = PEER_multipleTesting_perm_rsid_err
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {output} 1> {log.out} 2> {log.err}
        """

# Convert nominal SNP position to rsIDs
rule convert_nominal_rsids:
    input:
        rules.nominal_Filter.output
    output:
        nom_rsid
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix']
    log:
        out = nom_rsid_out,
        err = nom_rsid_err
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {output} 1> {log.out} 2> {log.err}
        """

# Populate config file for response QTL workflow
rule update_reQTL:
    input:
        vcf = rules.filterVCFvariants.output.vcf,
        rna = rules.indexQuant_All.output.bed,
        genoPC = rules.genoPCA.output.pcs,
        genoPCkneedle = rules.genoPCkneedle.output
    output:
        'config/config_reQTL_final.yaml'
    log:
        out = 'output/logs/update_reQTL.out',
        err = 'output/logs/update_reQTL.err'
    run:
        # Open original config file to get values for RNAKitBatch, RNASequencingBatch, genoBatch, DNAKitBatch
        with open('config/config_QTLtools_eQTL.yaml', 'r') as f:
            config = yaml.safe_load(f)

        # Open genoPCkneedle file to get number of geno PCs
        with open('output/covar/genoPCkneedle.txt', 'r') as f:
            genoPCkneedle = int(f.readlines()[0])


        with open('config/config_reQTL.yaml', 'r') as f:
            reqtl_config = yaml.safe_load(f)
        reqtl_config['eQTL_dir'] = 'output/qtl/'
        reqtl_config['vcf'] = input.vcf
        reqtl_config['rna'] = input.rna
        reqtl_config['genoPC'] = input.genoPC
        reqtl_config['genoPCkneedle'] = genoPCkneedle
        reqtl_config['RNAKitBatch'] = config['RNAKitBatch']
        reqtl_config['RNASequencingBatch'] = config['RNASequencingBatch']
        reqtl_config['genoBatch'] = config['genoBatch']
        reqtl_config['DNAKitBatch'] = config['DNAKitBatch']

        # Write to new updated config info
        with open('config/config_reQTL_final.yaml', 'w') as f:
            yaml.dump(reqtl_config, f)

# Implement conditional pass?