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

rule_all_inputs = [#'output/qc/multiqc_report.html',
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
nomFilter_prefix = 'output/qtl/' + filePrefix + '_nom1Mb'
nomFilter = 'output/qtl/' + filePrefix + '_nom1Mb_chr{chr}.csv'
peerMultipleTestingFinal = 'output/qtl/' + filePrefix + '_perm1Mb_FDR.csv'
peerMultipleTestingSig = 'output/qtl/' + filePrefix + '_perm1Mb_sig.csv'


if config['iteratePEER'] == 'TRUE':
    nomFilter_final = [expand(nomFilter, condition = '{condition}', Nk = '{Nk}', chr = c) for c in range(1, 23)]
else:
    nomFilter_final = [expand(nomFilter, condition = '{condition}', chr = c) for c in range(1, 23)]

# Log files
PEER_eQTL_out = 'output/logs/' + filePrefix + '_eQTL.out'
PEER_eQTL_err = 'output/logs/' + filePrefix + '_eQTL.err'
PEER_nominal_eQTL_out = 'output/logs/' + filePrefix + '_nominal_eQTL.out'
PEER_nominal_eQTL_err = 'output/logs/' + filePrefix + '_nominal_eQTL.err'
get_nomThreshold_out = 'output/logs/' + filePrefix + '_nomThreshold.out'
get_nomThreshold_err = 'output/logs/' + filePrefix + '_nomThreshold.err'
nomFilter_out = 'output/logs/' + filePrefix + '_nomFilter.out'
nomFilter_err = 'output/logs/' + filePrefix + '_nomFilter.err'
PEER_multipleTesting_perm_out_final = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm.out'
PEER_multipleTesting_perm_err_final = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm.err'
PEER_multipleTestingSig_out = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm_sig.out'
PEER_multipleTestingSig_err = 'output/logs/' + filePrefix + '_eQTL_multipleTesting_perm_sig.err'

if config['iteratePEER'] == 'TRUE':
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(nomThreshold, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerMultipleTestingSig, condition = ['CTL', 'FNF'], Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(nomFilter, condition = ['CTL', 'FNF'], Nk = n, chr = c) for n in range(1, Nk + 1, config['iterateBy']) for c in range(1, 23)]])
else:
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(nomThreshold, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(nomFilter, condition = ['CTL', 'FNF'], chr = c) for c in range(1, 23)]])
    rule_all_inputs.extend([[expand(peerMultipleTestingSig, condition = ['CTL', 'FNF'])]])
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

# Filtering nominal variants based on nominal threshold and split to separate files by chromosome
rule nominal_Filter:
    input:
        nomData = rules.PEER_nominal_eQTL.output,
        nomThreshold = rules.get_nomThreshold.output
    output:
        nomFilter_final
    params:
        version = config['pythonVersion'],
        prefix = nomFilter_prefix
    log:
        out = nomFilter_out,
        err = nomFilter_err
    shell:
        """
        module load python/{params.version}
        python3 scripts/correct_nomQTLs.py {input.nomData} {input.nomThreshold} {params.prefix} 1> {log.out} 2> {log.err}
        """

# Add multiple testing correction to permutation pass lead variants/eGenes 
rule PEER_multipleTesting_perm:
    input:
        qtlResult = peerQTL_perm
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
        Rscript scripts/correctQTLs.R {input.qtlResult} {output} 1> {log.out} 2> {log.err}
        """

# Get significant eGenes
rule sep_eGenes:
    input:
        rules.PEER_multipleTesting_perm.output   
    output:
        peerMultipleTestingSig
    params:
        version = config['Rversion'],
        threshold = config['FDRthreshold']
    log:
        out = PEER_multipleTestingSig_out,
        err = PEER_multipleTestingSig_err
    shell:
        """
        module load r/{params.version}
        Rscript scripts/separate_eGenes.R {input} {params.threshold} {output} 1> {log.out} 2> {log.err}
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