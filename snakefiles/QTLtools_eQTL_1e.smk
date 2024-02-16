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
nomSplit = 'output/qtl/' + filePrefix + '_nom1Mb_chr{chr}.txt'
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
nomFilter_out = 'output/logs/' + filePrefix + '_nomFilter_chr{chr}.out'
nomFilter_err = 'output/logs/' + filePrefix + '_nomFilter_chr{chr}.err'
nomSplit_out = 'output/logs/' + filePrefix + '_nomSplit_chr{chr}.out'
nomSplit_err = 'output/logs/' + filePrefix + '_nomSplit_chr{chr}.err'
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
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --nominal 1 --window 1000000 --std-err --out {output} 1> {log.out} 2> {log.err}
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
        Rscript scripts/eQTL/get_nomThreshold.R {input.permData} {params.FDRthreshold} {output} 1> {log.out} 2> {log.err}
        """

rule splitNom:
    input:
        rules.PEER_nominal_eQTL.output
    output:
        nomSplit
    log:
        out = nomSplit_out,
        err = nomSplit_err
    shell:
        """
        prefix=`basename {input} .txt`
        if [ {wildcards.chr} == 1 ]; then
            awk '{{if ($2 == "chr1"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr1.txt
        elif [ {wildcards.chr} == 2 ]; then
            awk '{{if ($2 == "chr2"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr2.txt
        elif [ {wildcards.chr} == 3 ]; then
            awk '{{if ($2 == "chr3"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr3.txt
        elif [ {wildcards.chr} == 4 ]; then
            awk '{{if ($2 == "chr4"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr4.txt
        elif [ {wildcards.chr} == 5 ]; then
            awk '{{if ($2 == "chr5"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr5.txt
        elif [ {wildcards.chr} == 6 ]; then
            awk '{{if ($2 == "chr6"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr6.txt
        elif [ {wildcards.chr} == 7 ]; then
            awk '{{if ($2 == "chr7"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr7.txt
        elif [ {wildcards.chr} == 8 ]; then
            awk '{{if ($2 == "chr8"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr8.txt
        elif [ {wildcards.chr} == 9 ]; then
            awk '{{if ($2 == "chr9"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr9.txt
        elif [ {wildcards.chr} == 10 ]; then
            awk '{{if ($2 == "chr10"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr10.txt
        elif [ {wildcards.chr} == 11 ]; then
            awk '{{if ($2 == "chr11"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr11.txt
        elif [ {wildcards.chr} == 12 ]; then
            awk '{{if ($2 == "chr12"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr12.txt
        elif [ {wildcards.chr} == 13 ]; then
            awk '{{if ($2 == "chr13"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr13.txt
        elif [ {wildcards.chr} == 14 ]; then
            awk '{{if ($2 == "chr14"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr14.txt
        elif [ {wildcards.chr} == 15 ]; then
            awk '{{if ($2 == "chr15"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr15.txt
        elif [ {wildcards.chr} == 16 ]; then
            awk '{{if ($2 == "chr16"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr16.txt
        elif [ {wildcards.chr} == 17 ]; then
            awk '{{if ($2 == "chr17"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr17.txt
        elif [ {wildcards.chr} == 18 ]; then
            awk '{{if ($2 == "chr18"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr18.txt
        elif [ {wildcards.chr} == 19 ]; then
            awk '{{if ($2 == "chr19"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr19.txt
        elif [ {wildcards.chr} == 20 ]; then
            awk '{{if ($2 == "chr20"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr20.txt
        elif [ {wildcards.chr} == 21 ]; then
            awk '{{if ($2 == "chr21"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr21.txt
        elif [ {wildcards.chr} == 22 ]; then
            awk '{{if ($2 == "chr22"){{print;}}}}' {input} > output/qtl/${{prefix}}_chr22.txt
        fi

        """


# Adding nominal filter to split nominal results
rule nominal_Filter:
    input:
        nomData = rules.splitNom.output,
        nomThreshold = rules.get_nomThreshold.output
    output:
        nomFilter
    params:
        version = config['pythonVersion'],
        prefix = nomFilter_prefix
    log:
        out = nomFilter_out,
        err = nomFilter_err
    shell:
        """
        module load python/{params.version}
        python3 scripts/eQTL/correct_nomQTLs.py {input.nomData} {input.nomThreshold} {params.prefix} {wildcards.chr} 1> {log.out} 2> {log.err}
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
        Rscript scripts/eQTL/correctQTLs.R {input.qtlResult} {output} 1> {log.out} 2> {log.err}
        """

# Get significant eGenes
rule sig_eGenes:
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
        Rscript scripts/eQTL/separate_eGenes.R {input} {params.threshold} {output} 1> {log.out} 2> {log.err}
        """ 
