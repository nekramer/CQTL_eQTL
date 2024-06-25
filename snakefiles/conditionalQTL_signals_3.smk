#!/usr/bin/env python3
# Nicole Kramer
# This snakemake workflow is the second part to conditional QTL analysis. conditionalQTL_top needs to be run first
# to identify independent eGene signals and their top variants. This workflow will get conditional nominal p-values for each separate signal and
# join them with their minor allele frequencies based on the input genotyping panel.
# It uses the same samplesheet from the first workflow and references the files made based on samplesheet 'name' column.

import pandas as pd
import os, subprocess

## Load config file
configfile: "config/config_conditionalQTL.yaml"

## Read in samplesheet containing eQTL files and corresponding geno, pheno, cov, and nom threshold information
eqtl_samplesheet = pd.read_csv(config["samplesheet"], sep = ",")

## Set up dictionary based on samplesheet to reference files based on wildcard
eqtl_dict = {}
for index, row in eqtl_samplesheet.iterrows():
    eqtl_dict[row["name"]] = {"pheno_data": row["pheno_data"], "geno_data": row["geno_data"], "cov_data": row["cov_data"], "nom_thresholds": row["nom_thresholds"], "nom_data_dir": row["nom_data_dir"], "nomFile_prefix": row["nomFile_prefix"]}

## Set up dictionary of eGenes with multiple signals and their corresponding top variants
multi_signal_genes = {}

for key in eqtl_dict:
    cond_results = pd.read_csv('output/qtl/' + key + '_cond1Mb_topSignals_rsID.csv')

    # multi signal genes
    cond_results['n_signal'] = cond_results.groupby('gene_id')['gene_id'].transform('count')
    multiple_signal_genes = cond_results.loc[cond_results['n_signal'] > 1]
    multiple_signal_genes_dict = multiple_signal_genes.groupby(['gene_id'])['variantID'].apply(lambda grp: list(grp.value_counts().index)).to_dict()
    multi_signal_genes[key] = multiple_signal_genes_dict


rule_all_nominal_file_inputs = []
rule_join_all_signals_nominal_inputs = {}
for key in multi_signal_genes:
    rule_join_all_signals_nominal_inputs[key] = []
    for gene in multi_signal_genes[key]:
        rule_all_nominal_file_inputs.extend([expand('output/cond/{name}_{gene}_{var}_conditionalQTL_nominal.txt', name = key, gene = gene, var = multi_signal_genes[key][gene])])
        rule_join_all_signals_nominal_inputs[key].extend(expand('output/cond/{name}_{gene}_{var}_conditionalQTL_nominal_reformat.csv', name = key, gene = gene, var = multi_signal_genes[key][gene]))

rule all:
    input:
        rule_all_nominal_file_inputs,
        [expand('output/qtl/{name}_nom1Mb_1signal_{chrom}.csv', name = eqtl_dict.keys(), chrom = "chr" + str(c)) for c in range(1, 23)],
        [expand('output/vcf/{name}_reformat_variantMAFs_{chrom}.done', name = eqtl_dict.keys(), chrom = "chr" + str(c)) for c in range(1, 23)],
        [expand('output/qtl/{name}_allSignals_nom1Mb_MAFs_{chrom}.csv', name = eqtl_dict.keys(), chrom = "chr" + str(c)) for c in range(1, 23)]


# This rule will identify eGenes with more than 1 signal and create dummy files containing their gene ids
# to be compatible with the QTLtools flag --include-phenotypes 
rule multiSignal_phenos:
    input:
        n_signals = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_nSignals.csv'
    output:
        'output/cond/{name}_multiSignal_phenos.done'
    params:
        version = config['Rversion']
    log:
        out = 'output/cond/{name}_multiSignal_phenos.out',
        err = 'output/cond/{name}_multiSignal_phenos.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/conditionalQTL/multiSignal_phenos_6.R {input.n_signals} 1> {log.out} 2> {log.err}
        touch {output}
        """

# This rule will write the top variants for signals to a single column text file for 
# input with GATK SelectVariants
rule make_topVars_list:
    input:
        condSignals = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_topSignals_rsID.csv'
    output:
        'output/cond/{name}.list'
    params:
        version = config['Rversion']
    log:
        out = 'output/cond/{name}_make_topVars_list.out',
        err = 'output/cond/{name}_make_topVars_list.out'
    shell:
        """
        module load r/{params.version} 
        Rscript scripts/conditionalQTL/make_topVars_list.R {input} {output} 1> {log.out} 2> {log.err}
        """


# This rule will filter the vcf for top variants for ease of getting
# variant genotypes when adding them to covariate files
rule split_geno_topVars:
    input:
        vcf = lambda wildcards: eqtl_dict[wildcards.name]["geno_data"],
        lead_list = rules.make_topVars_list.output
    output:
        v = 'output/cond/{name}_signalVars.vcf.gz',
        i = 'output/cond/{name}_signalVars.vcf.gz.tbi'
    params:
        version = config['gatkVersion']
    log:
        out = 'output/cond/{name}_split_geno_topVars.out',
        err = 'output/cond/{name}_split_geno_topVars.out'
    shell:
        """
        module load gatk/{params.version}
        gatk SelectVariants -V {input.vcf} --keep-ids {input.lead_list} -O {output.v} 1> {log.out} 2> {log.err}

        """

# This rule will make covariate files for each signal conditioning on other top variants
rule makeCovar_files:
    input:
        condSignals = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_topSignals_rsID.csv',
        covData = lambda wildcards: eqtl_dict[wildcards.name]["cov_data"],
        vcf_subset = rules.split_geno_topVars.output.v
    output:
        'output/cond/{name}_makeCovar_files.done'
    params:
        version = config['Rversion']
    log:
        out = 'output/cond/{name}_makeCovar_files.out',
        err = 'output/cond/{name}_makeCovar_files.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/conditionalQTL/makeCovar_files_8.R {input.condSignals} {input.covData} {input.vcf_subset} {wildcards.name} 1> {log.out} 2> {log.err}
        touch {output}
        """

# This rule will re-run QTLtools nominal pass with conditioned covariate files for each
# eGene with multiple signals
rule conditionalQTL_nominal:
    input:
        pheno = lambda wildcards: eqtl_dict[wildcards.name]["pheno_data"],
        cov = rules.makeCovar_files.output,
        multiSignal_phenos = rules.multiSignal_phenos.output,
        vcf = lambda wildcards: eqtl_dict[wildcards.name]["geno_data"]
    output:
        'output/cond/{name}_{gene}_{var}_conditionalQTL_nominal.txt'
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/cond/{name}_{gene}_{var}_conditionalQTL_nominal.out',
        err = 'output/cond/{name}_{gene}_{var}_conditionalQTL_nominal.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.pheno} \
            --cov output/covar/{wildcards.name}_{wildcards.gene}_{wildcards.var}.txt \
            --nominal 1 --window 1000000 --std-err \
            --include-phenotypes output/cond/{wildcards.gene}.txt --out {output} 1> {log.out} 2> {log.err}
        """

# This rule will reformat all the signal nominal results for each eGene/variant pair by adding a column indicating which
# signal a row's variant p-value belongs to.
rule reformat_conditionalQTL_nominal:
    input:
        condSignals = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_topSignals_rsID.csv',
        nomFile = rules.conditionalQTL_nominal.output, 
        nomThresholds = lambda wildcards: eqtl_dict[wildcards.name]["nom_thresholds"]
    output:
        'output/cond/{name}_{gene}_{var}_conditionalQTL_nominal_reformat.csv'
    params:
        version = config['Rversion']
    log:
        out = 'output/cond/{name}_{gene}_{var}_reformat_conditionalQTL_nominal.out',
        err = 'output/cond/{name}_{gene}_{var}_reformat_conditionalQTL_nominal.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/conditionalQTL/reformat_conditionalQTL_nominal_10.R {wildcards.name} {wildcards.gene} {wildcards.var} {input.condSignals} {input.nomThresholds} {output} 1> {log.out} 2> {log.err}
        """

# This rule will isolate the original nominal results for eGenes with 1 signal
# and join them with the conditional nominal p-values of eGenes with more than 1
# signal. This will be done by chromosome.
rule get_1signal_nom:
    input:
        condSignals = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_topSignals_rsID.csv'
    output:
        'output/qtl/{name}_nom1Mb_1signal_{chrom}.csv'
    params:
        version = config['Rversion'],
        chrom = lambda wildcards: wildcards.chrom,
        nom_data_dir = lambda wildcards: eqtl_dict[wildcards.name]["nom_data_dir"],
        nomFile_prefix = lambda wildcards: eqtl_dict[wildcards.name]["nomFile_prefix"]
    log:
        out = 'output/cond/{name}_{chrom}_get_1signal_nom.out',
        err = 'output/cond/{name}_{chrom}_get_1signal_nom.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/conditionalQTL/get_1signal_nom_11.R {input.condSignals} {wildcards.chrom} {params.nom_data_dir} {params.nomFile_prefix} {output} 1> {log.out} 2> {log.err}
        """

# This rule will join the reformatted nominal conditional results of eGenes with more than 1 signal
# with the isolated nominal results of eGenes with one signal. This will be done by chromosome.
rule join_all_signals_nominal:
    input:
        condSignals = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_topSignals_rsID.csv',
        cond_nom_results = lambda wildcards: rule_join_all_signals_nominal_inputs[wildcards.name],
        signal1_nom_results = rules.get_1signal_nom.output
    output:
        'output/qtl/{name}_nom1Mb_allSignals_{chrom}.csv'
    params:
        version = config['Rversion'],
        chrom = lambda wildcards: wildcards.chrom
    log:
        out = 'output/cond/{name}_{chrom}_join_all_signals_nominal.out',
        err = 'output/cond/{name}_{chrom}_join_all_signals_nominal.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/conditionalQTL/join_all_signals_nominal_12.R {input.condSignals} {wildcards.name} {wildcards.chrom} {input.signal1_nom_results} {output} 1> {log.out} 2> {log.err}
        """

# This rule will use the input genotyping data to get variant MAFs with PLINK --freq
rule get_variantMAFs:
    input:
        vcf = lambda wildcards: eqtl_dict[wildcards.name]["geno_data"]
    output:
         'output/vcf/{name}_get_variantMAFs.done'
    params:
        prefix = lambda wildcards: ld_prefixes[wildcards.name]
    log:
        out = 'output/cond/{name}_get_variantMAFs.out',
        err = 'output/cond/{name}_get_variantMAFs.err'
    shell:
        """
        module load plink
        if ! test -f output/vcf/{params.prefix}.frq; then
            plink --vcf {input.vcf} --freq --out output/vcf/{params.prefix} 1> {log.out} 2> {log.err}
        fi
        touch {output}
        """

# This rule will reformat the output of PLINK .frq MAFs to csv files that are split by chromosome
rule reformat_variantMAFs:
    input:
        rules.get_variantMAFs.output
    output:
        [expand('output/vcf/{{name}}_MAFs_{chrom}.csv', chrom = 'chr' + str(c)) for c in range(1, 23)] 
    params:
        version = config['pythonVersion'],
        prefix = lambda wildcards: ld_prefixes[wildcards.name]
    log:
        out = 'output/cond/{name}_reformat_variantMAFs.out',
        err = 'output/cond/{name}_reformat_variantMAFs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/reformat_variantMAFs.py output/vcf/{params.prefix}.frq {wildcards.name} 1> {log.out} 2> {log.err}
        """

# This rule will join the minor alleles and minor allele frequency of variants to each of the
# final nominal result files, split by chromosome
rule join_all_signals_nominal_MAFs:
    input:
        all_signals_nominal_chrom = lambda wildcards: 'output/qtl/' + wildcards.name + '_nom1Mb_allSignals_{chrom}.csv',
        mafs_chrom = lambda wildcards: 'output/vcf/' + wildcards.name + '_MAFs_' + wildcards.chrom + '.csv' 
    output:
        'output/qtl/{name}_allSignals_nom1Mb_MAFs_{chrom}.csv'
    params:
        version = config['pythonVersion'],
        filePrefix = lambda wildcards: 'output/qtl/' + wildcards.name + '_allSignals_nom1Mb'
    log:
        out = 'output/cond/{name}_{chrom}_join_all_signals_nominal_MAFs.out',
        err = 'output/cond/{name}_{chrom}_join_all_signals_nominal_MAFs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/join_nomMAFs.py {input.all_signals_nominal_chrom} {input.mafs_chrom} {params.filePrefix} {wildcards.chrom} 1> {log.out} 2> {log.err}
        """