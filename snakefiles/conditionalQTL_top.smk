#!/usr/bin/env python3
# Nicole Kramer
# 1/24/2024
# This snakemake workflow will take eQTL results from a permutation pass and
# performs conditional analysis with QTLtools. It will determine the number of independent signals
# per phenotype (i.e. eGene) and produce files detailing eGene signals and their corresponding top variants
import pandas as pd
import os

## Load config file
configfile: "config/config_conditionalQTL.yaml"

## Read in samplesheet containing eQTL files and corresponding geno, pheno, cov, and nom threshold information
eqtl_samplesheet = pd.read_csv(config["samplesheet"], sep = ",")

## Set up dictionary based on samplesheet to reference files based on wildcard
eqtl_dict = {}
for index, row in eqtl_samplesheet.iterrows():
    eqtl_dict[row["name"]] = {"perm_qtl": row["perm_qtl"], "perm_varID_col": str(row["perm_varID_col"]), "pheno_data": row["pheno_data"], "geno_data": row["geno_data"], "cov_data": row["cov_data"], "nom_thresholds": row["nom_thresholds"], "nom_data_dir": row["nom_data_dir"], "nomFile_prefix": row["nomFile_prefix"]}

rule all:
    input:
        expand('output/qtl/{name}_cond1Mb_topSignals_rsID.csv', name = eqtl_dict.keys())
        
        
# This rule will remove any header from the nom_thresholds files and 
# rewrite to a tab-delimited file to be compatible with QTLtools
rule removeHeader:
    input:
        lambda wildcards: eqtl_dict[wildcards.name]["nom_thresholds"]
    output:
        'output/qtl/{name}_nominal_thresholds_QTLtools.txt'
    params:
        Rversion = config["Rversion"]
    log:
        out = 'output/cond/{name}_removeHeader.out',
        err = 'output/cond/{name}_removeHeader.err'
    shell:
        """
        module load r/{params.Rversion}
        Rscript scripts/conditionalQTL/removeHeader_1.R {input} {output} 1> {log.out} 2> {log.err}
        """

# This rule will filter the phenotype data to match the phenotypes included in the 
# nominal threshold file, write to a bed file, and gzip and tabix it. QTLtools will error if 
# the phenotypes do not match. In this eQTL analysis, NA's were filtered from the eQTL 
# permutation pass, thus the nominal results file
# has fewer phenotypes than the original phenotype file.
rule filterPheno:
    input:
        pheno = lambda wildcards: eqtl_dict[wildcards.name]["pheno_data"],
        nomFilters = rules.removeHeader.output
    output:
        bedgz = 'output/qtl/{name}_pheno_filtered.bed.gz',
        bedtbi = 'output/qtl/{name}_pheno_filtered.bed.gz.tbi'
    params:
        Rversion = config["Rversion"],
        samtoolsVersion = config["samtoolsVersion"]
    log:
        out = 'output/cond/{name}_filterPheno.out',
        err = 'output/cond/{name}_filterPheno.err'
    shell:
        """
        module load r/{params.Rversion}
        module load samtools/{params.samtoolsVersion}

        bedFile=$(basename {output.bedgz} .gz)

        Rscript scripts/conditionalQTL/filterPheno_2.R {input.pheno} {input.nomFilters} output/qtl/${{bedFile}} && bgzip output/qtl/${{bedFile}} && tabix -p bed {output.bedgz} 1> {log.out} 2> {log.err}
        """

# This rule performs a conditional pass with the given data with QTLtools
rule conditionalQTL:
    input:
        vcf = lambda wildcards: eqtl_dict[wildcards.name]["geno_data"],
        pheno = rules.filterPheno.output.bedgz,
        cov = lambda wildcards: eqtl_dict[wildcards.name]["cov_data"],
        nomFilters = rules.removeHeader.output
    output:
        'output/qtl/{name}_conditional.txt'
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/cond/{name}_conditionalQTL.out',
        err = 'output/cond/{name}_conditionalQTL.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.pheno} --cov {input.cov} --mapping {input.nomFilters} --out {output} 1> {log.out} 2> {log.err}
        """

# This rule will reformat the results of the conditional pass and rejoin with information 
# from the original eQTL perm results (i.e. subset for sig eGenes and get gene symbols). 
# It will return results in a csv file.
rule reformatConditionalResults:
    input:
        cond = rules.conditionalQTL.output,
        perm = lambda wildcards: eqtl_dict[wildcards.name]["perm_qtl"]
    output:
        'output/qtl/{name}_cond1Mb.csv'
    params:
        Rversion = config["Rversion"]
    log:
        out = 'output/cond/{name}_reformatConditionalResults.out',
        err = 'output/cond/{name}_reformatConditionalResults.err'
    shell:
        """
        module load r/{params.Rversion}
        Rscript scripts/conditionalQTL/reformatConditionalResults_4.R {input.cond} {input.perm} {output} 1> {log.out} 2> {log.err}
        """

# This rule takes the reformatted results of the conditional pass, pulls out the top variants for each independent signal,
# and makes a separate file quantifying the number of independent signals per phenotype
rule getSignalTopVariants:
    input:
        rules.reformatConditionalResults.output
    output:
        signal_top = 'output/qtl/{name}_cond1Mb_topSignals.csv',
        n_signals = 'output/qtl/{name}_cond1Mb_nSignals.csv'
    params:
        Rversion = config["Rversion"]
    log:
        out = 'output/cond/{name}_getSignalTopVariants.out',
        err = 'output/cond/{name}_getSignalTopVariants.err'
    shell:
        """
        module load r/{params.Rversion}
        Rscript scripts/conditionalQTL/getSignalTopVariants_5.R {input} {output.signal_top} {output.n_signals} 1> {log.out} 2> {log.err}
        """

# This rule will subset the top signal results for the top signal variants that are missing rsIDs
# after joining with rsIDs from leads in permutation pass
rule determineMissing_rsIDs:
    input:
        top_signals = rules.getSignalTopVariants.output.signal_top,
        perm = lambda wildcards: eqtl_dict[wildcards.name]["perm_qtl"]
    output:
        'output/qtl/{name}_cond1Mb_missing_rsIDs.csv'
    params:
        Rversion = config["Rversion"]
    log:
        out = 'output/cond/{name}_determineMissing_rsIDs.out',
        err = 'output/cond/{name}_determineMissing_rsIDs.err'
    shell:
        """
        module load r/{params.Rversion}
        Rscript scripts/conditionalQTL/determineMissing_rsIDs_6.R {input.top_signals} {input.perm} {output} 1> {log.out} 2> {log.err}
        """

# This rule will get the rsIDs for any top variants of signals that are missing rsIDs
rule getMissing_rsIDs:
    input:
        rules.determineMissing_rsIDs.output
    output:
        'output/qtl/{name}_cond1Mb_missing_rsIDs_results.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/cond/{name}_getMissing_rsIDs.out',
        err = 'output/cond/{name}_getMissing_rsIDs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'lead' {output} 1> {log.out} 2> {log.err}
        """

# This rule will join back the results of rsIDs for missing rsIDs with the top signal results
rule joinMissing_rsIDs:
    input:
        top_signals = rules.getSignalTopVariants.output.signal_top,
        missing_rsID_results = rules.getMissing_rsIDs.output,
        perm = lambda wildcards: eqtl_dict[wildcards.name]["perm_qtl"]
    output:
        'output/qtl/{name}_cond1Mb_topSignals_rsID.csv'
    params:
        Rversion = config["Rversion"]
    log:
        out = 'output/cond/{name}_joinMissing_rsIDs.out',
        err = 'output/cond/{name}_joinMissing_rsIDs.err'
    shell:
        """
        module load r/{params.Rversion}
        Rscript scripts/conditionalQTL/joinMissing_rsIDs_8.R {input.top_signals} {input.perm} {input.missing_rsID_results} {output} 1> {log.out} 2> {log.err}
        """