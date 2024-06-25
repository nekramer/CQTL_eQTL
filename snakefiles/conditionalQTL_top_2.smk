#!/usr/bin/env python3
# Nicole Kramer
# This snakemake workflow will take eQTL results from a permutation pass and
# performs conditional analysis with QTLtools. It will determine the number of independent signals
# per phenotype (i.e. eGene) and produce files detailing eGene signal top variants, also getting their 
# rsIDs and LD buddies based on the study's genotype panel.
import pandas as pd
import os, subprocess

## Load config file
configfile: "config/config_conditionalQTL.yaml"

## Read in samplesheet containing eQTL files and corresponding geno, pheno, cov, and nom threshold information
eqtl_samplesheet = pd.read_csv(config["samplesheet"], sep = ",")

## Set up dictionary based on samplesheet to reference files based on wildcard
eqtl_dict = {}
for index, row in eqtl_samplesheet.iterrows():
    eqtl_dict[row["name"]] = {"perm_qtl": row["perm_qtl"], "perm_varID_col": str(row["perm_varID_col"]), "pheno_data": row["pheno_data"], "geno_data": row["geno_data"], "cov_data": row["cov_data"], "nom_thresholds": row["nom_thresholds"], "nom_data_dir": row["nom_data_dir"], "nomFile_prefix": row["nomFile_prefix"]}

## Set up dictionary from eqtl_dict to get the geno vcf file prefixes for generating PLINK LD reference files
ld_prefixes = {}
for key in eqtl_dict:
    ld_prefixes[key] = os.path.splitext(os.path.splitext(os.path.basename(eqtl_dict[key]["geno_data"]))[0])[0]

rule all:
    input:
        [expand('output/qtl/{name}_cond1Mb_topSignals_rsID_LD.csv', name = eqtl_dict.keys())],
        [expand('output/qtl/{name}_cond1Mb_nSignals.csv', name = eqtl_dict.keys())]
        
        
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

# This rule will get the rsIDs for signal top variants
rule get_rsIDs:
    input:
        rules.getSignalTopVariants.output.signal_top
    output:
         'output/qtl/{name}_cond1Mb_topSignals_rsID.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/cond/{name}_rsIDs.out',
        err = 'output/cond/{name}_rsIDs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'lead' {output} 1> {log.out} 2> {log.err}
        """

# This rule will make LD reference files in PLINK format based on the corresponding genotyping data
rule makeLDref:
    input:
        lambda wildcards: eqtl_dict[wildcards.name]["geno_data"]
    output:
        'output/vcf/{name}_ldref.done'
    log:
        out = 'output/cond/logs/{name}_makeLDref.out',
        err = 'output/cond/logs/{name}_makeLDref.err'
    params:
        prefix = lambda wildcards: ld_prefixes[wildcards.name]
    shell:
        """
        module load plink
        plink --vcf {input} --double-id --make-bed --out output/vcf/{params.prefix} 1> {log.out} 2> {log.err}
        touch {output}
        """

# This rule will use the LD reference generated in makeLDref and run PLINK to calculate LD buddies
# for top signal variants
rule get_LDbuddies:
    input:
        ldref = rules.makeLDref.output,
        topSignal_data = rules.get_rsIDs.output
    output:
        #ld = temp('output/cond/ld/{name}_{snp}.ld'),
        #nosex = temp('output/cond/ld/{name}_{snp}.nosex')
        touch('output/cond/ld/{name}.done')
    params:
        prefix = lambda wildcards: ld_prefixes[wildcards.name]
    log:
        #'output/cond/{name}_getLDbuddies{snp}.log'
        'output/cond/{name}_getLDbuddies.log'
    run:
        topSignal_variants = pd.read_csv(str(input.topSignal_data), usecols = ['variantID'])
        topSignal_variantIDs = topSignal_variants['variantID'].tolist()
        for variant in topSignal_variantIDs:
            #p2 = subprocess.run('plink --bfile output/vcf/' + params.prefix + ' --ld-snp ' + variant + ' --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out output/cond/ld/' + wildcards.name + '_' + variant, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            p = subprocess.Popen('module load plink; plink --bfile output/vcf/' + params.prefix + ' --ld-snp ' + variant + ' --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out output/cond/ld/' + wildcards.name + '_' + variant, shell=True, stdout=subprocess.PIPE)
            p.wait()

# This rule will use the output of get_LD buddies to join a top signal variant's ld buddy with information
# from the original result file
rule reformat_LDbuddies:
    input:
        buddies = rules.get_LDbuddies.output,
        topSignal_data = rules.get_rsIDs.output
    output:
        #'output/cond/ld/{name}_{snp}_ld.csv'
        touch('output/cond/ld/{name}_reformatLDbuddies.done')
    params:
        version = config['Rversion']
    log:
        out = 'output/cond/{name}_reformatLDbuddies.out',
        err = 'output/cond/{name}_reformatLDbuddies.err'
    run:
        topSignal_variants = pd.read_csv(str(input.topSignal_data), usecols = ['variantID'])
        topSignal_variantIDs = topSignal_variants['variantID'].tolist()
        for variant in topSignal_variantIDs:
            p = subprocess.Popen('module load r/' + params.version + '; Rscript scripts/reformat_LDbuddies.R ' + str(input.topSignal_data) + ' output/cond/ld/' + wildcards.name + '_' + variant + '.ld ' + variant + ' output/cond/ld/' + wildcards.name + '_' + variant + '_ld.csv', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            p.wait()

# This rule will join the reformatted outputs of separated LD buddies back into one file 
rule join_LDbuddies:
    input:
        #lambda wildcards: expand('output/cond/ld/{name}_{snp}_ld.csv', name = wildcards.name, snp = ld_snps[wildcards.name])
        reformatted_buddies = rules.reformat_LDbuddies.output,
        topSignal_data = rules.get_rsIDs.output
    output:
        'output/qtl/{name}_cond1Mb_topSignals_rsID_LD.csv'
    log:
        out = 'output/cond/{name}_joinLDbuddies.out',
        err = 'output/cond/{name}_joinLDbuddies.err'
    run:
        topSignal_variants = pd.read_csv(str(input.topSignal_data), usecols = ['variantID'])
        topSignal_variantIDs = topSignal_variants['variantID'].tolist()
        all_variants = []
        for variant in topSignal_variantIDs:
            data = pd.read_csv('output/cond/ld/' + wildcards.name + '_' + variant + '_ld.csv')
            all_variants.append(data)

        final_data = pd.concat(all_variants)
        final_data.to_csv(output[0], index = False)