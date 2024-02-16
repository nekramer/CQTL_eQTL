#!/usr/bin/env python3
# Nicole Kramer
# 1/31/2024
# This snakemake workflow is the second part to conditional QTL analysis. conditionalQTL_top needs to be run first
# to identify independent eGene signals and their top variants. This workflow will get LD buddies for any
# secondary top variants and get conditional nominal p-values for each separate signal. It uses the same samplesheet
# from the first workflow and references the files made based on samplesheet 'name' column.
import pandas as pd
import os, subprocess

## Load config file
configfile: "config/config_conditionalQTL.yaml"

## Read in samplesheet containing eQTL files and corresponding geno, pheno, cov, and nom threshold information
eqtl_samplesheet = pd.read_csv(config["samplesheet"], sep = ",")

## Set up dictionary based on samplesheet to reference files based on wildcard
eqtl_dict = {}
for index, row in eqtl_samplesheet.iterrows():
    eqtl_dict[row["name"]] = {"LD_perm": row["LD_perm"], "pheno_data": row["pheno_data"], "geno_data": row["geno_data"], "cov_data": row["cov_data"], "nom_thresholds": row["nom_thresholds"], "nom_data_dir": row["nom_data_dir"], "nomFile_prefix": row["nomFile_prefix"]}

## Set up dictionary from eqtl_dict to get the geno vcf file prefixes for generating PLINK LD reference files
ld_prefixes = {}
## Read in top conditional results and make dictionary of list of snps to get LD for
ld_snps = {}

## Set up dictionary of eGenes with multiple signals and their corresponding top variants
multi_signal_genes = {}

for key in eqtl_dict:
    # ld prefix
    ld_prefixes[key] = os.path.splitext(os.path.splitext(os.path.basename(eqtl_dict[key]["geno_data"]))[0])[0]
    # ld snps
    cond_results = pd.read_csv('output/qtl/' + key + '_cond1Mb_topSignals_rsID.csv')
    #secondary_cond_results = cond_results.loc[cond_results['signal'] > 0]
    ld_snps[key] = cond_results['variantID'].values.tolist()

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
        [expand('output/qtl/{name}_cond1Mb_topSignals_rsID_LD_rsID_final.csv', name = eqtl_dict.keys())],
        rule_all_nominal_file_inputs,
        [expand('output/qtl/{name}_nom1Mb_1signal_{chrom}.csv', name = eqtl_dict.keys(), chrom = "chr" + str(c)) for c in range(1, 23)],
        [expand('output/vcf/{name}_reformat_variantMAFs_{chrom}.done', name = eqtl_dict.keys(), chrom = "chr" + str(c)) for c in range(1, 23)],
        [expand('output/qtl/{name}_allSignals_nom1Mb_MAFs_{chrom}.csv', name = eqtl_dict.keys(), chrom = "chr" + str(c)) for c in range(1, 23)]


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
# for secondary signal top signal variants
rule get_LDbuddies:
    input:
        rules.makeLDref.output
    output:
        ld = temp('output/cond/ld/{name}_{snp}.ld'),
        nosex = temp('output/cond/ld/{name}_{snp}.nosex')
    params:
        prefix = lambda wildcards: ld_prefixes[wildcards.name]
    log:
        'output/cond/{name}_getLDbuddies{snp}.log'
    shell:
        """
        module load plink
        plink --bfile output/vcf/{params.prefix} --ld-snp {wildcards.snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out output/cond/ld/{wildcards.name}_{wildcards.snp}
        """

# This rule will use the output of get_LD buddies to join a top signal variant's ld buddy with information
# from the original result file
rule reformat_LDbuddies:
    input:
        buddies = rules.get_LDbuddies.output.ld,
        topSignal_data = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_topSignals_rsID.csv'
    output:
        'output/cond/ld/{name}_{snp}_ld.csv'
    params:
        version = config['Rversion']
    log:
        out = 'output/cond/{name}_reformatLDbuddies{snp}.out',
        err = 'output/cond/{name}_reformatLDbuddies{snp}.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/reformat_LDbuddies.R {input.topSignal_data} {input.buddies} {wildcards.snp} {output} 1> {log.out} 2> log.err
        """

# This rule will join the reformatted outputs of separated LD buddies back into one file 
rule join_LDbuddies:
    input:
        lambda wildcards: expand('output/cond/ld/{name}_{snp}_ld.csv', name = wildcards.name, snp = ld_snps[wildcards.name])
    output:
        'output/qtl/{name}_cond1Mb_topSignals_rsID_LD_signals.csv'
    log:
        out = 'output/cond/{name}_joinLDbuddies.out',
        err = 'output/cond/{name}_joinLDbuddies.err'
    run:
        all_variants = []
        for file in input:
            data = pd.read_csv(file)
            all_variants.append(data)

        final_data = pd.concat(all_variants)
        final_data.to_csv(output[0], index = False)

rule get_LDbuddy_rsIDs:
    input:
        rules.join_LDbuddies.output
    output:
        'output/qtl/{name}_cond1Mb_topSignals_rsID_LD_rsID_final.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/cond/{name}_getLDbuddy_rsIDs.out',
        err = 'output/cond/{name}_getLDbuddy_rsIDs.out'
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'ld' {output} 1> {log.out} 2> {log.err}
        """

# This rule will make the total top variant conditional LD file 
# by joining original LD results for first signal variants and
# the output of join_LDbuddies for secondary signal variants
# rule join_signal_LDbuddies:
#     input:
#         topSignal_data = lambda wildcards: 'output/qtl/' + wildcards.name + '_cond1Mb_topSignals_rsID.csv',
#         perm_LD = lambda wildcards: eqtl_dict[wildcards.name]["LD_perm"],
#         cond_signal_LD = rules.get_LDbuddy_rsIDs.output
#     output:
#         'output/qtl/{name}_cond1Mb_topSignals_rsID_LD_rsID_final.csv'
#     params:
#         version = config['Rversion']
#     log:
#         out = 'output/cond/{name}_join_signal_LDbuddies.out',
#         err = 'output/cond/{name}_join_signal_LDbuddies.err'
#     shell:
#         """
#         module load r/{params.version}
#         Rscript scripts/conditionalQTL/join_signal_LDbuddies_5.R {input.topSignal_data} {input.perm_LD} {input.cond_signal_LD} {output} 1> {log.out} 2> {log.err}
#         """

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