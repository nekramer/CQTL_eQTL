#!/usr/bin/env python3
# Nicole Kramer
import os
import pandas as pd
# This snakemake file tests for response eQTLs between the two conditions used in this study: PBS and FNF.
# It uses independent lead eSNP-eGene pairs identified through the part 1 of the conditional
# analysis workflow (conditional_QTL_top2.smk).

## Load config file
configfile: "config/config_reQTL.yaml"
vcf = config['vcf']

rna = config['rna']
CTL_Nk = str(config['CTL_PEERfactors']).strip()
FNF_Nk = str(config['FNF_PEERfactors']).strip()
threshold_prefix = str(int(float(config['FDRthreshold'])*100)).zfill(2)

peerCov = 'output/covar/ALL_PEER_k{Nk}_genoPC.csv'

rule all:
    input:
        expand('output/reQTL/{condition}_sig' + threshold_prefix + '_reQTLs_PEER_k{Nk}_genoPC.csv', condition = "CTL", Nk = CTL_Nk),
        expand('output/reQTL/{condition}_sig' + threshold_prefix + '_reQTLs_PEER_k{Nk}_genoPC.csv', condition = "FNF", Nk = FNF_Nk)

# Make lists of signal lead variant IDs for each condition
rule make_leadList:
    input:
        ctl_leads = lambda wildcards: expand(config['eQTL_dir'] + '{condition}_PEER_k{Nk}_genoPC_cond1Mb_topSignals_rsID.csv', condition = 'CTL', Nk = eval('CTL_Nk')),
        fnf_leads = lambda wildcards: expand(config['eQTL_dir'] + '{condition}_PEER_k{Nk}_genoPC_cond1Mb_topSignals_rsID.csv', condition = 'FNF', Nk = eval('FNF_Nk'))
    output:
        temp('output/reQTL/FNF_variants.list'),
        temp('output/reQTL/CTL_variants.list')
    log:
        out = 'output/logs/make_leadList.out',
        err = 'output/logs/make_leadList.err'
    shell:
        """
        cut -d "," -f2 {input.ctl_leads} > output/reQTL/CTL_variants.list
        sed -i '1d' output/reQTL/FNF_variants.list

        cut -d "," -f2 {input.fnf_leads} > output/reQTL/FNF_variants.list
        sed -i '1d' output/reQTL/FNF_variants.list
        """

# Filter the VCF file for significant lead variants (based on variantID for matching) 
rule subsetVCF_leadvar:
    input:
        leadvar_CTL = 'output/reQTL/CTL_variants.list',
        leadvar_FNF = 'output/reQTL/FNF_variants.list',
        vcf = vcf
    output:
        'output/qtl/CTLk' + str(CTL_Nk) + '_FNFk' + str(FNF_Nk) + '_genoPC' + fileExt + '_leadVars.vcf.gz'
    params:
        version = config['gatkVersion']
    log:
        out = 'output/logs/subsetVCF_leadvar.out',
        err = 'output/logs/subsetVCF_leadvar.err'
    shell:
        """
        module load gatk/{params.version}
        cat {input.leadvar_CTL} {input.leadvar_FNF} > output/reQTL/ALL_variants.list
        # Subset vcf for concatenated variants
        gatk SelectVariants -V {input.vcf} --keep-ids output/reQTL/ALL_variants.list -O {output} 1> {log.out} 2> {log.err}
        """

# Get PEER factors for ALL norm RNA counts
rule getPEER_ALL:
    input:
        rna
    output:
        factors = [expand('output/covar/ALL_PEERfactors_k{Nk}.txt', Nk = FNF_Nk)],
        var = [expand('output/covar/ALL_PEERfactors_k{Nk}_variance.txt', Nk = FNF_Nk)]
    params:
        Nk = config['FNF_PEERfactors']
    log:
        out = 'output/logs/ALL_getPEER.out',
        err = 'output/logs/ALL_getPEER.err'
    shell:
        """
        module load r/4.2.2
        Rscript scripts/PEERfactors.R {input} ALL {params.Nk} FALSE 1> {log.out} 2> {log.err}
        """

# Organize covariates for ALL RNA counts
rule makePEERcovar_geno_ALL:
    input:
        peer = rules.getPEER_ALL.output.factors
    output:
        [expand(peerCov, Nk = FNF_Nk)]
    params:
        version = config['Rversion'],
        genoPC = config['genoPC'],
        numgenoPCs = config['genoPCkneedle'],
        donorSamplesheet = config['donorSamplesheet'],
        dnaSamplesheet = config['dnaSamplesheet'],
        samplesheet = config['samplesheet'],
        RNAKitBatch = config['RNAKitBatch'],
        RNASequencingBatch = config['RNASequencingBatch'],
        genoBatch = config['genoBatch'],
        DNAKitBatch = config['DNAKitBatch'],
        Nk = FNF_Nk
    log:
        out = 'output/logs/makePEERcovar_geno_ALL.out',
        err = 'output/logs/makePEERcovar_geno_ALL.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/formatPEERcovariates_geno_ALL.R {input.peer} {params.donorSamplesheet} {params.dnaSamplesheet} {params.genoBatch} {params.DNAKitBatch} {params.samplesheet} {params.RNAKitBatch} {params.RNASequencingBatch} {params.genoPC} {params.numgenoPCs} {output} 1> {log.out} 2> {log.err}
        """

rule get_reQTLs:
    input:
        normExpression = rna,
        covariates = rules.makePEERcovar_geno_ALL.output,
        vcf = rules.subsetVCF_leadvar.output,
        eGene = lambda wildcards: expand('output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig.csv', condition = wildcards.condition, Nk = eval(str(wildcards.condition) + '_Nk'))
    output:
        csv = 'output/reQTL/{condition}_sig' + threshold_prefix + '_reQTLs_PEER_k{Nk}_genoPC' + fileExt + '.csv'
    params:
        version = config['Rversion'],
        threshold = config['FDRthreshold']
    log:
        out = "output/logs/{condition}k{Nk}_getreQTLs.out",
        err = "output/logs/{condition}k{Nk}_getreQTLs.err"
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/testInteraction.R {input.normExpression} {input.covariates} {input.vcf} {input.eGene} {params.threshold} {output.csv} 1> {log.out} 2> {log.err}
        """
