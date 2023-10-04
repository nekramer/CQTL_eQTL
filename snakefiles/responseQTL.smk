#!/usr/bin/env python3

## Load config file
configfile: "config/config_reQTL_final.yaml"
vcf = config['vcf']
rna = config['rna']
Nk = config['PEERfactors']

batches = ['RNAKitBatch', 'RNASequencingBatch', 'genoBatch', 'DNAKitBatch']
peerCov = 'output/covar/ALL_PEERk{Nk}_genoPC'
fileExt = ''

for b in batches:
    b_include = config[b]
    if b_include == "TRUE":
        peerCov += '_{}'.format(b)
        fileExt += '_{}'.format(b)

peerCov += ".csv"

rule all:
    input:
        'output/reQTL/FNF_sig_reQTLs.rds',
        'output/reQTL/FNF_sig_reQTLs.csv'

# Separate eGenes that are only found in FNF
rule sep_eGenes:
    input:
        CTL_eQTL = config['eQTL_dir'] + 'CTL_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_FDR_rsids.csv',
        FNF_eQTL = config['eQTL_dir'] + 'FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_FDR_rsids.csv'   
    output:
        'output/reQTL/FNFonly_sig_eGenes.csv'
    params:
        version = config['Rversion'],
        correction = config['correction'],
        threshold = config['FDRthreshold']
    log:
        out = 'output/logs/sep_eGenes.out',
        err = 'output/logs/sep_eGenes.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/separate_eGenes.R {input.CTL_eQTL} {input.FNF_eQTL} {params.correction} {params.threshold} 1> {log.out} 2> {log.err}
        """ 

# Filter the VCF file for both sets of significant lead variants
rule subsetVCF_leadvar:
    input:
        eGenes_FNF = 'output/reQTL/FNFonly_sig_eGenes.csv',
        vcf = vcf
    output:
        'output/reQTL/FNF_leadVars.vcf.gz'
    params:
        version = config['gatkVersion']
    log:
        out = 'output/logs/subsetVCF_leadvar.out',
        err = 'output/logs/subsetVCF_leadvar.err'
    shell:
        """
        module load gatk/{params.version}
        # Grab variantID from each condition-specific eGenes list, removing header
        # 7th column in eGenes file is variantID
        cut -d "," -f7 {input.eGenes_FNF} > output/reQTL/FNF_variants.list
        sed -i '1d' output/reQTL/FNF_variants.list
        # Subset vcf for these
        gatk SelectVariants -V {input.vcf} --keep-ids output/reQTL/FNF_variants.list -O {output} 1> {log.out} 2> {log.err}
        """

# Get PEER factors for ALL norm RNA counts
rule getPEER_ALL:
    input:
        rna
    output:
        factors = [expand('output/covar/ALL_PEERfactors_k{Nk}.txt', Nk = Nk)],
        var = [expand('output/covar/ALL_PEERfactors_k{Nk}_variance.txt', Nk = Nk)]
    params:
        Nk = config['PEERfactors']
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
        [expand(peerCov, Nk = Nk)]
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
        Nk = Nk
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
        eGene = rules.sep_eGenes.output
    output:
        rds = 'output/reQTL/FNF_sig_reQTLs.rds',
        csv = 'output/reQTL/FNF_sig_reQTLs.csv'
    params:
        version = config['Rversion'],
        threshold = config['FDRthreshold']
    log:
        out = "output/logs/getreQTLs.out",
        err = "output/logs/getreQTLs.err"
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/testInteraction.R {input.normExpression} {input.covariates} {input.vcf} {input.eGene} {params.threshold} {output.rds} {output.csv} 1> {log.out} 2> {log.err}
        """

