#!/usr/bin/env python3
import os
import pandas as pd

## Load config file
configfile: "config/config_reQTL_final.yaml"
vcf = config['vcf']
ldref = config['ldref']
ldref_file = os.path.basename(ldref)
ldref_prefix = ldref_file[:re.search("_ALL_qc.vcf.gz", ldref_file).span()[0]]
rna = config['rna']
CTL_Nk = str(config['CTL_PEERfactors']).strip()
FNF_Nk = str(config['FNF_PEERfactors']).strip()

batches = ['RNAKitBatch', 'RNASequencingBatch', 'genoBatch', 'DNAKitBatch']
peerCov = 'output/covar/ALL_PEER_k{Nk}_genoPC'
fileExt = ''

for b in batches:
    b_include = config[b]
    if b_include == "TRUE":
        peerCov += '_{}'.format(b)
        fileExt += '_{}'.format(b)

peerCov += ".csv"

CTL_snps = pd.read_csv(config['eQTL_dir'] + 'CTL_PEER_k' + str(CTL_Nk) + '_genoPC' + fileExt + '_perm1Mb_sig.csv', usecols = ['variantID'])['variantID'].values.tolist()
FNF_snps = pd.read_csv(config['eQTL_dir'] + 'FNF_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '_perm1Mb_sig.csv', usecols = ['variantID'])['variantID'].values.tolist()

rule all:
    input:
        expand('output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD_rsID.csv', condition = "CTL", Nk = CTL_Nk),
        expand('output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD_rsID.csv', condition = "FNF", Nk = FNF_Nk),
        expand('output/ld/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_{snp}_ld.csv', condition = "CTL", Nk = CTL_Nk, snp = CTL_snps),
        expand('output/ld/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_{snp}_ld.csv', condition = "FNF", Nk = FNF_Nk, snp = FNF_snps),
        expand('output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD.csv', condition = "CTL", Nk = CTL_Nk),
        expand('output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD.csv', condition = "FNF", Nk = FNF_Nk),
        'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '.rds',
        'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '.csv',
        'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '_LD_rsID.csv'

# Get rsIDs for significant lead variants
rule sig_rsIDs:
    input:
        lambda wildcards: expand(config['eQTL_dir'] + '{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig.csv', condition = wildcards.condition, Nk = eval(str(wildcards.condition) + '_Nk'))
    output:
        'output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + 'sig_rsIDs.out',
        err = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + 'sig_rsIDs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/responseQTL/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'lead' {output} 1> {log.out} 2> {log.err}
        """

# Get rsIDs for nominal results
rule sig_rsIDs_nom:
    input:
        lambda wildcards: expand(config['eQTL_dir'] + '{condition}_PEER_k{Nk}_genoPC' + fileExt + '_nom1Mb_chr{chr}.csv', condition = wildcards.condition, Nk = eval(str(wildcards.condition) + '_Nk'), chr = wildcards.chr)
    output:
        'output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_nom1Mb_chr{chr}_rsID.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_nom1Mb_chr{chr}_rsIDs.out',
        err = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_nom1Mb_chr{chr}_rsIDs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'lead' {output} 1> {log.out} 2> {log.err}
        """

# Make LD reference file based on input vcf
rule makeLDref:
    input:
        ldref
    output:
        bed = 'output/vcf/' + ldref_prefix + '.bed',
        bim = 'output/vcf/' + ldref_prefix + '.bim',
        fam = 'output/vcf/' + ldref_prefix + '.fam'
    log:
        out = 'output/vcf/logs/makeLDref.out',
        err = 'output/vcf/logs/makeLDref.err'
    params:
        prefix = ldref_prefix
    shell:
        """
        module load plink
        plink --vcf {input} --double-id --make-bed --out output/vcf/{params.prefix} 1> {log.out} 2> {log.err}
        """

#Get LD buddies for significant lead variants
rule get_sigLDbuddies:
    input:
        leads = lambda wildcards: expand(config['eQTL_dir'] + '{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID.csv', condition = wildcards.condition, Nk = eval(str(wildcards.condition) + '_Nk')),
        ldref_bed = rules.makeLDref.output.bed,
        ldref_bim = rules.makeLDref.output.bim,
        ldref_fam = rules.makeLDref.output.fam
    output:
        'output/ld/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_{snp}_ld.ld'
    params:
        filePrefix = '{condition}_PEER_k{Nk}_genoPC' + fileExt
    log:
        'output/ld/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_get_sigLDbuddies_{snp}_ld.log'
    shell:
        """
        module load plink
        # Get file prefix of ldref for plink
        ld_path={input.ldref_bed}
        ldref_prefix=${{ld_path%%.*}}

        plink --bfile ${{ldref_prefix}} --ld-snp {wildcards.snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out output/ld/{params.filePrefix}_perm1Mb_sig_rsID_{wildcards.snp}_ld

        """

rule reformat_sigLDbuddies:
    input:
        buddies = rules.get_sigLDbuddies.output,
        leads = lambda wildcards: expand(config['eQTL_dir'] + '{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID.csv', condition = wildcards.condition, Nk = eval(str(wildcards.condition) + '_Nk'))
    output:
        'output/ld/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_{snp}_ld.csv'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_reformat_sigLDbuddies_{snp}.out',
        err = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_reformat_sigLDbuddies_{snp}.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/reformat_LDbuddies.R {input.leads} {input.buddies} {wildcards.snp} {output} 1> {log.out} 2> log.err
        """

rule join_sigLDbuddies:
    input:
        lambda wildcards: expand('output/ld/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_{snp}_ld.csv', condition = wildcards.condition, Nk = wildcards.Nk, snp = eval(str(wildcards.condition) + '_snps'))
    output:
        'output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD.csv'
    log:
        out = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_joinLD.out',
        err = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_joinLD.err'
    run:
        
        all_leads = []
        for file in input:
            data = pd.read_csv(file)
            all_leads.append(data)

        final_data = pd.concat(all_leads)
        final_data.to_csv(output[0], index = False)

rule get_LDbuddy_rsIDs:
    input:
        rules.join_sigLDbuddies.output
    output:
        'output/qtl/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD_rsID.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_get_LDbuddy_rsIDs.out',
        err = 'output/logs/{condition}_PEER_k{Nk}_genoPC' + fileExt + '_get_LDbuddy_rsIDs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'ld' {output} 1> {log.out} 2> {log.err}
        """

rule make_leadList:
    input:
       lambda wildcards: expand(config['eQTL_dir'] + '{condition}_PEER_k{Nk}_genoPC' + fileExt + '_perm1Mb_sig_rsID.csv', condition = 'FNF',Nk = eval('FNF_Nk'))
    output:
        temp('output/reQTL/FNF_variants.list')
    log:
        out = 'output/logs/make_leadList.out',
        err = 'output/logs/make_leadList.err'
    shell:
        """
        cut -d "," -f2 {input} > output/reQTL/FNF_variants.list
        sed -i '1d' output/reQTL/FNF_variants.list
        """

# Filter the VCF file for significant lead variants (based on variantID for matching)
rule subsetVCF_leadvar:
    input:
        leadvar_FNF = 'output/reQTL/FNF_variants.list',
        vcf = vcf
    output:
        'output/reQTL/FNF_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '_leadVars.vcf.gz'
    params:
        version = config['gatkVersion']
    log:
        out = 'output/logs/subsetVCF_leadvar.out',
        err = 'output/logs/subsetVCF_leadvar.err'
    shell:
        """
        module load gatk/{params.version}
        
        # Subset vcf for these
        gatk SelectVariants -V {input.vcf} --keep-ids {input.leadvar_FNF} -O {output} 1> {log.out} 2> {log.err}
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
        eGene = 'output/qtl/FNF_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '_perm1Mb_sig_rsID.csv'
    output:
        rds = 'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '.rds',
        csv = 'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '.csv'
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

rule join_reQTL_LD:
    input:
        reQTLs = rules.get_reQTLs.output.csv,
        ld = 'output/qtl/FNF_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD_rsID.csv'
    output:
        'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(FNF_Nk) + '_genoPC' + fileExt + '_LD_rsID.csv'
    params:
        version = config['Rversion']
    log:
        out = "output/logs/join_reQTL_LD.out",
        err = "output/logs/join_reQTL_LD.err"
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/join_reQTL_LD.R {input.reQTLs} {input.ld} {output} 1> {log.out} 2> {log.err}
        """
        

