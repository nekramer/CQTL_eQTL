#!/usr/bin/env python3

## Load config file
configfile: "config/config_reQTL_final.yaml"
vcf = config['vcf']
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_ALL_qc.vcf.gz", vcf_file).span()[0]]
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
        'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(Nk) + '_genoPC' + fileExt + '.rds',
        'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(Nk) + '_genoPC' + fileExt + '.csv'

# Get rsIDs for significant lead variants
rule sig_rsIDs:
    input:
        config['eQTL_dir'] + 'FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_sig.csv'
    output:
        'output/qtl/FNF_PEER_k' + str(Nk) + 'genoPC' + fileExt + '_perm1Mb_sig_rsID.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/logs/FNF_PEER_k' + str(Nk) + 'genoPC' + fileExt + 'sig_rsIDs.out',
        err = 'output/logs/FNF_PEER_k' + str(Nk) + 'genoPC' + fileExt + 'sig_rsIDs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'lead' {output} 1> {log.out} 2> {log.err}
        """

# Make LD reference file based on input vcf
rule makeLDref:
    input:
        vcf
    output:
        bed = 'output/vcf/' + vcf_prefix + '.bed',
        bim = 'output/vcf/' + vcf_prefix + '.bim',
        fam = 'output/vcf/' + vcf_prefix + '.fam'
    log:
        out = 'output/vcf/logs/makeLDref.out',
        err = 'output/vcf/logs/makeLDref.err'
    params:
        prefix = vcf_prefix
    shell:
        """
        module load plink
        plink --gzvcf {input} --make-bed --out {params.prefix} 1> {log.out} 2> {log.err}
        """

# Get LD buddies for significant lead variants
rule get_sigLDbuddies:
    input:
        leads = rules.sig_rsIDs.output,
        ldref_bed = rules.makeLDref.output.bed,
        ldref_bim = rules.makeLDref.output.bim,
        ldref_fam = rules.makeLDref.output.fam
    output:
        'output/ld/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + 'sig_rsID_{snp}_ld.ld'
    params:
        filePrefix = filePrefix
    log:
        out = 'output/logs/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_get_sigLDbuddies_{snp}.out',
        err = 'output/logs/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_get_sigLDbuddies_{snp}.err'
    shell:
        """
        module load plink
        # Get file prefix of ldref for plink
        ldref_prefix=${{{input.ldref_bed}%%.*}}

        # Iterate through snps by variantID to match with vcf file
        for snp in `awk -F ',' '{{print $8}}' {input.leads} | awk 'NR!=1 {{print}}'`
        do
            plink --bfile ${{ldref_prefix}} --ld-snp ${{snp}} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out output/ld/{params.filePrefix}_perm1Mb_sig_rsID_${{snp}}_ld 1> {log.out} 2> {log.err}
        done
        """

rule reformat_sigLDbuddies:
    input:
        buddies = rules.get_sigLDbuddies.output,
        leads = rules.sig_rsIDs.output
    output:
        'output/ld/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_sig_rsID_{snp}_ld.csv'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_reformat_sigLDbuddies_{snp}.out',
        err = 'output/logs/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_reformat_sigLDbuddies_{snp}.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/reformat_LDbuddies.R {input.leads} {input.buddies} {wildcards.snp} {output} 1> {log.out} 2> log.err
        """

rule join_sigLDbuddies:
    input:
        lambda wildcards: [expand('output/ld/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_sig_rsID_{snp}_ld.csv', snp = s) for s in pd.read_csv(rules.sig_rsIDs.output, usecols = ['variantID'])['variantID'].values]
    output:
        'output/qtl/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD.csv'
    log:
        out = 'output/logs/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_joinLD.out',
        err = 'output/logs/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_joinLD.err'
    run:
        
        all_leads = []
        for file in input:
            data = pd.read_csv(file)
            all_leads.append(data)

        final_data = pd.concat(all_leads)
        final_data.to_csv(output)

rule get_LDbuddy_rsIDs:
    input:
        rules.join_sigLDbuddies.output
    output:
        'output/qtl/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_sig_rsID_LD_rsID.csv'
    params:
        version = config['pythonVersion'],
        dbSNP_dir = config['dbSNP_dir'],
        dbSNP_prefix = config['dbSNP_prefix'],
        dbSNP_suffix = config['dbSNP_suffix']
    log:
        out = 'output/logs/FNF_PEER_k'+ str(Nk) + '_genoPC' + fileExt + '_get_LDbuddy_rsIDs.out',
        err = 'output/logs/FNF_PEER_k'+ str(Nk) + '_genoPC' + fileExt + '_get_LDbuddy_rsIDs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/get_rsids.py {input} {params.dbSNP_dir} {params.dbSNP_prefix} {params.dbSNP_suffix} 'ld' {output} 1> {log.out} 2> {log.err}
        """

rule make_leadList:
    input:
        rules.sig_rsIDs.output
    output:
        temp('output/reQTL/FNF_variants.list')
    log:
        out = 'output/logs/make_leadList.out',
        err = 'output/logs/make_leadList.err'
    shell:
        """
        cut -d "," -f7 {input} > output/reQTL/FNF_variants.list
        sed -i '1d' output/reQTL/FNF_variants.list
        """

# Filter the VCF file for significant lead variants
rule subsetVCF_leadvar:
    input:
        leadvar_FNF = 'output/reQTL/FNF_variants.list',
        vcf = vcf
    output:
        'output/reQTL/FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_leadVars.vcf.gz'
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
        eGene = rules.sig_rsIDs.output
    output:
        rds = 'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(Nk) + '_genoPC' + fileExt + '.rds',
        csv = 'output/reQTL/FNF_sig_reQTLs_PEER_k' + str(Nk) + '_genoPC' + fileExt + '.csv'
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

