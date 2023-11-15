#!/usr/bin/env python3
import os

## Load config file
configfile: "config/config_colocalization.yaml"
vcf = config['vcf']

pbs_nominal_qtls = config['pbs_nominal_qtls']
fnf_nominal_qtls = config['fnf_nominal_qtls']

rule all:
    input:
        [expand('output/coloc/' + os.path.basename(vcf).split('.')[0] + '_MAFs_chr{chr}.csv', chr = chrom) for chrom in range(1, 23)],
        [expand('output/coloc/' + os.path.basename(fnf_nominal_qtls) + '_MAFs_chr{chr}.csv', chr = chrom) for chrom in range(1,23)],
        [expand('output/coloc/' + os.path.basename(pbs_nominal_qtls) + '_MAFs_chr{chr}.csv', chr = chrom) for chrom in range(1,23)],
        'output/coloc/PBS_ALL_coloc.rda',
        'output/coloc/PBS_EUR_coloc.rda',
        'output/coloc/FNF_noresponse_ALL_coloc.rda',
        'output/coloc/FNF_noresponse_EUR_coloc.rda',
        'output/coloc/FNFresponse_ALL_coloc.rda',
        'output/coloc/FNFresponse_EUR_coloc.rda'

# Get MAFs with PLINK --freq for each variant
rule get_variantMAFs:
    input:
        vcf
    output:
        vcf.replace(".vcf.gz", ".frq")
    log:
        out = 'output/logs/get_variantMAFs.out',
        err = 'output/logs/get_variantMAFs.err'
    shell:
        """
        module load plink
        plink --vcf {input} --freq --out output/vcf/$(basename {output} .frq) 1> {log.out} 2> {log.err}
        """

# Reformat this file to a csv, split by chromosome
rule reformat_variantMAFs:
    input:
        rules.get_variantMAFs.output
    output:
        [expand('output/coloc/' + os.path.basename(vcf).split('.')[0] + '_MAFs_chr{chr}.csv', chr = chrom) for chrom in range(1, 23)]  
    params:
        version = config['pythonVersion'],
        filePrefix = os.path.basename(vcf).split('.')[0]
    log:
        out = 'output/logs/reformat_variantMAFs.out',
        err = 'output/logs/reformat_variantMAFs.err'
    shell:
        """
        module load python/{params.version}
        python3 scripts/colocalization/reformat_variantMAFs.py {input} {params.filePrefix} 1> {log.out} 2> {log.err}
        """

rule join_pbs_nomqtl_MAFs:
    input:
        chr_qtl = pbs_nominal_qtls + '_chr{chr}.csv',
        chr_maf = 'output/coloc/' + os.path.basename(vcf).split('.')[0] + '_MAFs_chr{chr}.csv'
    output:
        'output/coloc/' + os.path.basename(pbs_nominal_qtls) + '_MAFs_chr{chr}.csv'  
    log:
        out = 'output/logs/join_pbs_nomqtl_MAFs_chr{chr}.out',
        err = 'output/logs/join_pbs_nomqtl_MAFs_chr{chr}.err'
    params:
        version = config['pythonVersion'],
        filePrefix = 'output/coloc/' + os.path.basename(pbs_nominal_qtls) + '_MAFs'
    shell:
        """
        module load python/{params.version} 
        python3 scripts/colocalization/join_nomMAFs.py {input.chr_qtl} {input.chr_maf} {params.filePrefix} {wildcards.chr} 1> {log.out} 2> {log.err}
        """

rule join_fnf_nomqtl_MAFs:
    input:
        chr_qtl = fnf_nominal_qtls + '_chr{chr}.csv',
        chr_maf = 'output/coloc/' + os.path.basename(vcf).split('.')[0] + '_MAFs_chr{chr}.csv'
    output:
        'output/coloc/' + os.path.basename(fnf_nominal_qtls) + '_MAFs_chr{chr}.csv'  
    log:
        out = 'output/logs/join_fnf_nomqtl_MAFs_chr{chr}.out',
        err = 'output/logs/join_fnf_nomqtl_MAFs_chr{chr}.err'
    params:
        version = config['pythonVersion'],
        filePrefix = 'output/coloc/' + os.path.basename(fnf_nominal_qtls) + '_MAFs'
    shell:
        """
        module load python/{params.version} 
        python3 scripts/colocalization/join_nomMAFs.py {input.chr_qtl} {input.chr_maf} {params.filePrefix} {wildcards.chr} 1> {log.out} 2> {log.err}
        """

# Separate response out of all FN-f eGenes
rule sep_FNF:
    input:
        fnf_all = config['fnf_lead'],
        fnf_response = config['response_qtls']
    output:
        'output/qtl/' + os.path.basename(config['fnf_lead']) + '_noresponse.csv',
    log:
        out = 'output/logs/sep_FNF.out',
        err = 'output/logs/sep_FNF.err'
    params:
        version = config['Rversion']
    shell:
        """
        module load r/{params.version}
        Rscript scripts/colocalization/separate_response.R {input.fnf_all} {input.fnf_response} {output} 1> {log.out} 2> {log.err}
        """

# Perform colocalizations (PBS, FN-f(noresponse), response)
rule coloc_PBS:
    input:
        lambda wildcards: expand('output/coloc/' + os.path.basename(pbs_nominal_qtls) + '_MAFs_chr{chr}.csv', chr = chrom) for chrom in range(1, 23)
    output:
        gwas_all = 'output/coloc/PBS_ALL_coloc.rda',
        gwas_eur = 'output/coloc/PBS_EUR_coloc.rda'
    log:
        out_all = 'output/logs/coloc_PBS_ALL.out',
        err_all = 'output/logs/coloc_PBS_ALL.err',
        out_eur = 'output/logs/coloc_PBS_EUR.out',
        err_eur = 'output/logs/coloc_PBS_EUR.err'
    params:
        eGenes = config['pbs_lead'],
        version = config['Rversion'],
        filePrefix = 'output/coloc/' + os.path.basename(pbs_nominal_qtls) + '_MAFs',
        eqtlN = config['eqtlN'],
        gwas_path = config['gwas_path'],
        ld_info = config['pbs_ld']
    shell:
        """ 
        module load r/{params.version}

        # ALL and EUR GWAS LD
        Rscript scripts/colocalization/colocalizations.r {params.eGenes} {params.ld_info} {params.gwas_path} ALL {params.filePrefix} {params.eqtlN} {output.gwas_all} 1> {log.out_all} 2> {log.err_all}
        Rscript scripts/colocalization/colocalizations.r {params.eGenes} {params.ld_info} {params.gwas_path} EUR {params.filePrefix} {params.eqtlN} {output.gwas_eur} 1> {log.out_eur} 2> {log.err_eur}
        """

rule coloc_FNF_noresponse:
    input:
        [lambda wildcards: expand('output/coloc/' + os.path.basename(fnf_nominal_qtls) + '_MAFs_chr{chr}.csv', chr = chrom) for chrom in range(1, 23)],
        eGenes = rules.sep_FNF.output
    output:
        gwas_all = 'output/coloc/FNF_noresponse_ALL_coloc.rda',
        gwas_eur = 'output/coloc/FNF_noresponse_EUR_coloc.rda'
    log:
        out_all = 'output/logs/coloc_FNF_noresponse_ALL.out',
        err_all = 'output/logs/coloc_FNF_noresponse_ALL.err',
        out_eur = 'output/logs/coloc_FNF_noresponse_EUR.out',
        err_eur = 'output/logs/coloc_FNF_noresponse_EUR.err'
    params:
        version = config['Rversion'],
        filePrefix = 'output/coloc/' + os.path.basename(fnf_nominal_qtls) + '_MAFs',
        eqtlN = config['eqtlN'],
        gwas_path = config['gwas_path'],
        ld_info = config['fnf_ld']
    shell:
        """ 
        module load r/{params.version}

        # ALL and EUR GWAS LD
        Rscript scripts/colocalization/colocalizations.r {input.eGenes} {params.ld_info} {params.gwas_path} ALL {params.filePrefix} {params.eqtlN} {output.gwas_all} 1> {log.out_all} 2> {log.err_all}
        Rscript scripts/colocalization/colocalizations.r {input.eGenes} {params.ld_info} {params.gwas_path} EUR {params.filePrefix} {params.eqtlN} {output.gwas_eur} 1> {log.out_eur} 2> {log.err_eur}
        """

rule coloc_FNFresponse:
    input:
        lambda wildcards: expand('output/coloc/' + os.path.basename(fnf_nominal_qtls) + '_MAFs_chr{chr}.csv', chr = chrom) for chrom in range(1, 23)
    output:
        gwas_all = 'output/coloc/FNFresponse_ALL_coloc.rda',
        gwas_eur = 'output/coloc/FNFresponse_EUR_coloc.rda'
    log:
        out_all = 'output/logs/coloc_FNFresponse_ALL.out',
        err_all = 'output/logs/coloc_FNFresponse_ALL.err',
        out_eur = 'output/logs/coloc_FNFresponse_EUR.out',
        err_eur = 'output/logs/coloc_FNFresponse_EUR.err'
    params:
        eGenes = config['response_qtls'],
        version = config['Rversion'],
        filePrefix = 'output/coloc/' + os.path.basename(fnf_nominal_qtls) + '_MAFs',
        eqtlN = config['eqtlN'],
        gwas_path = config['gwas_path'],
        ld_info = config['fnf_ld']
    shell:
        """ 
        module load r/{params.version}

        # ALL and EUR GWAS LD
        Rscript scripts/colocalization/colocalizations.r {params.eGenes} {params.ld_info} {params.gwas_path} ALL {params.filePrefix} {params.eqtlN} {output.gwas_all} 1> {log.out_all} 2> {log.err_all}
        Rscript scripts/colocalization/colocalizations.r {params.eGenes} {params.ld_info} {params.gwas_path} EUR {params.filePrefix} {params.eqtlN} {output.gwas_eur} 1> {log.out_eur} 2> {log.err_eur}
        """