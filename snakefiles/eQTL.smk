#!/usr/bin/env python3

# Parse minor allele filters
filterNum = config['minorAllele'].split(":")[0]
filterType = config['minorAllele'].split(":")[1]
if filterType == 'freq':
    filterFlag = '-q'

elif filterType == 'count':
    filterFlag = '-c'


if config['iteratePEER'] == 'TRUE':
    factors = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}.txt', Nk = n) for n in range(1, Nk + 1, config['iterateBy'])],
    var = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}_variance.txt', Nk = n) for n in range(1, Nk + 1, config['iterateBy'])]
else:
    include: "../rules/PEER_kneedle.smk"
    factors = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}.txt', Nk = Nk)],
    var = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}_variance.txt', Nk = Nk)]
    
include: "VCFprocessing.smk"

include: "RNAprocessing.smk"

# Rename geno sample IDs in VCF to donor IDs
rule renameVCFdonors:
    input:
        vcf = rules.updateConfig.output.v,
        index = rules.updateConfig.output.i
    output:
        vcf = 'output/vcf/' + vcf_prefix + '_rename.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_rename.vcf.gz.tbi'
    params:
        samtoolsVersion = config['samtoolsVersion'],
        pythonVersion = config['pythonVersion'],
        prefix = vcf_prefix
    log:
        pythonOut = 'output/logs/renameVCFdonors_py.out',
        pythonErr = 'output/logs/renameVCFdonors_py.err',
        bcftoolsErr = 'output/logs/renameVCFdonors_bcftools.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        module load python/{params.pythonVersion}
        bcftools query -l {input.vcf} > donors.txt
        python3 scripts/renameVCFdonors.py donors.txt 1> {log.pythonOut} 2> {log.pythonErr}
        bcftools reheader -s samples.txt {input.vcf} > output/vcf/{params.prefix}_rename.vcf 2> {log.bcftoolsErr}
        bgzip output/vcf/{params.prefix}_rename.vcf && tabix -p vcf {output.vcf}
        rm donors.txt
        rm samples.txt
        """

# Subset VCF to only include genotypes from donors in samplesheet
rule subsetVCFdonors:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        index = rules.renameVCFdonors.output.index
    output:
        vcf = 'output/vcf/' + vcf_prefix + '_subset.vcf'
    params:
        donors = ",".join(samples['Donor'].unique().tolist()),
        samtoolsVersion = config['samtoolsVersion'],
        pythonVersion = config['pythonVersion'],
        prefix = vcf_prefix
    log:
        pythonOut = 'output/logs/subsetVCFdonors_py.out',
        pythonErr = 'output/logs/subsetVCFdonors_py.err',
        bcftoolsErr = 'output/logs/subsetVCFdonors_bcftools.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        module load python/{params.pythonVersion}
        bcftools query -l {input.vcf} > donors.txt
        python3 scripts/subsetVCFdonors.py donors.txt {params.donors} 1> {log.pythonOut} 2> {log.pythonErr}
        bcftools view -S subset.txt {input.vcf} > output/vcf/{params.prefix}_subset.vcf 2> {log.bcftoolsErr}
        rm donors.txt
        rm subset.txt
        """

## Rule for MAF filtering and het num filtering on VCF
rule filterVCFvariants:
    input:
        vcf = rules.subsetVCFdonors.output.vcf
    output:
        vcf = 'output/vcf/' + vcf_prefix + '_qtl.recode.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_qtl.recode.vcf.gz.tbi'
    params:
        samtoolsVersion = config['samtoolsVersion'],
        gatkVersion = config['gatkVersion'],
        vcftoolsVersion = config['vcftoolsVersion'],
        filterFlag = filterFlag,
        filterNum = filterNum,
        minHets = config['minHets'],
        prefix = vcf_prefix
    log:
        mafFilter_err = 'output/logs/mafFilter.err',
        hetFilter_out = 'output/logs/hetFilter.out',
        hetFilter_err = 'output/logs/hetFilter.err',
        vcfFilter_out = 'output/logs/vcfFilter.out',
        vcfFilter_err = 'output/logs/vcfFilter.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        module load gatk/{params.gatkVersion}
        module load vcftools/{params.vcftoolsVersion}
        bcftools view {params.filterFlag} {filterNum}:minor {input.vcf} > output/vcf/{params.prefix}_minorAlleleFilter.vcf 2> {log.mafFilter_err}
        gatk VariantFiltration -V output/vcf/{params.prefix}_minorAlleleFilter.vcf -O output/vcf/{params.prefix}_hetFilter.vcf.gz --filter-expression "vc.getHetCount() >= {params.minHets}" --filter-name "minhets" 1> {log.hetFilter_out} 2> {log.hetFilter_err}
        vcftools --gzvcf output/vcf/{params.prefix}_hetFilter.vcf.gz --keep-filtered minhets --recode --out output/vcf/{params.prefix}_qtl 1> {log.vcfFilter_out} 2> {log.vcfFilter_err}
        bgzip output/vcf/{params.prefix}_qtl.recode.vcf && tabix -p vcf {output.vcf}
        """

# bcftools stats on final VCF
rule bcftools_stats:
    input:
        vcf = rules.filterVCFvariants.output.vcf,
        index = rules.filterVCFvariants.output.index
    output:
        'output/qc/' + vcf_prefix + '_qtl_stats.txt'
    params:
        version = config['samtoolsVersion']
    shell:
        """
        module load samtools/{params.version}
        bcftools stats {input.vcf} > {output}
        """

rule align:
    input:
        R1 = rules.trim.output.trim1,
        R2 = rules.trim.output.trim2,
        vcf = rules.filterVCFvariants.output.vcf
    output:
        bam = 'output/{group}/align/{group}.Aligned.sortedByCoord.out.bam',
        log = 'output/{group}/align/{group}.Log.final.out'
    threads: 8
    log:
        out = "output/logs/align_{group}.out",
        err = "output/logs/align_{group}.err"
    params:
        genomeDir = config['genomeDir'],
        starVersion = config['starVersion']
    shell:
        'module load star/{params.starVersion} &&'
        'mkdir -p output/{wildcards.group}/align &&'
        'star --runThreadN {threads} '
        '--genomeDir {params.genomeDir} '
        '--readFilesCommand zcat ' 
        '--readFilesIn {input.R1} {input.R2} '
        '--outFileNamePrefix output/{wildcards.group}/align/{wildcards.group}. ' 
        '--outSAMtype BAM SortedByCoordinate '
        '--outFilterType BySJout '
        '--outFilterMultimapNmax 20 ' 
        '--alignSJoverhangMin 8 ' 
        '--alignSJDBoverhangMin 1 '
        '--outFilterMismatchNmax 999 ' 
        '--outFilterMismatchNoverReadLmax 0.04 ' 
        '--alignIntronMin 20 ' 
        '--alignIntronMax 1000000 '
        '--alignMatesGapMax 1000000 '
        '--waspOutputMode SAMtag '
        '--varVCFfile <(zcat {input.vcf}) 1> {log.out}'

rule index:
    input:
        rules.align.output.bam
    output:
        bai = 'output/{group}/align/{group}.Aligned.sortedByCoord.out.bam.bai'
    threads: 8
    params:
        samtoolsVersion = config['samtoolsVersion']
    log:
        out = "output/logs/index_{group}.out",
        err = "output/logs/index_{group}.err"
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools index -@ {threads} {input} {output} 1> {log.out} 2> {log.err}
        """

# Add read groups to bam files for verifybamid, using donor name as group
rule addReadGroups:
    input:
        bam = rules.align.output.bam,
        bai = rules.index.output.bai
    output:
        bam = 'output/{group}/align/{group}.Aligned.sortedByCoord.RG.bam',
        bai = 'output/{group}/align/{group}.Aligned.sortedByCoord.RG.bam.bai'
    params:
        picardVersion = config['picardVersion']
    log:
        err1 = 'output/logs/addReadGroups_{group}.err',
        err2 = 'output/logs/addReadGroups_{group}_index.err'
    shell:
        """
        module load picard/{params.picardVersion}
        module load samtools

        # Parse group for donor name to add to read group
        IFS="_" read -r -a array <<< {wildcards.group}
        donor=${{array[1]}}

        picard AddOrReplaceReadGroups -I {input.bam} -O {output.bam} --RGSM ${{donor}} --RGPL ILLUMINA --RGLB lib1 --RGPU unit1 2> {log.err1}
        # Index
        samtools index -@ {threads} {output.bam} {output.bai} 2> {log.err2}
        """

# rule verifybamid:
#     input:
#         vcf = rules.filterVCFvariants.output.vcf,
#         bam = rules.addReadGroups.output.bam,
#         bai = rules.addReadGroups.output.bai
#     output:
#         selfSM = 'output/qc/{group}_verifybamid.selfSM',
#         selfRG = 'output/qc/{group}_verifybamid.selfRG',
#         bestRG = 'output/qc/{group}_verifybamid.bestRG',
#         bestSM = 'output/qc/{group}_verifybamid.bestSM',
#         depthRG = 'output/qc/{group}_verifybamid.depthRG',
#         depthSM = 'output/qc/{group}_verifybamid.depthSM'
#     params:
#         verifybamid = config['verifybamid']
#     log:
#         err = 'output/logs/verifybamid_{group}.err'
#     shell:
#         """
#         {params.verifybamid} --vcf {input.vcf} --bam {input.bam} --best --out output/qc/{wildcards.group}_verifybamid 2> {log.err}
#         """

rule quant:
    input:
        trim1 = rules.trim.output.trim1,
        trim2 = rules.trim.output.trim2
    output:
        "output/quant/{group}/quant.sf"
    params:
        version = config['salmonVersion'],
        index = config['salmon'],
        gcFlag = config['gcBias'],
        seqFlag = config['seqBias']
    log:
        out = 'output/logs/quant_{group}.out',
        err = 'output/logs/quant_{group}.err'
    shell:
        """
        module load salmon/{params.version}

        if [ {params.gcFlag} == "TRUE" ] && [ {params.seqFlag} == "TRUE" ]; then
            salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 --seqBias --gcBias 1> {log.out} 2> {log.err}
        elif [ {params.gcFlag} == "TRUE" ] && [ {params.seqFlag} != "TRUE" ]; then
            salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 --gcBias 1> {log.out} 2> {log.err}
        elif [ {params.gcFlag} != "TRUE"] && [ {params.seqFlag} == "TRUE" ]; then
            salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 --seqBias 1> {log.out} 2> {log.err}
        else
            salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 1> {log.out} 2> {log.err}
        fi
        """

# rule multiqc:
#     input:
#         [expand('output/qc/{group}_{R}_fastqc.zip', group = key, R = ['R1', 'R2']) for key in read1],
#         [expand('output/{group}/trim/{group}_{R}.fastq.gz_trimming_report.txt', group = key, R =['R1', 'R2']) for key in read1],
#         [expand('output/{group}/align/{group}.Log.final.out', group = key) for key in read1],
#         [expand('output/quant/{group}/quant.sf', group = key) for key in read1],
#         [expand('output/qc/{group}_verifybamid.selfSM', group = key) for key in read1],
#         rules.bcftools_stats.output
#     output:
#         'output/qc/multiqc_report.html'
#     params:
#         version = config['multiqcVersion']
#     log:
#         out = 'output/logs/multiqc.out',
#         err = 'output/logs/multiqc.err'
#     shell:
#         """
#         module load multiqc/{params.version}
#         multiqc -f -o output/qc . 1> {log.out} 2> {log.err}
#         """
    
rule quantNorm_Condition:
    input:
        [expand("output/quant/{group}/quant.sf", group = key) for key in read1]
    output:
        "output/normquant/{condition}_CPMadjTMM_invNorm.bed"
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/{condition}_quantNorm.out',
        err = 'output/logs/{condition}_quantNorm.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/quantNorm.R {params.samplesheet} {wildcards.condition} {input} 1> {log.out} 2> {log.err}
        """

rule quantNorm_All:
    input:
        [expand("output/quant/{group}/quant.sf", group = key) for key in read1]
    output:
        "output/normquant/ALL_CPMadjTMM_invNorm.bed"
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/ALL_quantNorm.out',
        err = 'output/logs/ALL_quantNorm.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/quantNorm.R {params.samplesheet} ALL {input} 1> {log.out} 2> {log.err}
        """

rule indexQuant_Condition:
    input:
        lambda wildcards: ['output/normquant/{condition}_CPMadjTMM_invNorm.bed'.format(condition=wildcards.condition)]
    output:
        bed = 'output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz',
        index = 'output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz.tbi'
    params:
        version = config['samtoolsVersion']
    log:
        out = 'output/logs/{condition}_indexQuant.out',
        err = 'output/logs/{condition}_indexQuant.err'
    shell:
        """
        module load samtools/{params.version}
        bgzip {input} && tabix -p bed {output.bed} 1> {log.out} 2> {log.err}
        """

rule indexQuant_All:
    input:
        rules.quantNorm_All.output
    output:
        bed = 'output/normquant/ALL_CPMadjTMM_invNorm.bed.gz',
        index = 'output/normquant/ALL_CPMadjTMM_invNorm.bed.gz.tbi'
    params:
         version = config['samtoolsVersion']
    log:
        out = 'output/logs/ALL_indexQuant.out',
        err = 'output/logs/ALL_indexQuant.err'
    shell:
        """
        module load samtools/{params.version}
        bgzip {input} && tabix -p bed {output.bed} 1> {log.out} 2> {log.err}
        """

rule getPEER:
    input:
        lambda wildcards: ['output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition=wildcards.condition)]
    output:
        factors = 'output/covar/{condition}_PEERfactors_k{Nk}.txt',
        var = 'output/covar/{condition}_PEERfactors_k{Nk}_variance.txt'
    params:
        Nk = config['PEERfactors'],
        iteratePEER = config['iteratePEER'],
        iterateBy = config['iterateBy']
    log:
        out = 'output/logs/{condition}_{Nk}_getPEER.out',
        err = 'output/logs/{condition}_{Nk}_getPEER.err'
    shell:
        """
        module load r/4.2.2
        Rscript scripts/PEERfactors.R {input} {wildcards.condition} {wildcards.Nk} FALSE {params.iterateBy} 1> {log.out} 2> {log.err}
        """