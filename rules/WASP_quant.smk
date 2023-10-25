# Filter tagged WASP reads
rule WASPfilter:
    input:
        R = rules.align.output.bam
    output:
        bam = 'output/{group}/align/{group}.Aligned.sortedByCoord.WASP.bam',
        bai = 'output/{group}/align/{group}.Aligned.sortedByCoord.WASP.bam.bai'
    threads: 2
    log:
        err1 = 'output/logs/WASPfilter_{group}_grep.err',
        err2 = 'output/logs/WASPfilter_{group}_samtoolsView.err',
        err3 = 'output/logs/WASPfilter_{group}_samtoolsIndex.err'
    shell:
        """
        module load samtools
        # Add header
        samtools view -H {input.R} > output/{wildcards.group}/align/{wildcards.group}.Aligned.sortedByCoord.WASP.sam
        
        # Grep for WASP-passing reads
        samtools view {input.R} | grep 'vW:i:1' >> output/{wildcards.group}/align/{wildcards.group}.Aligned.sortedByCoord.WASP.sam 2> {log.err1}

        # Compress and index
        samtools view -bS output/{wildcards.group}/align/{wildcards.group}.Aligned.sortedByCoord.WASP.sam > {output.bam} 2> {log.err2}
        samtools index -@ {threads} {output.bam} {output.bai} 2> {log.err3}
        """


# Add read groups to bam files, using donor name as group
rule addReadGroups:
    input:
        bam = rules.WASPfilter.output.bam,
        bai = rules.WASPfilter.output.bai
    output:
        bam = 'output/{group}/align/{group}.Aligned.sortedByCoord.WASP.RG.bam',
        bai = 'output/{group}/align/{group}.Aligned.sortedByCoord.WASP.RG.bam.bai'
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

# Convert filtered bams back to fastqs
rule convertBams:
    input:
        bam = rules.addReadGroups.output.bam,
        bai = rules.addReadGroups.output.bai
    output:
        R1 = 'output/{group}/align/{group}_WASP.RG_R1.fq',
        R2 = 'output/{group}/align/{group}_WASP.RG_R2.fq'
    params:
        samtoolsVersion = config['samtoolsVersion']
    log:
        out = 'output/logs/convertBams_{group}.out',
        err = 'output/logs/convertBams_{group}.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools bam2fq -1 {output.R1} -2 {output.R2} -n {input.bam} 1> {log.out} 2> {log.err}
        """

# Quant fastqcs with salmon
rule quant:
    input:
        fq1 = rules.convertBams.output.R1,
        fq2 = rules.convertBams.output.R2
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
            salmon quant --writeUnmappedNames -l A -1 {input.fq1} -2 {input.fq2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 --seqBias --gcBias 1> {log.out} 2> {log.err}
        elif [ {params.gcFlag} == "TRUE" ] && [ {params.seqFlag} != "TRUE" ]; then
            salmon quant --writeUnmappedNames -l A -1 {input.fq1} -2 {input.fq2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 --gcBias 1> {log.out} 2> {log.err}
        elif [ {params.gcFlag} != "TRUE"] && [ {params.seqFlag} == "TRUE" ]; then
            salmon quant --writeUnmappedNames -l A -1 {input.fq1} -2 {input.fq2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 --seqBias 1> {log.out} 2> {log.err}
        else
            salmon quant --writeUnmappedNames -l A -1 {input.fq1} -2 {input.fq2} -i {params.index} -o output/quant/{wildcards.group} --threads 2 1> {log.out} 2> {log.err}
        fi
        """