#!/usr/bin/env python3
import pandas as pd
import os

## Load config file
configfile: "config/config_RNA_signal.yaml"

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep = ",")

## Convert samplesheet columns to strings
samples = samples.astype(str)

## Concatenate Sequencing_Directory to Read1 and Read2 for full read paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Group Seq_Reps
samples['id'] = samples[['Proj', 'Donor']].agg('_'.join, axis=1) + '_R_' + samples[['Condition', 'Time', 'Tech_Rep']].agg('_'.join, axis=1)

## Extract grouped read1 and read2s
read1 = samples.groupby(['id'])['Read1'].apply(list).to_dict()
read2 = samples.groupby(['id'])['Read2'].apply(list).to_dict()

## Merge according to mergeBy parameter, define merge name (mn)
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

## Build dictionary of merged BAM files
mergeSamples = samples.groupby('mn')['id'].apply(list).to_dict()

# Unique sample names in dictionary
mergeSamples_dedup = dict()
for key in mergeSamples:
	unique_values = list(set(mergeSamples[key]))
	mergeSamples_dedup[key] = unique_values



rule all:
    input:
        [expand("output/{group}/signal/{group}.bw", group = key) for key in read1],
        [expand("output/mergeAlign/{mergeName}_{ext}", mergeName=key, ext=['sorted.bam', 'sorted.bam.bai', 'stats.txt']) for key in mergeSamples_dedup],
        [expand("output/mergeSignal/{mergeName}.bw", mergeName=key) for key in mergeSamples_dedup]


rule signal:
	input:
		bam = lambda wildcards: ['output/{group}/align/{group}.Aligned.sortedByCoord.out.bam'.format(group=wildcards.group)]
	output:
		"output/{group}/signal/{group}.bw"
	log:
		err = 'output/logs/signal_{group}.err',
		out = 'output/logs/signal_{group}.out'
	params:
		version = config['deeptoolsVersion']
	shell:
		"""
		module load deeptools/{params.version};
		bamCoverage -b {input.bam} -o {output} 1> {log.out} 2> {log.err}
        """

rule mergeAlign:
	input:
		bams = lambda wildcards: ["output/{group}/align/{group}.Aligned.sortedByCoord.out.bam".format(group=value) for value in mergeSamples_dedup[wildcards.mergeName]],
		bais = lambda wildcards: ["output/{group}/align/{group}.aligned.sortedByCoord.out.bam.bai".format(group=value) for value in mergeSamples_dedup[wildcards.mergeName]]
	output:
		bam = "output/mergeAlign/{mergeName}_sorted.bam",
		bai = "output/mergeAlign/{mergeName}_sorted.bam.bai",
		stats = "output/mergeAlign/{mergeName}_stats.txt"
	log:
		err = 'output/logs/mergeAlign_{mergeName}.err',
		out = 'output/logs/mergeAlign_{mergeName}.out'
	params:
		version = config['samtoolsVersion']
	shell:
		"""
		module load samtools/{params.version};
		samtools merge {output.bam} {input.bams} 1>> {log.out} 2>> {log.err};
		samtools flagstat {output.bam} > {output.stats} 2>> {log.err};
		samtools index {output.bam} 1>> {log.out} 2>> {log.err}
		"""

rule mergeSignal:
	input:
		bam = rules.mergeAlign.output.bam
	output:
		"output/mergeSignal/{mergeName}.bw"
	log:
		err = 'output/logs/mergeSignal_{mergeName}.err',
		out = 'output/logs/mergeSignal_{mergeName}.out'
	params:
		version = config['deeptoolsVersion']
	shell:
		"""
		module load deeptools/{params.version};
		bamCoverage -b {input.bam} -o {output} 1> {log.out} 2> {log.err}
		"""