# Snakefile
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running mothur 16S pipeline with Leave-One-Out for OptiFit and diagnosis prediction



# Code for creating list of sample and sequence names after filtering out names of mock samples (not used in this study).
import pandas as pd
import re
data = pd.read_csv("data/metadata/SraRunTable.txt")
names = data["Sample Name"].tolist()
regex = re.compile(r'\d+')
sampleNames = [i for i in names if regex.match(i)]
sequenceNames = data[data["Sample Name"].isin(sampleNames)]["Run"].tolist()

# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
	input:
		"data/learning/results/confusion_matrix.tsv"
	shell:
		'''
		if $(ls | grep -q "mothur.*logfile"); then
			mkdir -p logs/mothur/
			mv mothur*logfile logs/mothur/
		fi
		'''





##################################################################
#
# Part 1: Download Data
#
##################################################################

# Remove mock samples from run data
rule prepareSRARunTable:
	input:
		sra="data/metadata/SraRunTable.txt" # Output from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP062005&o=acc_s%3Aa RunInfo
	output:
		sraNoMock="data/metadata/SraRunTable_no_mock.txt"
	shell:
		"awk '!/mock/' {input.sra} > {output.sraNoMock}"


# Download 16S SRA sequences from SRP062005.
checkpoint getSRASequences:
	input:
		script="code/bash/getSRAFiles.sh",
		sra=rules.prepareSRARunTable.output.sraNoMock
	output:
		dir=directory("data/raw") # Setting output as directory because output files are unknown (samples with 1 read file are removed)
	conda:
		"envs/sra_tools.yaml"
	shell:
		"bash {input.script} {input.sra}"


# Defining a function that pulls the names of all the SRA sequences from getSRASequences after the checkpoint finishes.
def readNames(wildcards):
    checkpoint_output = checkpoints.getSRASequences.get(**wildcards).output.dir
    return expand("data/raw/{readName}.fastq.gz",
    	readName=glob_wildcards(os.path.join(checkpoint_output, "{readName}.fastq.gz")).readName)


# Retrieve tidied metadata from https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015
rule getMetadata:
	input:
		script="code/bash/getMetadata.sh"
	output:
		metadata="data/metadata/metadata.tsv"
	shell:
		"bash {input.script}"





##################################################################
#
# Part 2: Generate Reference Files
#
##################################################################

# Downloading and formatting SILVA and RDP reference databases. The v4 region is extracted from
# SILVA database for use as reference alignment.
rule get16SReferences:
	input:
		script="code/bash/mothurReferences.sh"
	output:
		silvaV4="data/references/silva.v4.align",
		rdpFasta="data/references/trainset14_032015.pds.fasta",
		rdpTax="data/references/trainset14_032015.pds.tax"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script}"





##################################################################
#
# Part 3: Running Mothur
#
##################################################################

# Using SRA Run Selector RunInfo table to create mothur files file.
rule makeFilesFile:
	input:
		script="code/R/makeFilesFile.R",
		sra=rules.prepareSRARunTable.output.sraNoMock,
		seqs=readNames
	output:
		files="data/process/glne.files"
	conda:
		"envs/r.yaml"
	shell:
		"Rscript {input.script} {input.sra} {input.seqs}"


# Preclustering and preparing sequences for leave one out analysis.
rule preclusterSequences:
	input:
		script="code/bash/mothurPrecluster.sh",
		files=rules.makeFilesFile.output.files,
		refs=rules.get16SReferences.output
	output:
		fasta="data/process/precluster/glne.precluster.fasta",
		count="data/process/precluster/glne.precluster.count_table",
		tax="data/process/precluster/glne.precluster.taxonomy"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.files} {input.refs}"


# Removing one sample at a time and generating cluster files separately for that sample and for
# the remaining data.
rule leaveOneOut:
	input:
		script="code/bash/mothurLOO.sh",
		precluster=rules.preclusterSequences.output
	params:
		sample="{sample}"
	output:
		inFasta="data/process/loo/{sample}/{sample}.in.fasta",
		inDist="data/process/loo/{sample}/{sample}.in.dist",
		inCount="data/process/loo/{sample}/{sample}.in.count_table",
		outFasta="data/process/loo/{sample}/{sample}.out.fasta",
		outDist="data/process/loo/{sample}/{sample}.out.dist",
		outList="data/process/loo/{sample}/{sample}.out.list",
		looSubShared="data/process/loo/{sample}/{sample}.out.opti_mcc.0.03.subsample.shared" # Used in ML pipeline
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.precluster} {params.sample}"


# Using OptiFit to cluster the output files from the leave-one-out rule
rule clusterOptiFit:
	input:
		script="code/bash/mothurOptiFit.sh",
		loo=rules.leaveOneOut.output
	output:
		optifitSubShared="data/process/optifit/{sample}/{sample}.optifit_mcc.0.03.subsample.shared" # Used in ML pipeline
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.loo}"





##################################################################
#
# Part 4: Running ML Model
#
##################################################################

# Using OptiFit to cluster the output files from the leave-one-out rule
rule predictDiagnosis:
	input:
		script="code/learning/main.R",
		looSubShared=rules.leaveOneOut.output.looSubShared,
		optifitSubShared=rules.clusterOptiFit.output.optifitSubShared,
		metadata=rules.getMetadata.output.metadata
	params:
		model="L2_Logistic_Regression",
		outcome="dx"
	output:
		cvauc="data/learning/output/cv_results_{sample}.csv",
		prediction="data/learning/output/prediction_results_{sample}.csv"
	conda:
		"envs/r.yaml"
	shell:
		"Rscript {input.script} {input.looSubShared} {input.optifitSubShared} {input.metadata} {params.model} {params.outcome}"


# Collating all ML pipeline results and constructing confusion matrix
rule makeConfusionMatrix:
	input:
		script="code/R/makeConfusionMatrix.R",
		metadata=rules.getMetadata.output.metadata,
		output=expand(rules.predictDiagnosis.output,
			sample = sampleNames)
	params:
		dxDiffThresh=0.05, # Threshold for wanting to investigate health data because prediction scores are too close
		classThresh=0.5 # Threshold for calling normal based on prediction values
	output:
		results="data/learning/results/model_results.tsv",
		confusion="data/learning/results/confusion_matrix.tsv"
	conda:
		"envs/r.yaml"
	shell:
		"Rscript {input.script} {input.metadata} {input.output} {params.dxDiffThresh} {params.classThresh}"





##################################################################
#
# Part 5: Cleaning
#
##################################################################

# Resets directory by deleting all files created by this workflow.
rule clean:
	shell:
		"""
		echo PROGRESS: Removing all workflow output.
		rm -rf data/raw/ data/references/ data/process/ data/learning/ data/metadata/metadata.tsv data/metadata/SraRunTable_no_mock.txt
		"""
