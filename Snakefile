# Snakefile
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running mothur 16S pipeline with Leave-One-Out for OptiFit and diagnosis prediction


# NOTE: This will work for now but will need to create function to pull sample names from files file
# Can use input function similar to readNames below to populate 'expand()' when aggregating data to run model
# Function for creating list of sample names.
import pandas as pd
sampleNames = pd.read_csv("data/metadata/SraRunTable.txt")["Sample Name"].tolist()


# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
	input:
		# "test.txt",
		expand("data/process/optifit/{sample}/{sample}.optifit_mcc.0.03.subsample.shared",
			sample = sampleNames)
	shell:
		"""
		mkdir -p logs/mothur/
		mv mothur*logfile logs/mothur/
		"""





##################################################################
#
# Part 1: Download Data
#
##################################################################

# Download 16S SRA sequences from SRP062005.
checkpoint getSRASequences:
	input:
		script="code/bash/getSRAFiles.sh",
		sra="data/metadata/SraRunTable.txt"
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
		sra="data/metadata/SraRunTable.txt",
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


# NOTE: Will need to adjust leaveOneOut and clusterOptiFit scripts to deal with samples that
# don't have 10000 reads in them. Use count table?



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
		outSubShared="data/process/loo/{sample}/{sample}.out.opti_mcc.0.03.subsample.shared" # Used in ML pipeline
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

# # Load R and Rtidyverse modules
# rule Model:
# 	input:
# 		Rscript code/learning/main.R "L2_Logistic_Regression" "dx" {num}
# 		metadata=rules.getMetadata.output.metadata,





# ##################################################################
# #
# # Part 5: Cleaning
# #
# ##################################################################

# # Resets directory by deleting all files created by this workflow.
# rule clean:
# 	shell:
# 		"""
# 		echo PROGRESS: Removing all workflow output.
# 		rm -rf data/references/ data/process/
# 		"""
