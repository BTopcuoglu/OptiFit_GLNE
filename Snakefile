# Snakefile
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan


# Purpose: # Snakemake file for running mothur 16S pipeline with Leave-One-Out for OptiFit

# # Path to config
# configfile: "config/config.yaml"

# Function for aggregating list of sample numbers.
numSamples = [line.rstrip('\n') for line in open('data/sample_names.txt')]

# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
	input:
		"test.txt"
	shell:
		"""
		mkdir -p logs/mothur/
		mv mothur*logfile logs/mothur/
		"""








##################################################################
#
# Part 1: Grab data and generate sample files
#
##################################################################


# Rule to get glne007.files ans sra data from Neil

# Rule to generate sample_names.txt file from glne007.files file
# awk '{print $1}' data/glne007.files >> "data/sample_names.txt"

numSamples = [line.rstrip('\n') for line in open('data/sample_names.txt')]

##################################################################
#
# Part 2: Generate Reference and Mock Control Files
#
##################################################################

rule getSequences:
	input:
		scipt:"code/bash/getFiles.sh"
	output:
		?


# Downloading and formatting SILVA and RDP reference databases. The v4 region is extracted from
# SILVA database for use as reference alignment.
rule get16SReferences:
	input:
		script="code/bash/mothurReferences.sh"
	output:
		silvaV4="data/mothur/references/silva.v4.align",
		rdpFasta="data/mothur/references/trainset14_032015.pds.fasta",
		rdpTax="data/mothur/references/trainset14_032015.pds.tax"
	conda:
		"envs/glne.yaml"
	shell:
		"bash {input.script}"


##################################################################
#
# Part 3: Generate Contigs for all the samples
#
##################################################################

# Generating master OTU shared file.
rule makeContigs:
	input:
		script="code/bash/mothurContigs.sh",
		refs=rules.get16SReferences.output,
		files_file="data/glne007.files"
	output:
	#You'll want these to be the names of the output files for the files used for leave one out
		sample_shared=expand("data/process/{num}.sample.shared", num=numSamples),
		all_shared=expand("data/process/all_but_{num}.subsampled.shared", num=numSamples)
	conda:
		"envs/glne.yaml"
	shell:
		"bash {input.script} data/process/baxter/ {input.files_file} {input.refs}"



rule leaveOneOut:
	input:
		script:"code/bash/mothurLOO.sh"



rule OptiFIt:
	input:
		script:"code/bash/mothurOptiFit.sh"


# Load R and Rtidyverse modules
rule Model:
	input:
		Rscript code/learning/main.R "L2_Logistic_Regression" "dx" {num}

##################################################################
#
# Part 6: Cleaning
#
##################################################################

# Resets directory by deleting all files created by this workflow.
rule clean:
	shell:
		"""
		echo PROGRESS: Removing all workflow output.
		rm -rf data/references/ data/process/
		"""
