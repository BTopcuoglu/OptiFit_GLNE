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
		echo done
		# mkdir -p logs/mothur/
		# mv mothur*logfile logs/mothur/
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


# Using SRA Run Selector RunInfo table to create mothur files file
rule makeMothurFilesFile:
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
		seqs=readNames
	output:
		"test.txt"
	shell:
		"echo {input.seqs}; touch test.txt"

# # Generating master OTU shared file.
# rule makeContigs:
# 	input:
# 		script="code/bash/mothurContigs.sh",
# 		refs=rules.get16SReferences.output,
# 		files_file="data/glne007.files",
# 		seqs=readNames
# 	output:
# 	#You'll want these to be the names of the output files for the files used for leave one out
# 		sample_shared=expand("data/process/{num}.sample.shared", num=numSamples),
# 		all_shared=expand("data/process/all_but_{num}.subsampled.shared", num=numSamples)
# 	conda:
# 		"envs/glne.yaml"
# 	shell:
# 		"bash {input.script} data/process/baxter/ {input.files_file} {input.refs}"



# rule leaveOneOut:
# 	input:
# 		script:"code/bash/mothurLOO.sh"



# rule OptiFIt:
# 	input:
# 		script:"code/bash/mothurOptiFit.sh"


# # Load R and Rtidyverse modules
# rule Model:
# 	input:
# 		Rscript code/learning/main.R "L2_Logistic_Regression" "dx" {num}

# ##################################################################
# #
# # Part 6: Cleaning
# #
# ##################################################################

# # Resets directory by deleting all files created by this workflow.
# rule clean:
# 	shell:
# 		"""
# 		echo PROGRESS: Removing all workflow output.
# 		rm -rf data/references/ data/process/
# 		"""
