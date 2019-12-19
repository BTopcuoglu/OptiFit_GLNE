# Snakefile
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

# Purpose: # Snakemake file for running mothur 16S pipeline with Leave-One-Out for OptiFit

# # Path to config
# configfile: "config/config.yaml"

# # Function for aggregating list of sample numbers.
# numSamples = [line.rstrip('\n') for line in open('data/sample_names.txt')]

# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
	input:
		# "test.txt",
		expand("data/process/loo/{sample}/{sample}.in.fasta",
			sample = 2003650)
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
		outList="data/process/loo/{sample}/{sample}.out.list"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.precluster} {params.sample}"



# rule OptiFIt:
# 	input:
# 		script:"code/bash/mothurOptiFit.sh"





##################################################################
#
# Part 4: Running ML Model
#
##################################################################

# # Load R and Rtidyverse modules
# rule Model:
# 	input:
# 		Rscript code/learning/main.R "L2_Logistic_Regression" "dx" {num}





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
