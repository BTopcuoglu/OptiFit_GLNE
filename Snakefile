# Snakefile
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################################################################
#
# Part 1: Grab data and generate sample files
#
##################################################################

# Snakemake file for running mothur 16S pipeline with Leave-One-Out for OptiFit

# Rule to get glne007.files ans sra data from Neil

# Rule to generate sample_names.txt file from glne007.files file
# awk '{print $1}' data/glne007.files >> "data/sample_names.txt"

numSamples = [line.rstrip('\n') for line in open('data/sample_names.txt')]

##################################################################
#
# Part 2: Generate Reference and Mock Control Files
#
##################################################################

# Downloading and formatting SILVA and RDP reference databases. The v4 region is extracted from
# SILVA database for use as reference alignment.
rule get16SReferences:
	input:
		script="code/bash/mothurReferences.sh"
	output:
		silvaV4="data/mothur/references/silva.v4.align",
		rdpFasta="data/mothur/references/trainset16_022016.pds.fasta",
		rdpTax="data/mothur/references/trainset16_022016.pds.tax"
	conda:
		"envs/glne.yaml"
	shell:
		"bash {input.script}"


##################################################################
#
# Part 3: Generate Shared Files
#
##################################################################

# Generating master OTU shared file.
rule make16SShared:
	input:
		script="code/bash/mothurShared.sh",
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
