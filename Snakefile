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
sampleNames = set([i for i in names if regex.match(i)])
sequenceNames = set(data[data["Sample Name"].isin(sampleNames)]["Run"].tolist())

# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
	input:
		"data/learning/summary/confusion_matrix.tsv",
		expand("data/process/opticlust/loo/{sample}/{sample}.in.opti_mcc.0.03.subsample.shared",
			sample = sampleNames)
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

# Download 16S SRA sequences from SRP062005.
rule getSRASequences:
	input:
		script="code/bash/getSRAFiles.sh"
	params:
		sequence="{sequence}"
	output:
		read1="data/raw/{sequence}_1.fastq.gz",
		read2="data/raw/{sequence}_2.fastq.gz"
	conda:
		"envs/sra_tools.yaml"
	shell:
		"bash {input.script} {params.sequence}"


# Retrieve tidied metadata from https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015
rule getMetadata:
	input:
		script="code/bash/getMetadata.sh"
	output:
		metadata="data/metadata/metadata.tsv"
	shell:
		"bash {input.script}"


# Downloading and formatting mothur SILVA and RDP reference databases. The v4 region is extracted from
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
# Part 2: Running Mothur
#
##################################################################

# Using SRA Run Selector RunInfo table to create mothur files file.
rule makeFilesFile:
	input:
		script="code/R/makeFilesFile.R",
		sra="data/metadata/SraRunTable.txt", # Output from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP062005&o=acc_s%3Aa RunInfo
		sequences=expand(rules.getSRASequences.output,
			sequence = sequenceNames)
	output:
		files="data/process/glne.files"
	conda:
		"envs/r.yaml"
	shell:
		"Rscript {input.script} {input.sra} {input.sequences}"


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





##################################################################
#
# Part 3: OptiFit Leave One Out (LOO)
#
##################################################################

# Removing one sample at a time and generating cluster files separately for that sample and for
# the remaining data.
rule leaveOneOut:
	input:
		script="code/bash/mothurLOO.sh",
		precluster=rules.preclusterSequences.output
	params:
		sample="{sample}"
	output:
		inFasta="data/process/optifit/loo/{sample}/{sample}.in.fasta",
		inDist="data/process/optifit/loo/{sample}/{sample}.in.dist",
		inCount="data/process/optifit/loo/{sample}/{sample}.in.count_table",
		outFasta="data/process/optifit/loo/{sample}/{sample}.out.fasta",
		outDist="data/process/optifit/loo/{sample}/{sample}.out.dist",
		outList="data/process/optifit/loo/{sample}/{sample}.out.list",
		looSubShared="data/process/optifit/loo/{sample}/{sample}.out.opti_mcc.0.03.subsample.shared" # Used in ML pipeline
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
		sampleSubShared="data/process/optifit/shared/{sample}/{sample}.optifit_mcc.0.03.subsample.shared" # Used in ML pipeline
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.loo}"





##################################################################
#
# Part 4: OptiClust Leave One Out (LOO)
#
##################################################################

# Clustering all of the samples together using OptiClust and generating subsampled shared file.
rule clusterOptiClust:
	input:
		script="code/bash/mothurOptiClust.sh",
		precluster=rules.preclusterSequences.output
	output:
		subShared="data/process/opticlust/shared/glne.opticlust.opti_mcc.0.03.subsample.shared"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.precluster}"


# Removing one group at a time from the OptiClust-generated shared file.
rule leaveOneOutOptiClust:
	input:
		script="code/bash/mothurOptiClustLOO.sh",
		subShared=rules.clusterOptiClust.output.subShared
	params:
		sample="{sample}"
	output:
		sampleSubShared="data/process/opticlust/loo/{sample}/{sample}.in.opti_mcc.0.03.subsample.shared", # Used in ML pipeline
		looSubShared="data/process/opticlust/loo/{sample}/{sample}.out.opti_mcc.0.03.subsample.shared" # Used in ML pipeline
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.subShared} {params.sample}"





##################################################################
#
# Part 5: Running ML Model
#
##################################################################

# Using OptiFit to cluster the output files from the leave-one-out rule
rule predictDiagnosis:
	input:
		script="code/learning/main.R",
		looSubShared=rules.leaveOneOut.output.looSubShared,
		optifitSubShared=rules.clusterOptiFit.output.sampleSubShared,
		metadata=rules.getMetadata.output.metadata
	params:
		model="L2_Logistic_Regression",
		outcome="dx"
	output:
		cvauc="data/learning/results/cv_results_{sample}.csv",
		prediction="data/learning/results/prediction_results_{sample}.csv"
	conda:
		"envs/r.yaml"
	shell:
		"Rscript {input.script} {input.looSubShared} {input.optifitSubShared} {input.metadata} {params.model} {params.outcome}"


# Collating all ML pipeline results and constructing confusion matrix
rule makeConfusionMatrix:
	input:
		script="code/R/makeConfusionMatrix.R",
		metadata=rules.getMetadata.output.metadata,
		results=expand(rules.predictDiagnosis.output,
			sample = sampleNames)
	params:
		dxDiffThresh=0.05, # Threshold for wanting to investigate health data because prediction scores are too close
		classThresh=0.5 # Threshold for calling normal based on prediction values
	output:
		results="data/learning/summary/model_results.tsv",
		confusion="data/learning/summary/confusion_matrix.tsv"
	conda:
		"envs/r.yaml"
	shell:
		"Rscript {input.script} {input.metadata} {input.results} {params.dxDiffThresh} {params.classThresh}"





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
		rm -rf data/raw/ data/references/ data/process/ data/learning/ data/metadata/metadata.tsv
		"""
