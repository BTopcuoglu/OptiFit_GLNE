# Snakefile
# Courtney R Armour
# William L. Close
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running mothur 16S pipeline with 
# OptiFit/OptiClust and comparing prediction performance

# Code for creating list of sample and sequence names after filtering out names of mock samples (not used in this study).
import pandas as pd
import re
data = pd.read_csv("data/metadata/SraRunTable.txt")
names = data["Sample Name"].tolist()
regex = re.compile(r'\d+')
sampleNames = set([i for i in names if regex.match(i)])
sequenceNames = set(data[data["Sample Name"].isin(sampleNames)]["Run"].tolist())

alogrithmNames = ['optifit','opticlust']

# Master rule for controlling workflow. Cleans up mothur log files when complete.
# rule all:
#     input:
#         expand("data/learning/summary/{alogrithm}/model_results.tsv",
#             alogrithm = alogrithmNames),
#         expand("data/learning/summary/{alogrithm}/confusion_matrix.tsv",
#             alogrithm = alogrithmNames),
#         "results/tables/fraction_reads_mapped.tsv"
#     shell:
#         '''
#         if $(ls | grep -q "mothur.*logfile"); then
#             mkdir -p logs/mothur/
#             mv mothur*logfile logs/mothur/
#         fi
#         '''

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
        tax="data/process/precluster/glne.precluster.taxonomy",
        dist="data/process/precluster/glne.precluster.dist"
    conda:
        "envs/mothur.yaml"
    shell:
        "bash {input.script} {input.files} {input.refs}"

##################################################################
#
# Part 3: Generate Data Splits
#
##################################################################

 rule make8020splits:
     input:
     output:
     conda:
     shell: