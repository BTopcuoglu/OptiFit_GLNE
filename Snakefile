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
split_nums = range(1,6) # number of test/train splits to make

# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
    input:
        #expand("data/process/opticlust/split_{num}/{partition}/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared", 
        #       num = split_nums,partition = ('train','test')),
        #expand("data/process/optifit/split_{num}/train/glne.precluster.pick.opti_mcc.0.03.subsample.shared",
        #       num = split_nums),
        #expand("data/process/optifit/split_{num}/test/glne.precluster.pick.subsample.renamed.fit.optifit_mcc.shared",
        #       num = split_nums),
        expand("data/learning/results/opticlust/prediction_results_split_{num}.csv",
               num = split_nums),
        expand("data/learning/results/optifit/prediction_results_split_{num}.csv",
               num = split_nums)
    # shell:
    #     '''
    #     if $(ls | grep -q "mothur.*logfile"); then
    #         mkdir -p logs/mothur/
    #         mv mothur*logfile logs/mothur/
    #     fi
    #     '''

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
# Part 3: Cluster OptiClust
#
##################################################################

# Creating subsampled shared file using all of the samples and OptiClust.
rule clusterOptiClust:
    input:
        script="code/bash/mothurOptiClust.sh",
        precluster=rules.preclusterSequences.output
    output:
        shared="data/process/opticlust/shared/glne.precluster.opti_mcc.0.03.subsample.shared"
    conda:
        "envs/mothur.yaml"
    shell:
        "bash {input.script} {input.precluster}"

##################################################################
#
# Part 4: Generate Data Splits
#
##################################################################

rule make8020splits:
    input:
        script="code/R/generate_splits.R",
        shared=rules.clusterOptiClust.output.shared
    params:
        n_splits=5
    output:
        splits=expand("data/process/splits/split_{num}.csv", num = split_nums)
    conda:
        "envs/R.yaml"
    shell:
        "Rscript {input.script} {input.shared} {params.n_splits}"

#view how many times each sample is in train/test set across splits
rule plot8020splits:
    input:
        script="code/R/view_splits.R",
        splits=rules.make8020splits.output.splits
    output:
        split_plot="analysis/view_splits.pdf"
    conda:
        "envs/R.yaml"
    shell:
        "Rscript {input.script} {input.splits}"
        
##################################################################
#
# Part N: Generate OptiClust data
#
##################################################################

rule generateOptiClustData:
    input: 
        script="code/bash/splitOptiClust.sh",
        shared=rules.clusterOptiClust.output.shared,
        split="data/process/splits/split_{num}.csv"
    output:
        train="data/process/opticlust/split_{num}/train/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared",
        test="data/process/opticlust/split_{num}/test/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared"
    conda:
        "envs/R.yaml"
    shell:
        "bash {input.script} {input.shared} {input.split}"

##################################################################
#
# Part N: Generate OptiFit data
#
##################################################################

# cluster the training datasets
rule clusterOptiFitData:
    input:
        script="code/bash/mothurClusterOptiFit.sh",
        precluster=rules.preclusterSequences.output,
        split="data/process/splits/split_{num}.csv"
    output:
        shared="data/process/optifit/split_{num}/train/glne.precluster.pick.opti_mcc.0.03.subsample.shared",
        reffasta="data/process/optifit/split_{num}/train/glne.precluster.pick.fasta",
        refdist="data/process/optifit/split_{num}/train/glne.precluster.pick.dist",
        reflist="data/process/optifit/split_{num}/train/glne.precluster.pick.opti_mcc.list"
    conda:
        "envs/mothur.yaml"
    shell:
        "bash {input.script} {input.precluster} {input.split}"

# fit test data to the training clusters
rule fitOptiFit:
    input:
        script="code/bash/mothurFitOptiFit.sh",
        precluster=rules.preclusterSequences.output,
        split="data/process/splits/split_{num}.csv",
        reffasta=rules.clusterOptiFitData.output.reffasta,
        refdist=rules.clusterOptiFitData.output.refdist,
        reflist=rules.clusterOptiFitData.output.reflist
    output:
        fit="data/process/optifit/split_{num}/test/glne.precluster.pick.subsample.renamed.fit.optifit_mcc.shared"
    conda:
        "envs/mothur.yaml"
    shell:
        "bash {input.script} {input.precluster} {input.split} {input.reffasta} {input.refdist} {input.reflist}"

##################################################################
#
# Part N: Run Models
#
##################################################################

rule runOptiClustModels:
    input:
        script="code/R/run_model.R",
        metadata=rules.getMetadata.output.metadata,
        train=rules.generateOptiClustData.output.train,
        test=rules.generateOptiClustData.output.test
    params:
        model="glmnet",
        outcome="dx"
    output:
        cv="data/learning/results/opticlust/cv_results_split_{num}.csv",
        model="data/learning/results/opticlust/model_split_{num}.rds",
        prediction="data/learning/results/opticlust/prediction_results_split_{num}.csv",
        preproc_train="data/learning/results/opticlust/preproc_train_split_{num}.csv",
        preproc_test="data/learning/results/opticlust/preproc_test_split_{num}.csv"        
    conda:
        "envs/R.yaml"
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {params.model} {params.outcome}"
        
rule runOptiFitModels:
    input:
        script="code/R/run_model.R",
        metadata=rules.getMetadata.output.metadata,
        train=rules.clusterOptiFitData.output.shared,
        test=rules.fitOptiFit.output.fit
    params:
        model="glmnet",
        outcome="dx"
    output:
        cv="data/learning/results/optifit/cv_results_split_{num}.csv",
        model="data/learning/results/optifit/model_split_{num}.rds",
        prediction="data/learning/results/optifit/prediction_results_split_{num}.csv",
        preproc_train="data/learning/results/optifit/preproc_train_split_{num}.csv",
        preproc_test="data/learning/results/optifit/preproc_test_split_{num}.csv"
    conda:
        "envs/R.yaml"
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {params.model} {params.outcome}"
        
##################################################################
#
# Part N: Analysis
#
##################################################################        

# rule mergePredictionResults:
#     input:
#         data/learning/results/optifit/prediction_results_split_{num}.csv