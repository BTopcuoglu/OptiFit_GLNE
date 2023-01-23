# Snakefile
# Courtney R Armour
# William L. Close
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running mothur 16S pipeline with
# OptiFit/OptiClust and comparing prediction performance

# Code for creating list of sample and sequence names after filtering out 
# names of mock samples (not used in this study)
import pandas as pd
import re
data = pd.read_csv("data/metadata/SraRunTable.txt")
names = data["Sample Name"].tolist()
regex = re.compile(r'\d+')
sampleNames = set([i for i in names if regex.match(i)])
sequenceNames = set(data[data["Sample Name"].isin(sampleNames)]["Run"].tolist())

#algorithms to test
alogrithmNames = ['optifit','opticlust']
# number of test/train splits to make
split_nums = range(1,101) 


# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
    input:
        # expand("data/process/gg/full/split_{num}/train/glne.precluster.fit.optifit_mcc.0.03.subsample.0.03.pick.shared",
        #        num = split_nums),
        # expand("data/process/gg/full/split_{num}/test/glne.precluster.fit.optifit_mcc.0.03.subsample.0.03.pick.shared",
        #        num = split_nums)
        # expand("results/ml/gg_full/preproc_train_split_{num}.csv",
        #        num = split_nums)
        # expand("results/ml/gg_subsample_8000/preproc_train_split_{num}.csv",
        #        num = split_nums)
        #
        # expand("results/ml/gg_subsample_8000/performance_split_{num}.csv",
        #        num = split_nums)
        # expand("data/process/vsearch/split_{num}/test/glne.vsearch.userLabel.pick.shared",
        #        num = split_nums),
        # expand("data/process/vsearch/split_{num}/train/glne.vsearch.userLabel.pick.shared",
        #        num = split_nums)
        # expand("results/ml/vsearch/preproc_test_split_{num}.csv",
        #        num=split_nums)
        expand("results/ml/vsearch/performance_split_{num}.csv",
               num=split_nums)
        # "data/process/opticlust/split_1/train/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared",
        # "data/process/opticlust/split_1/test/glne.precluster.pick.renamed.fit.optifit_mcc.0.03.subsample.shared",
        #"data/process/optifit/split_1/train/glne.precluster.pick.opti_mcc.0.03.subsample.shared",
        #"data/process/optifit/split_1/test/glne.precluster.pick.renamed.fit.optifit_mcc.shared",
        #"results/ml/opticlust/preproc_test_split_1.csv",
        #"results/ml/optifit/preproc_test_split_1.csv",
        #"results/ml/opticlust/cv_results_split_1.csv",
        #"results/ml/optifit/cv_results_split_1.csv",
        #"results/ml/opticlust/performance_split_1.csv",
        #"results/ml/optifit/performance_split_1.csv",
        # expand("data/process/optifit/split_{num}/test/glne.precluster.pick.renamed.fit.optifit_mcc.shared",
        #        num = split_nums)
        # expand("results/ml/opticlust/prediction_results_split_{num}.csv",
        #       num = split_nums),
        # expand("results/ml/optifit/prediction_results_split_{num}.csv",
        #        num = split_nums)        
    shell:
        '''
        if $(ls | grep -q "mothur.*logfile"); then
            mkdir -p logs/mothur/
            mv mothur*logfile logs/mothur/
        fi
        '''

rule results:
    input:
        "data/learning/summary/merged_predictions.csv",
        "data/learning/summary/merged_performance.csv",
        "data/learning/summary/merged_HP.csv",
        "data/learning/summary/merged_MCC.csv",
        "data/learning/summary/all_sens_spec.csv",
        #"results/tables/pct_class_correct.csv",
        #"results/tables/fracNonMapped.csv",
        "results/tables/pvalues.csv",
        #"results/tables/splitTogetherFrequency.csv"
        #"results/tables/opticlust_20_mcc.csv"
        
rule plots: 
    input:
        "results/figures/split_8020.png",
        "results/figures/hp_performance.png",
        "results/figures/avg_auroc.png",
        "results/figures/avg_roc.png"

rule submission:
    input:
        "submission/figures/fig2.png"
        
##################################################################
#
# Part 1: Download Data
#
##################################################################

# Download 16S SRA sequences from SRP062005.
rule download_sra_sequences:
    input:
        script="code/bash/download_sra_files.sh"
    params:
        sequence="{sequence}"
    output:
        read1="data/raw/{sequence}_1.fastq.gz",
        read2="data/raw/{sequence}_2.fastq.gz"
    resources:
        time_min="00:45:00"
    shell:
        "bash {input.script} {params.sequence}"

# Retrieve tidied metadata from https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015
rule download_metadata:
    input:
        script="code/bash/download_metadata.sh"
    output:
        metadata="data/metadata/metadata.tsv"
    shell:
        "bash {input.script}"


# Downloading and formatting SILVA, RDP, and greengenes reference databases. 
rule format_16Sreferences:
    input:
        script="code/bash/mothur_format_16Sreferences.sh",
        ggscript="code/R/format_gg_tax.R"
    output:
        silvaV4="data/references/silva.v4.align",
        silvatax="data/references/silva.v4.tax",
        rdpFasta="data/references/trainset16_022016.pds.fasta",
        rdpTax="data/references/trainset16_022016.pds.tax",
        ggV4="data/references/gg.fasta",
        ggTax="data/references/gg.tax"
    resources:  
        ncores=4,
        mem_mb=10000,
        time_min="01:00:00"
    shell:
        "bash {input.script} {resources.ncores}"

##################################################################
#
# Part 2: Running Mothur
#
##################################################################

# Using SRA Run Selector RunInfo table to create mothur files file.
rule makeFilesFile:
    input:
        script="code/R/make_files_file.R",
        sra="data/metadata/SraRunTable.txt", # Output from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP062005&o=acc_s%3Aa RunInfo
        sequences=expand(rules.download_sra_sequences.output,
            sequence = sequenceNames)
    output:
        files="data/process/glne.files"
    shell:
        "Rscript {input.script} {input.sra} {input.sequences}"

# Preclustering and preparing sequences for analysis.
rule preclusterSequences:
    input:
        script="code/bash/mothur_precluster_seqs.sh",
        files=rules.makeFilesFile.output.files,
        refs=rules.format_16Sreferences.output
    output:
        fasta="data/process/precluster/glne.precluster.fasta",
        count_table="data/process/precluster/glne.precluster.count_table",
        tax="data/process/precluster/glne.precluster.taxonomy",
        dist="data/process/precluster/glne.precluster.dist"
    resources:
        ncores=12,
        time_min="08:00:00",
        mem_mb=10000
    shell:
        "bash {input.script} {input.files} {input.refs} {resources.ncores}"

##################################################################
#
# Part 3: Cluster OptiClust
#
##################################################################

# Creating subsampled shared file using all of the samples and OptiClust 
# open-reference clustering
rule cluster_OptiClust_all:
    input:
        script="code/bash/mothur_cluster_OptiClust.sh",
        precluster=rules.preclusterSequences.output
    output:
        shared="data/process/opticlust/shared/glne.precluster.opti_mcc.0.03.subsample.shared",
        lst="data/process/opticlust/shared/glne.precluster.opti_mcc.list"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    shell:
        "bash {input.script} {input.precluster} {resources.ncores}"

##################################################################
#
# Part X: Cluster GG Closed-ref and fit data 
#
##################################################################


rule precluster_gg:
    input:
        script="code/bash/mothur_precluster_db.sh",
        gg_db=rules.format_16Sreferences.output.ggV4,
        gg_tax=rules.format_16Sreferences.output.ggTax
    output:
        fasta="data/process/precluster/gg.precluster.fasta",
        count_table="data/process/precluster/gg.precluster.count_table",
        dist="data/process/precluster/gg.precluster.dist",
        tax="data/process/precluster/gg.precluster.tax"
    resources:
        ncores=12,
        time_min="03:00:00",
        mem_mb=10000
    shell: 
        "bash {input.script} {input.gg_db} {input.gg_tax} {resources.ncores}"

# cluster gg database to OTUS
rule cluster_gg:
    input:
        gg_dist=rules.precluster_gg.output.dist,
        gg_count=rules.precluster_gg.output.count_table
    output:
        reflist="data/process/gg/gg.precluster.opti_mcc.list"
    params:
        outdir="data/process/gg/"
    resources:
        ncores=12,
        time_min="03:00:00",
        mem_mb=10000
    shell:
        """
        mkdir -p {params.outdir}
        
        mothur "#set.current(outputdir={params.outdir},processors={resources.ncores});
        set.seed(seed=19900205);
        cluster(column={input.gg_dist},count={input.gg_count});
        make.shared(list=current,count=current,label=0.03)"
        """
#fit sequences to gg OTUS to create shared file
rule cluster_fit_gg:
    input:
        precluster_fasta=rules.preclusterSequences.output.fasta,
        precluster_count=rules.preclusterSequences.output.count_table,
        precluster_tax=rules.preclusterSequences.output.tax,
        precluster_dist=rules.preclusterSequences.output.dist,
        reffasta=rules.precluster_gg.output.fasta,
        refdist=rules.precluster_gg.output.dist,
        reflist=rules.cluster_gg.output.reflist        
    output:
        #fit_shared="data/process/gg/shared/glne.precluster.fit.optifit_mcc.0.03.subsample.shared"
        fit_shared="data/process/gg/shared/full/glne.precluster.fit.optifit_mcc.0.03.subsample.shared"
    params:
        outdir="data/process/gg/shared/full/",
        subsize=10000
    resources:
        ncores=12,
        time_min="08:00:00",
        mem_mb=50000
    shell: 
        """
        mkdir -p {params.outdir}
        
        mothur "#set.current(outputdir={params.outdir},processors={resources.ncores});
        cluster.fit(fasta={input.precluster_fasta},count={input.precluster_count},reffasta={input.reffasta}, refcolumn={input.refdist}, reflist={input.reflist},method=open, printref=f);
        make.shared(list=current, count=current, label=0.03);
        sub.sample(shared=current, size={params.subsize})" 
        """
        
#fit sequences to gg OTUS to create shared file
# the smallest has 8180 seqs so subsample to 8000
# rule cluster_fit_gg_subsample:
#     input:
#         precluster_fasta=rules.preclusterSequences.output.fasta,
#         precluster_count=rules.preclusterSequences.output.count_table,
#         precluster_tax=rules.preclusterSequences.output.tax,
#         precluster_dist=rules.preclusterSequences.output.dist,
#         reffasta=rules.precluster_gg.output.fasta,
#         refdist=rules.precluster_gg.output.dist,
#         reflist=rules.cluster_gg.output.reflist        
#     output:
#         fit_shared="data/process/gg/shared/subsample_8000/glne.precluster.fit.optifit_mcc.0.03.subsample.shared"
#     params:
#         outdir="data/process/gg/shared/subsample_8000/",
#         subsize=8000
#     resources:
#         ncores=12,
#         time_min=500,
#         mem_mb=50000
#     shell: 
#         """
#         mkdir -p {params.outdir}
        
#         mothur "#set.current(outputdir={params.outdir},processors={resources.ncores});
#         cluster.fit(fasta={input.precluster_fasta},count={input.precluster_count},reffasta={input.reffasta}, refcolumn={input.refdist}, reflist={input.reflist},method=closed, printref=f);
#         remove.seqs(count=current, accnos=current);
#         make.shared(list=current, count=current, label=0.03);
#         sub.sample(shared=current, size={params.subsize})" 
#         """

##################################################################
#
# Part 4: Generate Data Splits
#
##################################################################

# randomly split the data 80% train and 20% test for however many splits 
# specified at the top of this file with the variable 'split_nums'
rule generate_splits_8020:
    input:
        script="code/R/generate_splits_8020.R",
        shared=rules.cluster_OptiClust_all.output.shared
    params:
        n_splits=len(split_nums)
    output:
        splits=expand("data/process/splits/split_{num}.csv", num = split_nums)
    shell:
        "Rscript {input.script} {input.shared} {params.n_splits}"

##################################################################
#
# Part N: Generate OptiClust data
#
##################################################################

# split the shared file from rules.cluster_OptiClust_all.output.shared
# based on thte splits generated in rule.generate_splits_8020 to produce
# a training and test shared file for each split
rule generate_data_OptiClust:
    input: 
        script="code/bash/generate_data_gg.sh",
        shared=rules.cluster_OptiClust_all.output.shared,
        split="data/process/splits/split_{num}.csv"
    output:
        train="data/process/opticlust/split_{num}/train/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared",
        test="data/process/opticlust/split_{num}/test/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    shell:
        "bash {input.script} {input.shared} {input.split} {resources.ncores}"

##################################################################
#
# Part N: Generate OptiFit data
#
##################################################################

# cluster the 80% training datasets using open reference OptiClust
rule cluster_OptiClust_80:
    input:
        script="code/bash/mothur_cluster_OptiClust_80.sh",
        precluster=rules.preclusterSequences.output,
        split="data/process/splits/split_{num}.csv"
    output:
        shared="data/process/optifit/split_{num}/train/glne.precluster.pick.opti_mcc.0.03.subsample.shared",
        reffasta="data/process/optifit/split_{num}/train/glne.precluster.pick.fasta",
        refdist="data/process/optifit/split_{num}/train/glne.precluster.pick.dist",
        reflist="data/process/optifit/split_{num}/train/glne.precluster.pick.opti_mcc.list",
        refCount="data/process/optifit/split_{num}/train/glne.precluster.pick.count_table"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    shell:
        "bash {input.script} {input.precluster} {input.split} {resources.ncores}"

# fit the 20% test data to the training clusters generated in rules.cluster_OptiClust_80
rule fit_OptiFit_20:
    input:
        script="code/bash/mothur_fit_OptiFit_20.sh",
        precluster=rules.preclusterSequences.output,
        split="data/process/splits/split_{num}.csv",
        reffasta=rules.cluster_OptiClust_80.output.reffasta,
        refdist=rules.cluster_OptiClust_80.output.refdist,
        reflist=rules.cluster_OptiClust_80.output.reflist
    output:
        fit="data/process/optifit/split_{num}/test/glne.precluster.pick.renamed.fit.optifit_mcc.0.03.subsample.shared",
        query_count='data/process/optifit/split_{num}/test/glne.precluster.pick.renamed.count_table',
        list='data/process/optifit/split_{num}/test/glne.precluster.pick.renamed.fit.optifit_mcc.list'
    resources:
        ncores=12,
        time_min="01:20:00",
        mem_mb=30000
    shell:
        "bash {input.script} {input.precluster} {input.split} {input.reffasta} {input.refdist} {input.reflist} {resources.ncores}"


##################################################################
#
# Part N: Generate GG data
#
##################################################################

#split shared file into train and test sets
rule generate_gg_data_full:
    input: 
        script="code/bash/generate_data_gg.sh",
        shared=rules.cluster_fit_gg.output.fit_shared,
        split="data/process/splits/split_{num}.csv"
    output:
        train="data/process/gg/full/split_{num}/train/glne.precluster.fit.optifit_mcc.0.03.subsample.0.03.pick.shared",
        test="data/process/gg/full/split_{num}/test/glne.precluster.fit.optifit_mcc.0.03.subsample.0.03.pick.shared"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    params:
        outdir="data/process/gg/full/split_{num}/"
    shell:
        "bash {input.script} {input.shared} {input.split} {params.outdir} {resources.ncores}"
        
# rule generate_gg_data_subsample_8000:
#     input: 
#         script="code/bash/generate_data_gg.sh",
#         shared=rules.cluster_fit_gg_subsample.output.fit_shared,
#         split="data/process/splits/split_{num}.csv"
#     output:
#         train="data/process/gg/subsample_8000/split_{num}/train/glne.precluster.fit.optifit_mcc.0.03.subsample.0.03.pick.shared",
#         test="data/process/gg/subsample_8000/split_{num}/test/glne.precluster.fit.optifit_mcc.0.03.subsample.0.03.pick.shared"
#     resources:
#         ncores=12,
#         time_min="01:00:00",
#         mem_mb=30000
#     params:
#         outdir="data/process/gg/subsample_8000/split_{num}/"
#     shell:
#         "bash {input.script} {input.shared} {input.split} {params.outdir} {resources.ncores}" 

##################################################################
#
# Part N: Preprocess Data to prepare for ML
#
##################################################################

rule preprocess_OptiClust:
    input:
        script="code/R/ml_preprocess_data.R",
        metadata=rules.download_metadata.output.metadata,
        train=rules.generate_data_OptiClust.output.train,
        test=rules.generate_data_OptiClust.output.test
    output:
        preprocTrain="results/ml/opticlust/preproc_train_split_{num}.csv",
        preprocTest="results/ml/opticlust/preproc_test_split_{num}.csv" 
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=50000
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {resources.ncores}"

rule preprocess_OptiFit:
    input: 
        script="code/R/ml_preprocess_data.R",
        metadata=rules.download_metadata.output.metadata,
        train=rules.cluster_OptiClust_80.output.shared,
        test=rules.fit_OptiFit_20.output.fit
    output:
        preprocTrain="results/ml/optifit/preproc_train_split_{num}.csv",
        preprocTest="results/ml/optifit/preproc_test_split_{num}.csv"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=50000
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {resources.ncores}"

rule preprocess_gg_full:
    input: 
        script="code/R/ml_preprocess_gg.R",
        metadata=rules.download_metadata.output.metadata,
        train=rules.generate_gg_data_full.output.train,
        test=rules.generate_gg_data_full.output.test
    output:
        preprocTrain="results/ml/gg_full/preproc_train_split_{num}.csv",
        preprocTest="results/ml/gg_full/preproc_test_split_{num}.csv"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=50000
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {resources.ncores}"
        
# rule preprocess_gg_subsample_8000:
#     input: 
#         script="code/R/ml_preprocess_gg.R",
#         metadata=rules.download_metadata.output.metadata,
#         train=rules.generate_gg_data_subsample_8000.output.train,
#         test=rules.generate_gg_data_subsample_8000.output.test
#     output:
#         preprocTrain="results/ml/gg_subsample_8000/preproc_train_split_{num}.csv",
#         preprocTest="results/ml/gg_subsample_8000/preproc_test_split_{num}.csv"
#     resources:
#         ncores=12,
#         time_min="01:00:00",
#         mem_mb=50000
#     shell:
#         "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {resources.ncores}"
        
##################################################################
#
# Part N: Run Models
#
##################################################################

rule runOptiClustModels:
    input:
        script="code/R/ml_run_model.R",
        train=rules.preprocess_OptiClust.output.preprocTrain,
        test=rules.preprocess_OptiClust.output.preprocTest
    params:
        model="rf",
        outcome="dx"
    output:
        performance="results/ml/opticlust/performance_split_{num}.csv",
        model="results/ml/opticlust/model_split_{num}.rds",
        prediction="results/ml/opticlust/prediction_results_split_{num}.csv",
        hp_performance="results/ml/opticlust/hp_split_{num}.csv"
    resources: 
        ncores=12,
        time_min="72:00:00",
        mem_mb=50000        
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.train} {input.test} {params.model} {params.outcome} {resources.ncores}"
        
rule runOptiFitModels:
    input:
        script="code/R/ml_run_model.R",
        metadata=rules.download_metadata.output.metadata,
        train=rules.preprocess_OptiFit.output.preprocTrain,
        test=rules.preprocess_OptiFit.output.preprocTest
    params:
        model="rf",
        outcome="dx"
    output:
        performance="results/ml/optifit/performance_split_{num}.csv",
        model="results/ml/optifit/model_split_{num}.rds",
        prediction="results/ml/optifit/prediction_results_split_{num}.csv",
        hp_performance="results/ml/optifit/hp_split_{num}.csv"
    resources: 
        ncores=12,
        time_min="72:00:00",
        mem_mb=50000  
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.train} {input.test} {params.model} {params.outcome} {resources.ncores}"
       
rule runGGmodels_full:
    input:
        script="code/R/ml_run_model.R",
        train=rules.preprocess_gg_full.output.preprocTrain,
        test=rules.preprocess_gg_full.output.preprocTest
    params:
        model="rf",
        outcome="dx"
    output:
        performance="results/ml/gg_full/performance_split_{num}.csv",
        model="results/ml/gg_full/model_split_{num}.rds",
        prediction="results/ml/gg_full/prediction_results_split_{num}.csv",
        hp_performance="results/ml/gg_full/hp_split_{num}.csv"
    resources: 
        ncores=12,
        time_min="72:00:00",
        mem_mb=50000  
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.train} {input.test} {params.model} {params.outcome} {resources.ncores}"

# rule runGGmodels_subsample_8000:
#     input:
#         script="code/R/ml_run_model.R",
#         train=rules.preprocess_gg_subsample_8000.output.preprocTrain,
#         test=rules.preprocess_gg_subsample_8000.output.preprocTest
#     params:
#         model="rf",
#         outcome="dx"
#     output:
#         performance="results/ml/gg_subsample_8000/performance_split_{num}.csv",
#         model="results/ml/gg_subsample_8000/model_split_{num}.rds",
#         prediction="results/ml/gg_subsample_8000/prediction_results_split_{num}.csv",
#         hp_performance="results/ml/gg_subsample_8000/hp_split_{num}.csv"
#     resources: 
#         ncores=12,
#         time_min="72:00:00",
#         mem_mb=50000  
#     shell:
#         "Rscript --max-ppsize=500000 {input.script} {input.train} {input.test} {params.model} {params.outcome} {resources.ncores}"

##################################################################
#
# Part N: VSEARCH
#
##################################################################        
# Running de novo clustering with vsearch.
# adapted from https://github.com/SchlossLab/Sovacool_OptiFit_mSphere_2022/blob/main/subworkflows/4_vsearch/Snakefile
dist_thresh = 0.03
perc_identity = 1 - dist_thresh  # to match mothur's 0.03 dissimilarity threshold

# hard-coded params, same as used in Pat's 2015 PeerJ paper and Kelly's 2022 paper
min_seq_length = 30
max_accepts = 16
max_rejects = 64
word_length = 8 # the default value of wordlength is already 8 but I'm paranoid

# deunique seqs & replace underscores with hyphens in fasta headers
#initial inputs are from the preclusterSequences step
rule deunique:
    input:
        fna=rules.preclusterSequences.output.fasta,
        count_table=rules.preclusterSequences.output.count_table 
    output:
        redund='data/process/vsearch/tmp/glne.precluster.redundant.fasta',
        rename='data/process/vsearch/tmp/glne.precluster.redundant.renamed.fasta',
        groups='data/process/vsearch/tmp/glne.precluster.redundant.groups',
        groups_rename='data/process/vsearch/tmp/glne.precluster.redundant.renamed.groups' 
    params:
        outdir="data/process/vsearch/tmp/",
    resources: 
        ncores=12,
        time_min="00:30:00",
        mem_mb=10000  
    shell:
        """
        mothur '#set.current(outputdir={params.outdir},processors={resources.ncores});
            deunique.seqs(fasta={input.fna}, count={input.count_table})'
            
        sed "s/_/-/g" < {output.redund} > {output.rename}
        sed "s/_/-/g" < {output.groups} > {output.groups_rename}
        """

rule degap:
    input:
        fna=rules.deunique.output.rename,
        group=rules.deunique.output.groups_rename
    output:
        fna="data/process/vsearch/tmp/glne.precluster.ng.fasta",
        count_table='data/process/vsearch/tmp/glne.precluster.ng.count_table',
        dist='data/process/vsearch/tmp/glne.precluster.ng.dist'
    params:
        outdir="data/process/vsearch/tmp/",
        cutoff=dist_thresh
    resources:
        ncores=12,
        time_min="02:00:00",
        mem_mb=10000     
    shell:
        """
        mothur '#set.current(outputdir={params.outdir},processors={resources.ncores});
            unique.seqs(fasta={input.fna});
            count.seqs(name=current,group={input.group});
            dist.seqs(fasta=current, cutoff={params.cutoff});
            degap.seqs(fasta=current);
            rename.file(fasta=current, count=current, column=current, prefix=glne.precluster.ng)'
        """

rule vsearch_sort:
    input:
        fna=rules.degap.output.fna
    output:
        fna="data/process/vsearch/tmp/glne.precluster.ng.sorted.fasta",
        uc="data/process/vsearch/tmp/glne.precluster.ng.sorted.uc"
    resources:
        time_min="00:20:00",
        mem_mb=10000   
    shell:
        """
        vsearch \
            --derep_fulllength {input.fna} \
            --sizeout \
            --minseqlength 30 \
            --threads 1 \
            --uc {output.uc} \
            --output {output.fna} \
            --strand both
        """
        
rule vsearch_de_novo:
    input: 
        query=rules.vsearch_sort.output.fna
    output:
        uc="data/process/vsearch/tmp/glne.vsearch.uc"
    params:
        perc_identity=perc_identity,
        min_seq_length=min_seq_length,
        max_accepts=max_accepts,
        max_rejects=max_rejects,
        word_length=word_length
    resources:
        ncores=12,
        time_min="02:00:00",
        mem_mb=20000 
    shell:
        """
        vsearch --cluster_smallmem {input.query} \
            --usersort \
            --uc {output.uc} \
            --threads {resources.ncores} \
            --id {params.perc_identity} \
            --minseqlength {params.min_seq_length} \
            --maxaccepts {params.max_accepts} \
            --maxrejects {params.max_rejects} \
            --wordlength {params.word_length} \
            --strand both
        """

rule uc_to_list:
    input:
        code='code/py/uc_to_list.py',
        uc=rules.vsearch_de_novo.output.uc
    output:
        list_file="data/process/vsearch/tmp/glne.vsearch.list"
    script:
        'code/py/uc_to_list.py'

rule vsearch_make_shared:
    input:
        list_file=rules.uc_to_list.output.list_file, 
        count_table=rules.degap.output.count_table
    output:
        shared="data/process/vsearch/shared/glne.vsearch.shared"
    params:
        outdir="data/process/vsearch/shared/"
    resources:
        mem_mb=20000 
    shell:
        """
        mothur "#set.current(outputdir={params.outdir},processors={resources.ncores});
            make.shared(list={input.list_file},count={input.count_table})"
        """

# rule sensspec_vsearch:
#     input:
#         list=rules.uc_to_list.output.list,
#         count_table=rules.degap.output.count_table,
#         dist=rules.degap.output.dist
#     output:
#         query_accnos="data/process/vsearch/tmp/glne.precluster.ng.accnos",
#         list_accnos="data/process/vsearch/tmp/glne.vsearch.userLabel.pick.accnos",
#         list="data/process/vsearch/tmp/glne.vsearch.userLabel.pick.list",
#         tsv="data/process/vsearch/tmp/glne.vsearch.userLabel.pick.sensspec"
#     params:
#         outdir="data/process/vsearch/tmp/",
#         label="userLabel",
#         cutoff=dist_thresh
#     resources:
#         ncores=12,
#         time_min="02:00:00",
#         mem_mb=20000 
#     shell:
#         """
#         mothur "#set.current(outputdir={params.outdir},processors={resources.ncores});
#             list.seqs(count={input.count_table});
#             get.seqs(list={input.list}, accnos=current);
#             list.seqs(list=current);
#             sens.spec(list=current, count=current, column={input.dist}, label={params.label}, cutoff={params.cutoff})"
#         """
        
rule generate_data_vsearch:
    input: 
        script="code/bash/generate_data_vsearch.sh",
        shared=rules.vsearch_make_shared.output.shared,
        split="data/process/splits/split_{num}.csv"
    output:
        train="data/process/vsearch/split_{num}/train/glne.vsearch.userLabel.pick.shared",
        test="data/process/vsearch/split_{num}/test/glne.vsearch.userLabel.pick.shared"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    params:
        outdir="data/process/vsearch/split_{num}/"
    shell:
        "bash {input.script} {input.shared} {input.split} {params.outdir} {resources.ncores}"
        
rule preprocess_vsearch:
    input: 
        script="code/R/ml_preprocess_vsearch.R",
        metadata=rules.download_metadata.output.metadata,
        train=rules.generate_data_vsearch.output.train,
        test=rules.generate_data_vsearch.output.test
    output:
        preprocTrain="results/ml/vsearch/preproc_train_split_{num}.csv",
        preprocTest="results/ml/vsearch/preproc_test_split_{num}.csv"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=50000
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {resources.ncores}"
       
rule run_models_vsearch:
    input:
        script="code/R/ml_run_model.R",
        train=rules.preprocess_vsearch.output.preprocTrain,
        test=rules.preprocess_vsearch.output.preprocTest
    params:
        model="rf",
        outcome="dx"
    output:
        performance="results/ml/vsearch/performance_split_{num}.csv",
        model="results/ml/vsearch/model_split_{num}.rds",
        prediction="results/ml/vsearch/prediction_results_split_{num}.csv",
        hp_performance="results/ml/vsearch/hp_split_{num}.csv"
    resources: 
        ncores=12,
        time_min="72:00:00",
        mem_mb=50000  
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.train} {input.test} {params.model} {params.outcome} {resources.ncores}"

##################################################################
#
# Part N: Merge Results
#
##################################################################        

rule quantifySplitTogetherFreq:
    input:
        splits=rules.generate_splits_8020.output.splits
    output:
        splitTogetherFreq="results/tables/splitTogetherFrequency.csv"
    script:
        "code/R/quantifySplitTogether.R"

# rule calcFracNonMapped:
#     input:
#         script="code/R/calc_nonmapped.R",
#         OFshared=expand("data/process/optifit/split_{num}/test/glne.precluster.pick.renamed.fit.optifit_mcc.0.03.subsample.shared",
#                              num = split_nums)
#     output:
#         optifitNonMapped="results/tables/fracNonMapped.csv"
#     shell:
#         "Rscript {input.script} {input.OFshared}"
    
rule mergePredictionResults:
    input:
        script="code/R/merge_predictions.R",
        opticlustPred=expand("results/ml/opticlust/prediction_results_split_{num}.csv",
                             num = split_nums),
        optifitPred=expand("results/ml/optifit/prediction_results_split_{num}.csv",
                           num = split_nums)
    output:
        mergedPrediction="data/learning/summary/merged_predictions.csv",
    shell:
        "Rscript {input.script} {input.opticlustPred} {input.optifitPred}"
        
rule mergePerformanceResults:
    input:
        script="code/R/mergePerformance.R",
        opticlustPerf=expand("results/ml/opticlust/performance_split_{num}.csv",
                           num = split_nums),
        optifitPerf=expand("results/ml/optifit/performance_split_{num}.csv",
                           num = split_nums),
        ggFullPerf=expand("results/ml/gg_full/performance_split_{num}.csv",
                          num = split_nums),
        vsearchPerf=expand("results/ml/vsearch/performance_split_{num}.csv",
                           num = split_nums)
    output:
        mergedPerf="results/ml/summary/merged_performance.csv"
    shell:
        "Rscript {input.script} {input.opticlustPerf} {input.optifitPerf} {input.ggFullPerf} {input.vsearchPerf}"

rule mergeHPperformance:
    input:
        script="code/R/mergeHP.R",
        opticlustHP=expand("results/ml/opticlust/hp_split_{num}.csv",
                           num = split_nums),
        optifitHP=expand("results/ml/optifit/hp_split_{num}.csv",
                         num = split_nums)
    output:
        mergedHP="data/learning/summary/merged_HP.csv"
    shell:
        "Rscript {input.script} {input.opticlustHP} {input.optifitHP}"
        
rule getMCCdata:
    input:
        script="code/R/get_mcc.R",
        opticlustSensspec="data/process/opticlust/shared/glne.precluster.opti_mcc.sensspec",
        optifitTrainSensspec=expand("data/process/optifit/split_{num}/train/glne.precluster.pick.opti_mcc.sensspec",
                                    num = split_nums),
        optifitTestSensspec=expand("data/process/optifit/split_{num}/test/glne.precluster.pick.renamed.fit.optifit_mcc.steps",
                                   num = split_nums)
    output:
        mergedMCC="data/learning/summary/merged_MCC.csv"
    shell:
        "Rscript {input.script} {input.opticlustSensspec} {input.optifitTrainSensspec} {input.optifitTestSensspec}"
         
rule get_sens_spec:
    input: 
        opticlust_pred=expand("results/ml/opticlust/prediction_results_split_{num}.csv",
                              num = split_nums),
        optifit_pred=expand("results/ml/optifit/prediction_results_split_{num}.csv",
                            num = split_nums)   
    output:
        allSensSpec="data/learning/summary/all_sens_spec.csv"
    script:
        "code/R/get_sens_spec.R"
        
# rule get_opticlust_20_mcc: 
#     input:
#         script="code/bash/get_opticlust_20_mcc.sh",
#         dist=rules.preclusterSequences.output.dist,
#         count_table=rules.preclusterSequences.output.count_table,
#         split="data/process/splits/split_{num}.csv"
#     output:
#         sensspec="data/process/opticlust/split_{num}/sub/glne.precluster.pick.opti_mcc.sensspec",
#         shared="data/process/opticlust/split_{num}/sub/glne.precluster.pick.opti_mcc.0.03.subsample.shared"
#     resources:
#         ncores=12,
#         time_min="01:00:00",
#         mem_mb=30000
#     shell:
#         "bash {input.script} {input.dist} {input.count_table} {input.split} {resources.ncores}"

# rule merge_opticlust_20_mcc:
#     input:
#         sensspec=expand("data/process/opticlust/split_{num}/sub/glne.precluster.pick.opti_mcc.sensspec",
#                         num = split_nums)
#     output:
#         opticlust_20_mcc="results/tables/opticlust_20_mcc.csv"
#     script:
#         "code/R/merge_opticlust_20_mcc.R"

# rule calc_pct_correct:
#     input:
#         opticlust_pred=expand("results/ml/opticlust/prediction_results_split_{num}.csv",
#                               num = split_nums),
#         optifit_pred=expand("results/ml/optifit/prediction_results_split_{num}.csv",
#                             num = split_nums)
#     output:
#         pct_correct="results/tables/pct_class_correct.csv"
#     script:
#         "code/R/get_pct_correct.R"    
        
rule calcPvalues:
    input:
        mergedPerf=rules.mergePerformanceResults.output.mergedPerf
    output:
        pvalues="results/tables/pvalues.csv"
    script:
        "code/R/calculate_pvalues.R"
        
##################################################################
#
# Part N: Plots
#
##################################################################        

rule plotHPperformance:
    input:
        hp=rules.mergeHPperformance.output.mergedHP
    output:
        hpPlot="results/figures/hp_performance.png"
    script:
        "code/R/plot_HP_performance.R"

#view how many times each sample is in train/test set across splits
rule plot8020splits:
    input:
        splits=rules.generate_splits_8020.output.splits
    output:
        split_plot="results/figures/split_8020.png"
    script:
        "code/R/plot_split_8020.R"
            
rule plotAUROC:
    input:
        perf=rules.mergePerformanceResults.output.mergedPerf,
        pvals=rules.calcPvalues.output.pvalues
    output:
        fig_auc="results/figures/avg_auroc.png"
    script:
        "code/R/plot_auroc.R"
        
rule plotAvgROC:
    input:
        senspec=rules.get_sens_spec.output.allSensSpec
    output:
        fig_roc="results/figures/avg_roc.png"
    script:
        "code/R/plot_avgROC.R"
        
##################################################################
#
# DOCUMENTS
#
##################################################################                

rule update_pages:
    input: 
        exploratory="exploratory/exploratory.html"
    output:
        pages="docs/exploratory.html"
    shell:
        "cp {input.exploratory} docs/"

rule createFig2:
    input:    
        fig2ab=rules.plotAUROC.output.fig_auc,
        fig2c=rules.plotAvgROC.output.fig_roc
    output:
        fig2="submission/figures/fig2.png"
    script:
        "code/R/create_fig2.R"
    
rule createManuscript:
    input: 
        script="submission/manuscript.Rmd",
        fig1a="exploratory/figures/figure1_a.pdf",
        fig1b="exploratory/figures/figure1_b.pdf",
        fig2="exploratory/figures/figure2.pdf",
        ref="submission/references.bib"
    output:
        pdf="submission/manuscript.pdf",
        tex="submission/manuscript.tex"
    shell:
        """
        R -e 'library(rmarkdown);render("submission/manuscript.Rmd", output_format="all")'
        """
  
# FIX: has to be run from teh submission location because of figure relative paths
# rule createManuscriptWord:
#     input: 
#         tex="submission/manuscript.tex"
#     shell:
#         """
#         pandoc -s submission/manuscript.tex -o submission/manuscript.docx
#         """

##################################################################
#
# CLEAN UP
#
##################################################################                

# remove prior output to force recreation of all files
rule clean:
    shell:
        """
        rm -rf data/raw
        rm -rf data/process
        rm -rf data/learning
        rm analysis/*
        rm logs/mothur/*
        rm logs/slurm/*
        """
        
onsuccess:
    print("🎉 completed successfully")

onerror:
    print("💥 errors occured")