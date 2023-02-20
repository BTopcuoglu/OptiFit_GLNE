# Snakefile

# Courtney R Armour
# William L. Close
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running mothur 16S pipeline with
# OptiFit/OptiClust and comparing machine learning classification performance

# Code for creating list of sample and sequence names after filtering out 
# names of mock samples (not used in this study)
import pandas as pd
import re
data = pd.read_csv("data/metadata/SraRunTable.txt")
names = data["Sample Name"].tolist()
regex = re.compile(r'\d+')
sampleNames = set([i for i in names if regex.match(i)])
sequenceNames = set(data[data["Sample Name"].isin(sampleNames)]["Run"].tolist())

# Variables
# methods to compare
methods = ["opticlust_denovo","optifit_self","optifit_gg","vsearch_denovo","vsearch_gg"]
# number of test/train splits to make
split_nums = range(1,101) 
# subsample size 
subsample_size = 10000
# results to merge
metrics=["performance","prediction","hp"]
#for plotting
group_colors = ["#20639b","#3caea3","#f5ad5b","#ed553b","#a989ba"]

# Master rule for controlling workflow. Cleans up mothur log files when complete.
rule all:
    input:
        expand("results/ml/summary/merged_{metric}.csv",
               metric = metrics),
        "results/tables/pvalues.csv"
    shell:
        '''
        if $(ls | grep -q "mothur.*logfile"); then
            mkdir -p logs/mothur/
            mv mothur*logfile logs/mothur/
        fi
        '''
        
##################################################################
#
# Part 1: Download and Prep Data and References
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

# Downloading and formatting SILVA, RDP, and GreenGenes reference databases. 
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
        

# Preclustering the Greengenes reference sequences
rule precluster_gg_reference:
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
        
# cluster gg reference to OTUS
rule cluster_gg_reference:
    input:
        gg_dist=rules.precluster_gg_reference.output.dist,
        gg_count=rules.precluster_gg_reference.output.count_table
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
        
##################################################################
#
# Part 2: Preclustering
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
# Part 3: OptiClust
#
##################################################################

# DE NOVO
rule cluster_opticlust_denovo:
    input:
        fasta=rules.preclusterSequences.output.fasta,
        count_table=rules.preclusterSequences.output.count_table,
        tax=rules.preclusterSequences.output.tax,
        dist=rules.preclusterSequences.output.dist
    output:
        shared="data/process/opticlust_denovo/opticlust_denovo.opti_mcc.shared"#,
        #lst="data/process/opticlust/shared/glne.precluster.opti_mcc.list"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    params:
        outdir="data/process/opticlust_denovo/",
        subsize=subsample_size
    shell:
        """
        mkdir -p {params.outdir}
        
        mothur "#set.current(outputdir={params.outdir}, processors={resources.ncores}, fasta={input.fasta}, count={input.count_table}, taxonomy={input.tax}, column={input.dist});
            cluster(column=current, count=current);
            make.shared(list=current, count=current, label=0.03);
            sub.sample(shared=current, label=0.03, size={params.subsize});
            rename.file(shared=current,prefix=opticlust_denovo)"
        """

# rename for consistency to use wildcards later
rule rename_opticlust_denovo_shared:
    input:
        shared="data/process/opticlust_denovo/opticlust_denovo.opti_mcc.shared"
    output:
        shared_renamed="data/process/opticlust_denovo/glne.opticlust_denovo.shared"        
    shell:
        "mv {input.shared} {output.shared_renamed}"
        
##################################################################
#
# Part 4: Generate Data Splits
#
##################################################################

# randomly split the data 80% train and 20% test for however many splits 
# specified at the top of this file with the variable 'split_nums'
# need the opticlust_denovo shared file to get sample list
rule generate_splits_8020:
    input:
        script="code/R/generate_splits_8020.R",
        shared=rules.cluster_opticlust_denovo.output.shared
    params:
        n_splits=len(split_nums)
    output:
        splits=expand("data/process/splits/split_{num}.csv", num = split_nums)
    shell:
        "Rscript {input.script} {input.shared} {params.n_splits}"
        
##################################################################
#
# Part 5: OptiFit
#
##################################################################

### SELF REFERENCE
# cluster 80% for reference OTUs
rule cluster_optifit_ref:
    input: #THE ORDER OF THESE MATTERS TO RUN CORRECTLY!
        script="code/bash/mothur_cluster_optifit_ref.sh",
        fasta=rules.preclusterSequences.output.fasta,
        count_table=rules.preclusterSequences.output.count_table,
        tax=rules.preclusterSequences.output.tax,
        dist=rules.preclusterSequences.output.dist,
        split="data/process/splits/split_{num}.csv"
    output:
        shared="data/process/optifit_self/split_{num}/train/glne.optifit_self.shared",
        reffasta="data/process/optifit_self/split_{num}/train/glne.precluster.pick.fasta",
        refdist="data/process/optifit_self/split_{num}/train/glne.precluster.pick.dist",
        reflist="data/process/optifit_self/split_{num}/train/glne.precluster.pick.opti_mcc.list",
        refCount="data/process/optifit_self/split_{num}/train/glne.precluster.pick.count_table"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    params:
        outdir="data/process/optifit_self/",
        subsize=subsample_size
    shell:
        "bash {input} {resources.ncores} {params}"

# fit 20% to 80% reference OTUs
rule fit_optifit_self: 
    input: #THE ORDER OF THESE MATTERS TO RUN CORRECTLY!
        script="code/bash/mothur_fit_optifit_self.sh",
        fasta=rules.preclusterSequences.output.fasta,
        count_table=rules.preclusterSequences.output.count_table,
        tax=rules.preclusterSequences.output.tax,
        dist=rules.preclusterSequences.output.dist,
        split="data/process/splits/split_{num}.csv",
        reffasta=rules.cluster_optifit_ref.output.reffasta,
        refdist=rules.cluster_optifit_ref.output.refdist,
        reflist=rules.cluster_optifit_ref.output.reflist
    output:
        fit="data/process/optifit_self/split_{num}/test/glne.optifit_self.shared"
    resources:
        ncores=12,
        time_min="01:20:00",
        mem_mb=30000
    params:
        outdir="data/process/optifit_self/",
        subsize=subsample_size
    shell:
        "bash {input} {resources.ncores} {params}"
 
# GG REFERENCE
rule cluster_optifit_gg:
    input:
        precluster_fasta=rules.preclusterSequences.output.fasta,
        precluster_count=rules.preclusterSequences.output.count_table,
        precluster_tax=rules.preclusterSequences.output.tax,
        precluster_dist=rules.preclusterSequences.output.dist,
        reffasta=rules.precluster_gg_reference.output.fasta,
        refdist=rules.precluster_gg_reference.output.dist,
        reflist=rules.cluster_gg_reference.output.reflist        
    output:
        shared="data/process/optifit_gg/glne.precluster.fit.optifit_mcc.0.03.subsample.shared"
    params:
        outdir="data/process/optifit_gg/",
        subsize=subsample_size
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

# rename for consistency to use wildcards later  
rule rename_optifit_gg_shared:
    input:
        shared="data/process/optifit_gg/glne.precluster.fit.optifit_mcc.0.03.subsample.shared"
    output:
        shared_renamed="data/process/optifit_gg/glne.optifit_gg.shared"        
    shell:
        "mv {input.shared} {output.shared_renamed}"
    
##################################################################
#
# Part 6: VSEARCH
#
##################################################################

# Running de novo and closed ref clustering with vsearch.
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
        redund='data/process/vsearch_prep/tmp/glne.precluster.redundant.fasta',
        rename='data/process/vsearch_prep/tmp/glne.precluster.redundant.renamed.fasta',
        groups='data/process/vsearch_prep/tmp/glne.precluster.redundant.groups',
        groups_rename='data/process/vsearch_prep/tmp/glne.precluster.redundant.renamed.groups' 
    params:
        outdir="data/process/vsearch_prep/tmp/",
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
        fna="data/process/vsearch_prep/tmp/glne.precluster.ng.fasta",
        count_table='data/process/vsearch_prep/tmp/glne.precluster.ng.count_table',
        dist='data/process/vsearch_prep/tmp/glne.precluster.ng.dist'
    params:
        outdir="data/process/vsearch_prep/tmp/",
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
        fna="data/process/vsearch_prep/tmp/glne.precluster.ng.sorted.fasta",
        uc="data/process/vsearch_prep/tmp/glne.precluster.ng.sorted.uc"
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

### DE NOVO
rule vsearch_de_novo:
    input: 
        query=rules.vsearch_sort.output.fna
    output:
        uc="data/process/vsearch_denovo/glne.vsearch_denovo.uc"
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
        
### GG REFERENCE

#get gg ref for use with vsearch
rule download_gg_97:
    output:
        tar='data/process/vsearch_gg/tmp/gg_13_8_otus.tar.gz',
        fasta='data/process/vsearch_gg/gg.fasta'
    shell:
        """
        source /etc/profile.d/http_proxy.sh  # required for internet on the Great Lakes cluster
        wget -N -P data/process/vsearch_gg/tmp/ ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
        tar -xzvf {output.tar} -C data/process/vsearch_gg/tmp/
        cp data/process/vsearch_gg/tmp/gg_13_8_otus/rep_set/97_otus.fasta data/process/vsearch_gg/gg.fasta
        """

#open ref cluster with vsearch using gg ref
# will remove non ref otus before running models
rule vsearch_open_ref:
    input:
        query=rules.vsearch_sort.output.fna,
        gg_ref=rules.download_gg_97.output.fasta
    output:
        closed='data/process/vsearch_gg/glne.closed.uc',
        unmatched='data/process/vsearch_gg/glne.unmatched.fasta',
        denovo='data/process/vsearch_gg/glne.de_novo.uc'
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
        vsearch --usearch_global {input.query} \
            --db {input.gg_ref} \
            --notmatched {output.unmatched} \
            --uc {output.closed} \
            --threads {resources.ncores} \
            --id {params.perc_identity} \
            --minseqlength {params.min_seq_length} \
            --maxaccepts {params.max_accepts} \
            --maxrejects {params.max_rejects} \
            --wordlength {params.word_length} \
            --strand both
        vsearch --cluster_smallmem {output.unmatched} \
            --usersort \
            --uc {output.denovo} \
            --threads {resources.ncores} \
            --id {params.perc_identity} \
            --minseqlength {params.min_seq_length} \
            --maxaccepts {params.max_accepts} \
            --maxrejects {params.max_rejects} \
            --wordlength {params.word_length} \
            --strand both
        """

### summarize both methods
rule uc_to_list:
    input:
        code='code/py/uc_to_list.py',
        uc="data/process/{methods}/glne.{method_set}.uc"
    output:
        list="data/process/{methods}/glne.{method_set}.list"
    wildcard_constraints:
        methods="vsearch_denovo|vsearch_gg"
    script:
        'code/py/uc_to_list.py'
        
        
# combine the closed and de novo vsearch lists from the vsearch open gg ref
rule combine_open_lists:
    input:
        list_closed='data/process/vsearch_gg/glne.closed.list',
        list_denovo='data/process/vsearch_gg/glne.de_novo.list',
        py_combine='code/R/combine_open_lists.R'
    output:
        list='data/process/vsearch_gg/glne.vsearch_gg.list'
    script:
        'code/R/combine_open_lists.R'

rule vsearch_make_shared:
    input:
        list_file="data/process/{methods}/glne.{methods}.list", 
        count_table=rules.degap.output.count_table
    output:
        shared="data/process/{methods}/glne.{methods}.shared"
    params:
        outdir="data/process/{methods}/",
        method="glne.{methods}",
        subsize=subsample_size
    wildcard_constraints:
        methods="vsearch_denovo|vsearch_gg"
    resources:
        mem_mb=20000,
        ncores=8
    shell:
        """
        mothur "#set.current(outputdir={params.outdir},processors={resources.ncores});
            make.shared(list={input.list_file},count={input.count_table});
            sub.sample(shared=current, size={params.subsize}));
            rename.file(shared=current,prefix={params.method})"
        """

rule sensspec_vsearch:
    input:
        list_file="data/process/{methods}/glne.{methods}.list",
        count_table=rules.degap.output.count_table, # data/process/vsearch_prep/tmp/glne.precluster.ng.count_table
        dist=rules.degap.output.dist # data/process/vsearch_prep/tmp/glne.precluster.ng.dist
    output:
        query_accnos="data/process/{methods}/glne.precluster.ng.accnos",
        list_accnos="data/process/{methods}/glne.{methods}.userLabel.pick.accnos",
        list_file="data/process/{methods}/glne.{methods}.userLabel.pick.list",
        tsv='data/process/{methods}/glne.{methods}.userLabel.pick.sensspec'
    params:
        outdir='data/process/{methods}/',
        label='userLabel',
        cutoff=dist_thresh
    resources:
        mem_mb=20000,
        ncores=8
    wildcard_constraints:
        methods="vsearch_denovo|vsearch_gg"
    shell:
        """
        mothur '#set.current(outputdir={params.outdir},processors={resources.ncores});
            list.seqs(count={input.count_table});
            get.seqs(list={input.list_file}, accnos=current);
            list.seqs(list=current);
            sens.spec(list=current, count=current, column={input.dist}, label={params.label}, cutoff={params.cutoff})
            '
        """
##################################################################
#
# Part 7: Generate data 
#
##################################################################

rule generate_data:
    input:
        script="code/bash/generate_data.sh",
        shared="data/process/{methods}/glne.{methods}.shared",
        split="data/process/splits/split_{num}.csv"
    output:
        train="data/process/{methods}/split_{num}/train/glne.{methods}.shared",
        test="data/process/{methods}/split_{num}/test/glne.{methods}.shared"
    wildcard_constraints:
        methods="opticlust_denovo|optifit_gg|vsearch_denovo|vsearch_gg"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=30000
    params:
        outdir="data/process/{methods}/split_{num}/"
    shell:
        "bash {input.script} {input.shared} {input.split} {params.outdir} {resources.ncores}"

##################################################################
#
# Part 8: Preprocess data
#
##################################################################

rule preproecess_data:
    input:
        script="code/R/ml_preprocess_data.R",
        metadata=rules.download_metadata.output.metadata,
        train="data/process/{methods}/split_{num}/train/glne.{methods}.shared",
        test="data/process/{methods}/split_{num}/test/glne.{methods}.shared"
    output:
        preprocTrain="results/ml/{methods}/preproc_train_split_{num}.csv",
        preprocTest="results/ml/{methods}/preproc_test_split_{num}.csv"
    resources:
        ncores=12,
        time_min="01:00:00",
        mem_mb=50000
    shell:
         "Rscript --max-ppsize=500000 {input.script} {input.metadata} {input.train} {input.test} {resources.ncores}"
    
##################################################################
#
# Part 9: Run models
#
##################################################################

rule run_models:
    input:
        script="code/R/ml_run_model.R",
        train="results/ml/{methods}/preproc_train_split_{num}.csv",
        test="results/ml/{methods}/preproc_test_split_{num}.csv"
    output:
        performance="results/ml/{methods}/performance_split_{num}.csv",
        model="results/ml/{methods}/model_split_{num}.rds",
        prediction="results/ml/{methods}/prediction_split_{num}.csv",
        hp_performance="results/ml/{methods}/hp_split_{num}.csv"
    params:
        model="rf",
        outcome="dx"
    resources: 
        ncores=12,
        time_min="72:00:00",
        mem_mb=50000 
    shell:
        "Rscript --max-ppsize=500000 {input.script} {input.train} {input.test} {params.model} {params.outcome} {resources.ncores}"
        
##################################################################
#
# Part 10: Merge Results
#
##################################################################

rule merge_results:
    input:
        results=expand("results/ml/{method}/{{metric}}_split_{num}.csv",
                       method = methods, num = split_nums)
    output:
        merged_results="results/ml/summary/merged_{metric}.csv"
    script:
        "code/R/merge_results.R"

# # merge MCC scores
rule merge_mcc:
    input:
        "data/process/opticlust_denovo/glne.precluster.opti_mcc.sensspec",
        expand("data/process/optifit_self/split_{num}/train/glne.precluster.pick.opti_mcc.sensspec",
                                           num = split_nums),
        expand("data/process/optifit_self/split_{num}/test/glne.precluster.pick.renamed.fit.optifit_mcc.steps",
                                           num = split_nums),
        "data/process/optifit_gg/glne.precluster.fit.optifit_mcc.sensspec",
        "data/process/vsearch_gg/glne.vsearch_gg.userLabel.pick.sensspec",
        "data/process/vsearch_denovo/glne.vsearch_denovo.userLabel.pick.sensspec"
    output:
        merged_mcc="results/ml/summary/merged_mcc.csv"
    script:
        "code/R/get_mcc.R"
         
rule get_sens_spec:
    input: 
        expand("results/ml/opticlust_denovo/prediction_split_{num}.csv",
               num = split_nums),
        expand("results/ml/optifit_self/prediction_split_{num}.csv",
               num = split_nums),
        expand("results/ml/optifit_gg/prediction_split_{num}.csv",
               num = split_nums),
        expand("results/ml/vsearch_denovo/prediction_split_{num}.csv",
               num = split_nums),
        expand("results/ml/vsearch_gg/prediction_split_{num}.csv",
               num = split_nums)   
    output:
        outfile="results/ml/summary/all_sens_spec.csv"
    script:
        "code/R/get_sens_spec.R"

##################################################################
#
# Part 12: Analysis
#
##################################################################

# calculate Pvalues comparing AUCs between methods
rule calculate_pvalues:
    input:
        merged_performance="results/ml/summary/merged_performance.csv"
    output:
        outfile="results/tables/pvalues.csv"
    resources:
        time_min="00:60:00"
    script:
        "code/R/calculate_pvalues.R"

##################################################################
#
# Part 12: Plots
#
##################################################################

# #view how many times each sample is in train/test set across splits
rule plot8020splits:
    input:
        splits=rules.generate_splits_8020.output.splits
    output:
        split_plot="results/figures/split_8020.png"
    params:
        outdir="results/figures/"
    script:
        "code/R/plot_split_8020.R"
        
# plot mcc scores
# score for OptiFit is the query and reference sequences together
rule plot_mcc:
    input:
        mcc=rules.merge_mcc.output.merged_mcc
    output:
        plot="results/figures/mcc_scores.png"
    params:
        colors=group_colors
    script:
        "code/R/plot_mcc.R"

#plot hyperparameter performance 
rule plotHPperformance:
    input:
        hp="results/ml/summary/merged_hp.csv"
    output:
        hp_plot="results/figures/hp_performance.png"
    params:
        outdir="results/figures/",
        colors=group_colors
    script:
        "code/R/plot_HP_performance.R"
            
rule plotAUROC:
    input:
        perf="results/ml/summary/merged_performance.csv",
        pvals=rules.calculate_pvalues.output.outfile
    output:
        fig_auc="results/figures/avg_auroc.png"
    params:
        colors=group_colors
    script:
        "code/R/plot_auroc.R"
        
rule plotAvgROC:
    input:
        senspec=rules.get_sens_spec.output.outfile
    output:
        fig_roc="results/figures/avg_roc.png"
    params:
        colors=group_colors
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
    
# rule createManuscript:
#     input: 
#         script="submission/manuscript.Rmd",
#         fig1a="exploratory/figures/figure1_a.pdf",
#         fig1b="exploratory/figures/figure1_b.pdf",
#         fig2="exploratory/figures/figure2.pdf",
#         ref="submission/references.bib"
#     output:
#         pdf="submission/manuscript.pdf",
#         tex="submission/manuscript.tex"
#     shell:
#         """
#         R -e 'library(rmarkdown);render("submission/manuscript.Rmd", output_format="all")'
#         """
  
# FIX: has to be run from teh submission location because of figure relative paths
# rule createManuscriptWord:
#     input: 
#         tex="submission/manuscript.tex"
#     shell:
#         """
#         pandoc -s submission/manuscript.tex -o submission/manuscript.docx
#         """


onsuccess:
    print("🎉 completed successfully")

onerror:
    print("💥 errors occured")