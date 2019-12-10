#!/bin/bash

# The input should be something like data/marine. This script assumes that there is a file called
# $DATA/sra_info.tsv, which is generated from hitting the "RunInfo Table" button that can be found
# on pages like https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP012016 (replace ERP012016 with
# the SRA Study identifier. The output will be a directory ($DATA) full of paired fastq.gz files
# for downstream processing

DATA=data/mothur/raw/baxter


# Start by cleanining up directory by removing pre-existing fastq files
rm -f $DATA/*fastq*


# Need to get the column that contains the run file names that start wtih SRR (SRA) or ERR (ENA). The
# column position is not consistent across datasets, so need to remove column header and use sed to
# isolate column
SRRS=`tail -n +2 $DATA/sra_info.tsv | sed -E "s/.*(.RR\S*)\t.*/\1/"`


# Loop through each run file and pull it down from the SRA. After downloaded, we want to split it into
# the R1 and R2 files. Finaly, we'll compress the files with gzip
for sample in $SRRS
do
  echo $sample
	prefetch $sample
    fastq-dump --split-files $sample -O $OUTDIR
done


# Rename files file to reflect where fastq files are
sed -ie 's/SRR/data\/process\/baxter\/SRR/g' data/process/baxter/glne007.files
