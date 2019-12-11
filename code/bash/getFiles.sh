#!/bin/bash

# Load packages
module load Bioinformatics
module load sratoolkit/2.9.6

# Other variables
export DATA=data/mothur/raw/baxter
export OUTDIR=data/process/baxter


# Making output dir
mkdir -p "${OUTDIR}"
mkdir -p "${DATA}"

##########################
# Grab files from Github #
##########################

# (Baxter, N et al 2016).  The SRA project ID is SRP062005.
# Download the glne007.files file that has sample names

curl -o "${OUTDIR}"/glne007.files https://raw.githubusercontent.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015/master/data/glne007.files

curl -o "${DATA}"/sra_info.tsv https://raw.githubusercontent.com/SchlossLab/Schloss_Singletons_XXXXX_2019/master/data/human/sra_info.tsv
##########################
# Extract sra files #
##########################

# The input should be something like data/marine. This script assumes that there is a file called
# $DATA/sra_info.tsv, which is generated from hitting the "RunInfo Table" button that can be found
# on pages like https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP012016 (replace ERP012016 with
# the SRA Study identifier. The output will be a directory ($DATA) full of paired fastq.gz files
# for downstream processing

# Start by cleanining up directory by removing pre-existing fastq files
rm -f "${OUTDIR}"/*fastq*


# Need to get the column that contains the run file names that start wtih SRR (SRA) or ERR (ENA). The
# column position is not consistent across datasets, so need to remove column header and use sed to
# isolate column
SRRS=`tail -n +2 "${DATA}"/sra_info.tsv | sed -E "s/.*(.RR\S*)\t.*/\1/"`


# Loop through each run file and pull it down from the SRA. After downloaded, we want to split it into
# the R1 and R2 files. Finaly, we'll compress the files with gzip
for sample in $SRRS
do
  echo $sample
	prefetch $sample
    fasterq-dump --split-files $sample -O "${OUTDIR}"
done

# Some SRR files only contain data for one sequence read. So there aren't problems down the road, we
# want to  make sure all files hav both reads, remove those with only one read
SINGLE_FILES=`ls $DATA/*fastq | cut -f 1 -d _ | sort | uniq -u | sed -E "s/$/*/"`

if [ $SINGLE_FILES ]
then
rm $SINGLE_FILES
fi


# Rename files file to reflect where fastq files are
sed -ie 's/SRR/data\/process\/baxter\/SRR/g' data/process/baxter/glne007.files
