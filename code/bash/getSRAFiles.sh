#! /bin/bash
# getSRAFiles.sh
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export SRAINFO=${1:?ERROR: Need to define SRAINFO.}

# Other variables
export OUTDIR=data/raw/



##################################
# Download and extract sra files #
##################################

# Making output dir
mkdir -p "${OUTDIR}"

# The input should be something like data/marine. This script takes a file called
# sra_info.tsv, which is generated from hitting the "RunInfo Table" button that can be found
# on the SRA page https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP062005. The output will 
# be a directory ($OUTDIR) full of paired fastq.gz files for downstream processing

# Start by cleanining up directory by removing pre-existing fastq files
rm -f "${OUTDIR}"/*fastq*


# Need to get the column that contains the run file names that start wtih SRR (SRA) or ERR (ENA). The
# column position is not consistent across datasets, so need to remove column header and use sed to
# isolate column
SRRS=$(tail -n +2 "${SRAINFO}" | sed -E "s/.*(.RR\S*)\t.*/\1/")


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
SINGLE_FILES=$(ls "${OUTDIR}"/*fastq | cut -f 1 -d _ | sort | uniq -u | sed -E "s/$/*/")

# If $SINGLE_FILES is set (not empty or ""), remove those files
if [ -n $SINGLE_FILES ]
then
	rm $SINGLE_FILES
fi
