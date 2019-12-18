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

# This script takes a file called SraRunTable.txt, which is generated from hitting the "RunInfo Table"
# button that can be found on the SRA page https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP062005.
# The output will be a directory ($OUTDIR) full of paired fastq.gz files for downstream processing


# Need to get the column that contains the run file names that start wtih SRR (SRA). The
# column position is not consistent across datasets, so need to remove column header and use sed to
# isolate column
SRRS=$(tail -n +2 "${SRAINFO}" | sed 's:.*\(SRR[0-9]*\).*:\1:')


# Loop through each run file and pull it down from the SRA. After downloaded, we want to split it into
# the R1 and R2 files. Finally, we'll compress the files with gzip
for SAMPLE in $(echo "${SRRS}")
do
	prefetch "${SAMPLE}"
    fastq-dump --split-files -O "${OUTDIR}" --gzip "${SAMPLE}"
done


# Some SRR files only contain data for one sequence read. So there aren't problems down the road, we
# want to  make sure all files hav both reads, remove those with only one read
SINGLE_FILES=$(find "${OUTDIR}" -name "*fastq.gz" | cut -f 1 -d _ | sort | uniq -u)

# If $SINGLE_FILES is set (not empty or ""), remove those files
if [ -n "${SINGLE_FILES}" ]
then
	rm "${SINGLE_FILES}"
fi

# Cleaning up SRA temp directories (puts them in current working directory)
find ./ -type d -maxdepth 1 -regex ".*/SRR.*" -exec rm -r {} \;
