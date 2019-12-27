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
SEQUENCENAME=${1:?ERROR: Need to define SEQUENCENAME.}

# Other variables
OUTDIR=data/raw/



##################################
# Download and extract sra files #
##################################

# Making output dir
mkdir -p "${OUTDIR}"

# Downloading named SRR sequence pair
prefetch "${SEQUENCENAME}"
fastq-dump --split-files -O "${OUTDIR}" --gzip "${SEQUENCENAME}"

# Cleaning up SRA temp directories (puts them in current working directory)
rm -r "${SEQUENCENAME}"/
