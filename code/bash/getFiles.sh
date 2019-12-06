#! /bin/bash
# getFiles.sh
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Load packages
module load Bioinformatics
module load sratoolkit/2.7.0

# Other variables
export WORKDIR=data/mothur/raw/baxter
export OUTDIR=data/process/baxter


# Making output dir
mkdir -p "${OUTDIR}"

##########################
# Grab files from Github #
##########################

# (Baxter, N et al 2016).  The SRA project ID is SRP062005.
# Download the glne007.files file that has sample names

curl -o "${OUTDIR}"/glne007.files https://raw.githubusercontent.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015/master/data/glne007.files

##########################
# Extract sra files #
##########################

# Convert to fasta files that will be used
for sample in data/mothur/raw/baxter/*.sra
do
	fastq-dump --split-files $sample -O "${OUTDIR}"
done

# Rename files file to reflect where fastq files are
sed -ie 's/SRR/data\/process\/baxter\/SRR/g' data/process/baxter/glne007.files
