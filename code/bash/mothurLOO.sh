#! /bin/bash
# mothurLOO.sh
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
# export GROUPS=${1:?ERROR: Need to define GROUPS.}
export FASTA=${1:?ERROR: Need to define FASTA.}
export COUNT=${2:?ERROR: Need to define COUNT.}
export TAXONOMY=${3:?ERROR: Need to define TAXONOMY.}
export SAMPLE=${4:?ERROR: Need to define SAMPLE.}

# Other variables
export OUTDIR=data/process/loo/
export TMP="${OUTDIR}"/"${SAMPLE}"/ # Tmp dir based on sample name to keep things separate during parallelization

# export NUM=2003650
# export GROUP=data/process/baxter/final/full.groups
# export FASTA=data/process/baxter/final/full.fasta
# export COUNT=data/process/baxter/final/full.count_table
# export TAXONOMY=data/process/baxter/final/full.taxonomy



###############################################
# Generate Individual Shared File for Sample  #
###############################################

# Make output dirs if they don't exist
mkdir -p "${TMP}"/

# Create cluster distance file for individual sample
mothur "#get.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${SAMPLE}", outputdir="${TMP}"/);
	dist.seqs(fasta=current, cutoff=0.03)"

# Renaming outputs of files generated from single sample
for FILE in $(find "${TMP}"/ -regex ".*precluster.*"); do

	# Uses path and suffix of input file to rename output file using $SAMPLE and 'in' to represent the file is for
	# the individual sample only.
	REPLACEMENT=$(echo "${FILE}" | sed "s:\(.*/\).*\.\(.*\):\1"${SAMPLE}".in.\2:")

	# Rename files using new format without hyphens
	mv "${FILE}" "${REPLACEMENT}"

done



########################################################
# Generate shared file for all samples but the one #
########################################################

# Cluster all sequences while leaving out the specified sample
mothur "#remove.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${SAMPLE}", outputdir="${TMP}"/);
	dist.seqs(fasta=current, cutoff=0.03);
	cluster(column=current, count=current);
	make.shared(list=current, count=current, label=0.03)"

# Renaming outputs of files generated after leaving the specified sample out
for FILE in $(find "${TMP}"/ -regex ".*precluster.*"); do

	# Uses path and suffix of input file to rename output file using $SAMPLE and 'out' to represent the file is for
	# all of the other samples after the specified sample has been left out.
	REPLACEMENT=$(echo "${FILE}" | sed "s:\(.*/\).*\.\(.*\):\1"${SAMPLE}".out.\2:")

	# Rename files using new format without hyphens
	mv "${FILE}" "${REPLACEMENT}"

done
