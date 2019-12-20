#! /bin/bash
# mothurOptiFit.sh
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export INFASTA=${1:?ERROR: Need to define INFASTA.}
export INDIST=${2:?ERROR: Need to define INDIST.}
export INCOUNT=${3:?ERROR: Need to define INCOUNT.}
export OUTFASTA=${4:?ERROR: Need to define OUTFASTA.}
export OUTDIST=${5:?ERROR: Need to define OUTDIST.}
export OUTLIST=${6:?ERROR: Need to define OUTLIST.}

# Other variables
export OUTDIR=data/process/optifit/
export SUBSIZE=10000 # Number of reads to subsample to, based on Baxter, et al., Genome Med, 2016



##################################################
# Clustering Leave-One-Out Outputs Using OptiFit #
##################################################

# Pulling sample name from file name
SAMPLE=$(echo "${INFASTA}" | sed "s:.*/\([0-9]*\).in.fasta:\1:")

# Creating subdirectory using sample name for storing output files
SUBDIR="${OUTDIR}"/"${SAMPLE}"/

# Create subdirectory
mkdir -p "${SUBDIR}"/

# Running OptiFit to cluster left out sample with reference clusters
mothur "#cluster.fit(fasta="${INFASTA}", column="${INDIST}", count="${INCOUNT}", reffasta="${OUTFASTA}", refcolumn="${OUTDIST}", reflist="${OUTLIST}", method=closed, outputdir="${SUBDIR}");
	make.shared(list=current, count=current, label=0.03);
	sub.sample(shared=current, label=0.03, size="${SUBSIZE}")"

# Cleaning up names
for FILE in $(find "${SUBDIR}"/ -type f); do

	# Removing extra name information
	REPLACEMENT=$(echo "${FILE}" | sed "s:in\.::")

	# Rename files using new format
	mv "${FILE}" "${REPLACEMENT}"

done
