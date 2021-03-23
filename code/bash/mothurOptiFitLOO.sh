#! /bin/bash
# mothurLOO.sh
# William L. Close
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
FASTA=${1:?ERROR: Need to define FASTA.} # Preclustered fasta file
COUNT=${2:?ERROR: Need to define COUNT.} # Preclustered count file
TAXONOMY=${3:?ERROR: Need to define TAXONOMY.} # Preclustered tax file
DIST=${4::?ERROR: Need to define DIST.} # Preclustered dist file
SAMPLE=${5:?ERROR: Need to define SAMPLE.} # Sample to be removed

# Other variables
OUTDIR=data/process/optifit/"${SAMPLE}"/out/ # Output dir based on sample name to keep things separate during parallelization/organized
NPROC=$(nproc) # Setting number of processors to use based on available resources
SUBSIZE=10000 # Number of reads to subsample to, based on Baxter, et al., Genome Med, 2016

############################################
# Generate Shared After Leaving Sample Out #
############################################

# Make output dir if it doesn't exist
mkdir -p "${OUTDIR}"/

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}")" ]; then
	rm "${OUTDIR}"/*
fi

# Cluster all sequences after leaving out the specified sample
mothur "#set.current(outputdir="${OUTDIR}"/, processors="${NPROC}");
	remove.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${SAMPLE}", column="${DIST}");
	cluster(column=current, count=current);
	make.shared(list=current, count=current, label=0.03);
	sub.sample(shared=current, label=0.03, size="${SUBSIZE}")"
