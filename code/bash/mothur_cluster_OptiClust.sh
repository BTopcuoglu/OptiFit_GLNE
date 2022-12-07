#! /bin/bash
# mothurOptiClust.sh
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
DIST=${4:?ERROR: Need to define DIST.} # Preclustered distance file
NPROC=${5:?ERROR: Need to define NPROC.} # number or processors to use

# Other variables
OUTDIR=data/process/opticlust/shared/ # Output dir since the same shared file will be used for all leave one out steps
#NPROC=$(nproc) # Setting number of processors to use based on available resources
SUBSIZE=10000 # Number of reads to subsample to, based on Baxter, et al., Genome Med, 2016



#############################################
# Create Master Shared File Using OptiClust #
#############################################

# Making output dir if it doesn't already exist
mkdir -p "${OUTDIR}"/

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}")" ]; then
	rm "${OUTDIR}"/*
fi

# Cluster all sequences into master shared file
mothur "#set.current(outputdir="${OUTDIR}"/, processors="${NPROC}", fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}", column="${DIST}");
	cluster(column=current, count=current);
	make.shared(list=current, count=current, label=0.03);
	sub.sample(shared=current, label=0.03, size="${SUBSIZE}")"
