#! /bin/bash
# mothurOptiClustLOO.sh
# William L. Close
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
SHARED=${1:?ERROR: Need to define SHARED.} # Master OptiClust shared file generated using all samples
SAMPLE=${2:?ERROR: Need to define SAMPLE.} # Sample to be left out

# Other variables
OUTDIR=data/process/opticlust/"${SAMPLE}"/ # Output dir based on sample name to keep things separate during parallelization/organized
NPROC=$(nproc) # Setting number of processors to use based on available resources



############################
# Removing Specific Sample #
############################

# Make output dirs if they don't exist
mkdir -p "${OUTDIR}"/ "${OUTDIR}"/in/ "${OUTDIR}"/out/

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}"/out/)" ]; then
	rm "${OUTDIR}"/*
fi

# Create OptiClust shared file without specified sample
mothur "#set.current(outputdir="${OUTDIR}"/out/, processors="${NPROC}");
	remove.groups(shared="${SHARED}", groups="${SAMPLE}")"



########################################
# Removing All Except Specified Sample #
########################################

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}"/in/)" ]; then
	rm "${OUTDIR}"/*
fi

# Create OptiClust shared file with only the specified sample
mothur "#set.current(outputdir="${OUTDIR}"/in/, processors="${NPROC}");
	get.groups(shared="${SHARED}", groups="${SAMPLE}")"
