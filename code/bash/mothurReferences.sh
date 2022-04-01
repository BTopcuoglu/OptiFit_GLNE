#! /bin/bash
# mothurReferences.sh
# William L. Close
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Other variables
OUTDIR=data/references/ # Directory for storing mothur reference files
TMP="${OUTDIR}"/tmp/
#NPROC=$(nproc) # Setting number of processors to use based on available resources
NPROC=${1:?ERROR: Need to define NPROC.}


####################################
# Preparing Mothur Reference Files #
####################################

echo PROGRESS: Preparing mothur reference files.

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}")" ]; then
	rm -rf "${OUTDIR}"/*
fi

# Making reference output directory
mkdir -p "${TMP}"/



echo PROGRESS: Preparing SILVA database v4 sequence alignment files.

# Downloading the prepared SILVA database from the mothur website
wget -N -P "${TMP}"/ https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v132.tgz

# Decompressing the database
tar xvzf "${TMP}"/silva.seed_v132.tgz -C "${TMP}"/

# Using mothur to pull out the v4 region from bacterial sequences
mothur "#set.current(outputdir="${TMP}"/, processors="${NPROC}");
	get.lineage(fasta="${TMP}"/silva.seed_v132.align, taxonomy="${TMP}"/silva.seed_v132.tax, taxon=Bacteria);
	pcr.seqs(fasta=current, start=11894, end=25319, keepdots=F)"

# Renaming the output file and moving it from the tmp dir to the output dir
mv "${TMP}"/silva.seed_v132.pick.pcr.align "${OUTDIR}"/silva.v4.align

echo PROGRESS: Preparing Ribosomal Database Project taxonomy files.

# Downloading the prepared RDP database from the mothur website
wget -N -P "${TMP}"/ https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.pds.tgz

# Decompressing the database
tar xvzf "${TMP}"/trainset16_022016.pds.tgz -C "${TMP}"/

# Move the taxonomy files out of the tmp dir
mv "${TMP}"/trainset16_022016.pds/trainset16_022016.pds.* "${OUTDIR}"/

# Cleaning up reference dir
rm -rf "${TMP}"/
