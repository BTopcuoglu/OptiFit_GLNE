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
NPROC=$(nproc) # Setting number of processors to use based on available resources



####################################
# Preparing Mothur Reference Files #
####################################

echo PROGRESS: Preparing mothur reference files.

# Making reference output directory
mkdir -p "${TMP}"/

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}")" ]; then
	rm "${OUTDIR}"/*
fi


echo PROGRESS: Preparing SILVA database v4 sequence alignment files.

# Downloading the prepared SILVA database from the mothur website
# This version is from v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# This version is used to maintain methods from Baxter NT, et al. Microbiota-based model improves the sensitivity
# of fecal immunochemical test for detecting colonic lesions. Genome Med. 2016;8(1):37.
wget -N -P "${TMP}"/ http://mothur.org/w/images/1/15/Silva.seed_v123.tgz

# Decompressing the database
tar xvzf "${TMP}"/Silva.seed_v123.tgz -C "${TMP}"/

# Using mothur to pull out the v4 region from bacterial sequences
mothur "#set.current(outputdir="${TMP}"/, processors="${NPROC}");
	get.lineage(fasta="${TMP}"/silva.seed_v123.align, taxonomy="${TMP}"/silva.seed_v123.tax, taxon=Bacteria);
	pcr.seqs(fasta=current, start=11894, end=25319, keepdots=F)"

# Renaming the output file and moving it from the tmp dir to the output dir
mv "${TMP}"/silva.seed_v123.pick.align "${OUTDIR}"/silva.seed.align
mv "${TMP}"/silva.seed_v123.pick.pcr.align "${OUTDIR}"/silva.v4.align



echo PROGRESS: Preparing Ribosomal Database Project taxonomy files.

# Downloading the prepared RDP database from the mothur website
# For more information see http://blog.mothur.org/2015/05/27/RDP-v14-reference_files/
wget -N -P "${TMP}"/ http://mothur.org/w/images/8/88/Trainset14_032015.pds.tgz

# Decompressing the database
tar xvzf "${TMP}"/Trainset14_032015.pds.tgz -C "${TMP}"/

# Move the taxonomy files out of the tmp dir
mv "${TMP}"/trainset14_032015.pds/trainset14_032015* "${OUTDIR}"/

# Cleaning up reference dir
rm -rf "${TMP}"/
