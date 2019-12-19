#! /bin/bash
# mothurReferences.sh
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Other variables
export OUTDIR=data/references/ # Directory for storing mothur reference files



####################################
# Preparing Mothur Reference Files #
####################################

echo PROGRESS: Preparing mothur reference files.

# Making reference output directory
mkdir -p "${OUTDIR}"/tmp/


echo PROGRESS: Preparing SILVA database v4 sequence alignment files.

# Downloading the prepared SILVA database from the mothur website
# This version is from v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# This version is used to maintain methods from Baxter NT, et al. Microbiota-based model improves the sensitivity
# of fecal immunochemical test for detecting colonic lesions. Genome Med. 2016;8(1):37.
wget -N -P "${OUTDIR}"/tmp/ http://mothur.org/w/images/1/15/Silva.seed_v123.tgz

# Decompressing the database
tar xvzf "${OUTDIR}"/tmp/Silva.seed_v123.tgz -C "${OUTDIR}"/tmp/

# Using mothur to pull out bacterial sequences and remove sequence gaps
mothur "#get.lineage(fasta="${OUTDIR}"/tmp/silva.seed_v123.align, taxonomy="${OUTDIR}"/tmp/silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta="${OUTDIR}"/tmp/silva.seed_v123.pick.align, processors=8)"

# Renaming the output file and moving it from the tmp dir to the output dir
mv "${OUTDIR}"/tmp/silva.seed_v123.pick.align "${OUTDIR}"/silva.seed.align

# Using mothur to only keep sequences from the v4 region of the 16S rRNA DNA region
mothur "#pcr.seqs(fasta="${OUTDIR}"/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"

# Renaming the final v4 SILVA reference file
mv "${OUTDIR}"/silva.seed.pcr.align "${OUTDIR}"/silva.v4.align



echo PROGRESS: Preparing Ribosomal Database Project taxonomy files.

# Downloading the prepared RDP database from the mothur website
# For more information see http://blog.mothur.org/2015/05/27/RDP-v14-reference_files/
wget -N -P "${OUTDIR}"/tmp/ http://mothur.org/w/images/8/88/Trainset14_032015.pds.tgz

# Decompressing the database
tar xvzf "${OUTDIR}"/tmp/Trainset14_032015.pds.tgz -C "${OUTDIR}"/tmp/

# Move the taxonomy files out of the tmp dir
mv "${OUTDIR}"/tmp/trainset14_032015.pds/trainset14_032015* "${OUTDIR}"/

# Cleaning up reference dir
rm -rf "${OUTDIR}"/tmp/
