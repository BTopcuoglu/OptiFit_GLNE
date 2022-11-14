#! /bin/bash
# Courtney Armour
# Schloss Lab
# University of Michigan

OUTDIR=data/references/ # Directory for storing mothur reference files
TMP="${OUTDIR}"/tmp/
NPROC=${1:?ERROR: Need to define NPROC.}

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
# remove readmes
rm "${TMP}"/README*

# Using mothur to pull out the v4 region from bacterial sequences
mothur "#set.current(outputdir="${TMP}"/, processors="${NPROC}");
	get.lineage(fasta="${TMP}"/silva.seed_v132.align, taxonomy="${TMP}"/silva.seed_v132.tax, taxon=Bacteria);
	pcr.seqs(fasta=current, start=11894, end=25319, keepdots=F)"

# Renaming the output file and moving it from the tmp dir to the output dir
mv "${TMP}"/silva.seed_v132.pick.pcr.align "${OUTDIR}"/silva.v4.align
mv "${TMP}"/silva.seed_v132.tax "${OUTDIR}"/silva.v4.tax

echo PROGRESS: Preparing Ribosomal Database Project taxonomy files.

# Downloading the prepared RDP database from the mothur website
wget -N -P "${TMP}"/ https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.pds.tgz

# Decompressing the database
tar xvzf "${TMP}"/trainset16_022016.pds.tgz -C "${TMP}"/

# remove readmes
rm "${TMP}"/README*

# Move the taxonomy files out of the tmp dir
mv "${TMP}"/trainset16_022016.pds/trainset16_022016.pds.* "${OUTDIR}"/

echo PROGRESS: Preparing Greengenes Database files.

#download the 
wget -N -P "${TMP}"/ https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_13_8_99.taxonomy.tgz
wget -N -P "${TMP}"/ https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_13_8_99.refalign.tgz

tar -xzvPf "${TMP}"/gg_13_8_99.taxonomy.tgz -C "${TMP}"/
rm "${TMP}"/README*
tar -xzvPf "${TMP}"/gg_13_8_99.refalign.tgz -C "${TMP}"/
rm "${TMP}"/README*

Rscript code/R/format_gg_tax.R 

# filter to just bacteria, align to SILVA, and select v4 region
mothur "#set.current(outputdir="${TMP}"/, processors="${NPROC}");
	get.lineage(fasta="${TMP}"/gg_13_8_99.fasta, taxonomy="${TMP}"/gg.tax, taxon=Bacteria);
	align.seqs(candidate=current, template="${TMP}"/silva.seed_v132.align);
	screen.seqs(fasta=current, maxambig=0, maxhomop=8, start=11894, end=25319);	
    pcr.seqs(fasta=current, start=11894, end=25319, keepdots=F)"
	
mv "${TMP}"/gg_13_8_99.pick.good.pcr.align "${OUTDIR}"/gg.fasta
mv "${TMP}"/gg.tax "${OUTDIR}"/gg.tax

# Cleaning up reference dir
rm -rf "${TMP}"/
rm "${OUTDIR}"/README*