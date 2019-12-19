#! /bin/bash
# mothurContigs.sh
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export FILESFILE=${1:?ERROR: Need to define FILESFILE.} # File listing sample groups and sequence files
export SILVAV4=${2:?ERROR: Need to define SILVAV4.}
export RDPFASTA=${3:?ERROR: Need to define RDPFASTA.}
export RDPTAX=${4:?ERROR: Need to define RDPTAX.}

# Other variables
export OUTDIR=data/process/precluster/
export TMP="${OUTDIR}"/tmp/


###################
# Run QC Analysis #
###################

echo PROGRESS: Creating contigs for all the samples

# Making output dir
mkdir -p "${TMP}"/

mothur "#make.contigs(file="${FILESFILE}", outputdir="${TMP}"/);
	screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference="${SILVAV4}");
	screen.seqs(fasta=current, count=current, start=1968, end=11550);
	filter.seqs(fasta=current, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=T);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference="${RDPFASTA}", taxonomy="${RDPTAX}", cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"


# Moving and renaming important files for downstream use
mv "${TMP}"/*.contigs.good.groups "${OUTDIR}"/glne.precluster.groups
mv "${TMP}"/*.precluster.pick.pick.fasta "${OUTDIR}"/glne.precluster.fasta
mv "${TMP}"/*.precluster.denovo.vsearch.pick.pick.count_table "${OUTDIR}"/glne.precluster.count_table
mv "${TMP}"/*.precluster.pick.pds.wang.pick.taxonomy "${OUTDIR}"/glne.precluster.taxonomy


###############
# Cleaning Up #
###############

echo PROGRESS: Cleaning up working directory.

# Deleting unneccessary files
rm -r "${TMP}"/
