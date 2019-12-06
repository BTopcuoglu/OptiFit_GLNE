#! /bin/bash
# mothurContigs.sh
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
#export SAMPLEDIR=${1:?ERROR: Need to define SAMPLEDIR.}
#export FILE=${2:?ERROR: Need to define SAMPLEDIR.}
#export SILVAV4=${3:?ERROR: Need to define SILVAV4.}
#export RDPFASTA=${4:?ERROR: Need to define RDPFASTA.}
#export RDPTAX=${5:?ERROR: Need to define RDPTAX.}

# Other variables
export OUTDIR=data/process/baxter/intermediate
export FINAL= data/process/baxter/final



###################
# Run QC Analysis #
###################

echo PROGRESS: Creating contigs for all the samples

# Making output dir
mkdir -p "${OUTDIR}"

# Convert to fasta files that will be used
for sample in data/mothur/raw/baxter/*.sra
do
	fastq-dump --split-files $sample -O data/process/baxter
done

# Making contigs from fastq files, aligning reads to references, removing any non-bacterial sequences, calculating distance matrix, making shared file.

mothur "#make.contigs(file=data/process/baxter/glne007.files, outputdir="${OUTDIR}");
	screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference="${SILVAV4}");
	screen.seqs(fasta=current, count=current, start=13862, end=23444, maxhomop=8);
	filter.seqs(fasta=current, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=T);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference="${RDPFASTA}", taxonomy="${RDPTAX}", cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
	remove.groups(fasta=current, count=current, taxonomy=current, groups=mock1-mock2-mock5-mock6-mock7)"

###############
# Cleaning Up #
###############

echo PROGRESS: Cleaning up working directory.

# Moving important output files to the final directory for future use

cp "${OUTDIR}"/glne007.contigs.good.groups "${FINAL}"/full.groups

cp "${OUTDIR}"/data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta "${FINAL}"/full.fasta

cp "${OUTDIR}"/data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table "${FINAL}"/full.count_table

cp "${OUTDIR}"/data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy "${FINAL}"/full.taxonomy
