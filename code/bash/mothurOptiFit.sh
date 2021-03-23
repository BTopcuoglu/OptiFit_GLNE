#! /bin/bash
# mothurOptiFit.sh
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
SAMPLE=${5:?ERROR: Need to define SAMPLE.} # Sample to be fit to reference files
REFFASTA=${6:?ERROR: Need to define OUTFASTA.} # Reference fasta file with sample removed
REFDIST=${7:?ERROR: Need to define OUTDIST.} # Reference dist file with sample removed
REFLIST=${8:?ERROR: Need to define OUTLIST.} # Reference list file with sample removed


# Other variables
OUTDIR=data/process/optifit/"${SAMPLE}"/in/ # Output dir based on sample name to keep things separate during parallelization/organized
NPROC=$(nproc) # Setting number of processors to use based on available resources
SUBSIZE=10000 # Number of reads to subsample to, based on Baxter, et al., Genome Med, 2016



##################################################
# Clustering Leave-One-Out Outputs Using OptiFit #
##################################################

# Make output dir if it doesn't exist
mkdir -p "${OUTDIR}"/

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}")" ]; then
	rm "${OUTDIR}"/*
fi

# Selecting sample to be fit, subsampling input files to account for 'closed' clustering method, and fitting sample to previous shared file
mothur "#set.current(outputdir="${OUTDIR}"/, processors="${NPROC}");
	get.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${SAMPLE}", column="${DIST}");
	sub.sample(fasta=current, count=current, taxonomy=current, size="${SUBSIZE}");
	rename.seqs(fasta=current,count=current);
	cluster.fit(fasta=current, column=current, count=current, reffasta="${REFFASTA}", refcolumn="${REFDIST}", reflist="${REFLIST}", method=closed);
	make.shared(list=current, count=current, label=0.03);
	list.seqs(list=current)"  # creates accnos file of seqs in OTUs
