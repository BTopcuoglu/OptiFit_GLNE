#! /bin/bash
# mothurOptiFit.sh
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export INFASTA=${1:?ERROR: Need to define INFASTA.}
export INDIST=${2:?ERROR: Need to define INDIST.}
export INCOUNT=${3:?ERROR: Need to define INCOUNT.}
export OUTFASTA=${4:?ERROR: Need to define OUTFASTA.}
export OUTDIST=${5:?ERROR: Need to define OUTDIST.}
export OUTLIST=${6:?ERROR: Need to define OUTLIST.}

# Other variables
export OUTDIR=data/process/optifit/



##################################################
# Clustering Leave-One-Out Outputs Using OptiFit #
##################################################

# Running OptiFit to cluster left out sample with reference clusters
mothur "#cluster.fit(fasta="${INFASTA}", column="${INDIST}", count="${INCOUNT}", reffasta="${OUTFASTA}", refcolumn="${OUTDIST}", reflist="${OUTLIST}", method=closed, outputdir="${OUTDIR}");
	make.shared(list=current, count=current, label=0.03)"
