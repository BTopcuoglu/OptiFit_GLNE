#! /bin/bash
# mothurOpitFit.sh
# Courtney R. Armour
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
SPLIT=${5:?ERROR: Need to define SPLIT.} #file defining test/train split
NUM=`basename $SPLIT .csv`
NPROC=${6:?ERROR: Need to define NPROC.} # number or processors to use
OUTDIR=${7:?ERROR: Need to define OUTDIR}
SUBSIZE=${8:?ERROR: Need to define SUBSIZE.}

OUTDIR=${OUTDIR}/${NUM}/train/

###################
### GROUP SETUP ###
###################

#initialize train arrays
trainIDS=()

#assign IDS to train array
{
    read
    while IFS="," read -a line; do    
        if [[ "${line[1]}" == "train" ]]; then
            trainIDS+="${line[0]}-"
        fi
    done 
}< "$SPLIT"

echo "training: ${trainIDS%?}" # %? removes last character (extra "-")

############################################
# Generate Shared After Leaving Sample Out #
############################################

# Make output dir if it doesn't exist
mkdir -p "${OUTDIR}"/

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}")" ]; then
	rm "${OUTDIR}"/*
fi

# Cluster all sequences for the training set
mothur "#set.current(outputdir="${OUTDIR}"/, processors="${NPROC}");
	get.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups=${trainIDS%?}, column="${DIST}");
	cluster(column=current, count=current);
	make.shared(list=current, count=current, label=0.03);
	sub.sample(shared=current, label=0.03, size="${SUBSIZE}")"
    
#rename file for consistency 
mv data/process/optifit_self/${NUM}/train/glne.precluster.pick.opti_mcc.0.03.subsample.shared data/process/optifit_self/${NUM}/train/glne.optifit_self.shared
