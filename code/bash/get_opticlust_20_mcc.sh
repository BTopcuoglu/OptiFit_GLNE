#! /bin/bash

DIST=${1:?ERROR: Need to define DIST.}
COUNT=${2:?ERROR: Need to define COUNT.}
SPLIT=${3:?ERROR: Need to define SPLIT.} 
NPROC=${4:?ERROR: Need to define NPROC.} 

NUM=`basename $SPLIT .csv`

OUTDIR=data/process/opticlust/$NUM/
SUBSIZE=10000

###################
### GROUP SETUP ###
###################

#initialize train and test arrays
trainIDS=()
testIDS=()

{
    read
    while IFS="," read -a line; do    
        if [[ "${line[1]}" == "train" ]]; then
            trainIDS+="${line[0]}-"
        else
            testIDS+="${line[0]}-"
        fi
    done 
}< "$SPLIT"

#echo "training: ${trainIDS%?}" # %? removes last character (extra "-")
echo "testing: ${testIDS%?}"

##############
### MOTHUR ###
##############

# Make output dirs if they don't exist
mkdir -p "${OUTDIR}"/sub/

### TESTING DATA
# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}"/sub)" ]; then
	rm "${OUTDIR}"/sub/*
fi

# Create OptiClust shared file for training data
 mothur "#set.current(outputdir="${OUTDIR}"/sub/,processors="${NPROC}");
 	get.groups(column="${DIST}", count="${COUNT}", groups=${testIDS%?});
    cluster(column=current, count=current);
    make.shared(list=current, count=current, label=0.03);
	sub.sample(shared=current, label=0.03, size="${SUBSIZE}")"