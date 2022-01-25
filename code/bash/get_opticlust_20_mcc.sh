#! /bin/bash

LIST=${1:?ERROR: Need to define LIST.}
DIST=${2:?ERROR: Need to define DIST.}
COUNT=${3:?ERROR: Need to define COUNT.}
SPLIT=${4:?ERROR: Need to define SPLIT.} 
NPROC=${5:?ERROR: Need to define NPROC.} 

NUM=`basename $SPLIT .csv`

OUTDIR=data/process/opticlust/$NUM/

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

echo "training: ${trainIDS%?}" # %? removes last character (extra "-")
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
 	get.groups(list="${LIST}", column="${DIST}", count="${COUNT}", groups=${testIDS%?});
    sens.spec(list=current,column=current,count=current)"