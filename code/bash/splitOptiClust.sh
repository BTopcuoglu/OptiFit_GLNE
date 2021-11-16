#! /bin/bash
# splitOptiClust.sh
# Courtney Armour
# Schloss Lab
# University of Michigan

# Purpose: use split file to create train/test shared files

#argvals=("$@")

SHARED=${1:?ERROR: Need to define SHARED.}
SPLIT=${2:?ERROR: Need to define SPLITS.} 
NUM=`basename $SPLIT .csv`

#SHARED=${argvals[0]}
#SPLITS=(${argvals[@]:1}) 

OUTDIR=data/process/opticlust/$NUM/ # Output dir based on sample name to keep things separate during parallelization/organized
NPROC=$(nproc) # Setting number of processors to use based on available resources

###################
### GROUP SETUP ###
###################

#initialize train and test arrays
trainIDS=()
testIDS=()

#assign IDS to train or test array

#for f in ${SPLITS[@]}
#do
    #echo "processing $f"
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
#done
echo "training: ${trainIDS%?}" # %? removes last character (extra "-")
echo "testing: ${testIDS%?}"

##############
### MOTHUR ###
##############

# Make output dirs if they don't exist
mkdir -p "${OUTDIR}"/ "${OUTDIR}"/train/ "${OUTDIR}"/test/

### TRAINING DATA
# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}"/train)" ]; then
	rm "${OUTDIR}"/train/*
fi

# # Create OptiClust shared file for training data
 mothur "#set.current(outputdir="${OUTDIR}"/train/, processors="${NPROC}");
 	get.groups(shared="${SHARED}", groups=${trainIDS%?})"
     
### TESTING DATA
# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}"/test)" ]; then
	rm "${OUTDIR}"/test/*
fi

# Create OptiClust shared file for training data
 mothur "#set.current(outputdir="${OUTDIR}"/test/, processors="${NPROC}");
 	get.groups(shared="${SHARED}", groups=${testIDS%?})"