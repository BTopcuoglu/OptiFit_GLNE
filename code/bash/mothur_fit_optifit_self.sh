#! /bin/bash
# mothurFitOptiFit.sh
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
SPLIT=${5:?ERROR: Need to define SPLIT.} # Split of test and training sets
REFFASTA=${6:?ERROR: Need to define TRAINFASTA.} # Reference fasta file with test data removed
REFDIST=${7:?ERROR: Need to define TRAINDIST.} # Reference dist file with test data removed
REFLIST=${8:?ERROR: Need to define TRAINLIST.} # Reference list file with test data removed
NPROC=${9:?ERROR: Need to define NPROC.} # number or processors to use
OUTDIR=${10:?ERROR: Need to define OUTDIR.} # Output dir based on split number 
SUBSIZE=${11:?ERROR: Need to define SUBSIZE.} # Number of reads to subsample to, based on Baxter, et al., Genome Med, 2016

# Other variables
NUM=`basename $SPLIT .csv`

OUTDIR=${OUTDIR}/${NUM}/test/
echo $OUTDIR

###################
### GROUP SETUP ###
###################

#initialize test arrays
testIDS=()

#assign IDS to test array
{
    read
    while IFS="," read -a line; do    
        if [[ "${line[1]}" == "test" ]]; then
            testIDS+="${line[0]}-"
        fi
    done 
}< "$SPLIT"

echo "testing: ${testIDS%?}" # %? removes last character (extra "-")

##########################################
# Clustering the test data Using OptiFit #
##########################################

# Make output dir if it doesn't exist
mkdir -p "${OUTDIR}"/

# Removing old files if they exist
if [ -n "$(ls -A "${OUTDIR}")" ]; then
	rm "${OUTDIR}"/*
fi

# Selecting test data to be fit, subsampling input files to account for 'closed' clustering method, and fitting sample to previous shared file
mothur "#set.current(outputdir="${OUTDIR}"/, processors="${NPROC}");
	get.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}", groups=${testIDS%?});
	rename.seqs(fasta=current,count=current);
	cluster.fit(fasta=current, count=current, reffasta="${REFFASTA}", refcolumn="${REFDIST}", reflist="${REFLIST}", method=open, printref=f);
	make.shared(list=current, count=current, label=0.03);
	sub.sample(shared=current, size="${SUBSIZE}")"
    
mv data/process/optifit_self/${NUM}/test/glne.precluster.pick.renamed.fit.optifit_mcc.0.03.subsample.shared data/process/optifit_self/${NUM}/test/glne.optifit_self.shared
	