#! /bin/bash
# mothurLOO.sh
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
#export GROUPS=${1:?ERROR: Need to define GROUPS.}
#export FASTA=${2:?ERROR: Need to define FASTA.}
#export COUNT=${3:?ERROR: Need to define COUNT.}
#export TAXONOMY=${4:?ERROR: Need to define TAXONOMY.}
#export NUM=${5:?ERROR: Need to define NUM.}

# Other variables
export WORKDIR=data/process/baxter/final
export OUTDIR=data/process/
export NUM=2003650
export GROUP=data/process/baxter/final/full.groups
export FASTA=data/process/baxter/final/full.fasta
export COUNT=data/process/baxter/final/full.count_table
export TAXONOMY=data/process/baxter/final/full.taxonomy
########################################################
# Generate shared file for only one sample  #
########################################################


# Now let's extract fasta, taxonomy and count for the removed group and build subsampled shared for the removed sample.

mothur "#get.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${NUM}");
dist.seqs(fasta=current, cutoff=0.03)"


# Change the name of the files generated to represent that they only have 1 SAMPLE


mv data/process/baxter/final/full.pick.fasta "${OUTDIR}"/sample."${NUM}".fasta

mv data/process/baxter/final/full.pick.count_table "${OUTDIR}"/sample."${NUM}".count_table

mv data/process/baxter/final/full.pick.taxonomy "${OUTDIR}"/sample."${NUM}".taxonomy

mv data/process/baxter/final/full.pick.dist "${OUTDIR}"/sample."${NUM}".dist

########################################################
# Generate shared file for all samples but the one #
########################################################

# The generated subsampled file will have all the samples except the left-out-one.

# Now let's remove that 1 sample from rest of the samples by removing from original count, fasta and taxa files.

mothur "#remove.groups(fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${NUM}");
dist.seqs(fasta=current, cutoff=0.03);
cluster(column=current, count=current);
make.shared(list=current, count=current, label=0.03);
sub.sample(shared=current, label=0.03)"

mv data/process/baxter/final/full.pick.fasta "${OUTDIR}"/without."${NUM}".fasta

mv data/process/baxter/final/full.pick.count_table "${OUTDIR}"/without."${NUM}".count_table

mv data/process/baxter/final/full.pick.taxonomy "${OUTDIR}"/without."${NUM}".taxonomy

mv data/process/baxter/final/full.pick.dist "${OUTDIR}"/without."${NUM}".dist

mv data/process/baxter/final/full.pick.opti_mcc.list "${OUTDIR}"/without.opti_mcc."${NUM}".list

mv data/process/baxter/final/full.pick.opti_mcc.steps "${OUTDIR}"/without.opti_mcc."${NUM}".steps

mv data/process/baxter/final/full.pick.opti_mcc.sensspec "${OUTDIR}"/without.opti_mcc."${NUM}".sensspec

mv data/process/baxter/final/full.pick.opti_mcc.shared "${OUTDIR}"/without.opti_mcc."${NUM}".shared

mv data/process/baxter/final/full.pick.opti_mcc.0.03.subsample.shared "${OUTDIR}"/without.opti_mcc.0.03."${NUM}".subsample.shared
