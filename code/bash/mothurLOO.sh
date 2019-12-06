#! /bin/bash
# mothurLOO.sh
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export GROUPS=${1:?ERROR: Need to define GROUPS.}
export FASTA=${2:?ERROR: Need to define FASTA.}
export COUNT=${3:?ERROR: Need to define COUNT.}
export TAXONOMY=${4:?ERROR: Need to define TAXONOMY.}
export NUM=${5:?ERROR: Need to define NUM.}

# Other variables
export WORKDIR=data/process/baxter/final
export OUTDIR=data/process/

########################################################
# Generate shared file for only one sample  #
########################################################


# Now let's extract fasta, taxonomy and count for the removed group and build subsampled shared for the removed sample.

mothur "#get.groups(groups="${GROUPS}", fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${NUM}");
dist.seqs(fasta=current, cutoff=0.03);
cluster(column=current, count=current);
make.shared(list=current, count=current, label=0.03);
sub.sample(shared=current, label=0.03)"


# Change the name of the files generated to represent that they only have 1 SAMPLE

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.dist data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick."${NUM}".dist

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.list data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.sample.list

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.steps data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.sample.steps

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.sensspec data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.sample.sensspec

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.shared data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.sample.shared

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.0.03.subsample.shared data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.0.03.subsample.sample.shared

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.pick.taxonomy data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.pick.sample.taxonomy

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.pick.count_table data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.pick.sample.count_table

mv data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.fasta data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.sample.fasta



########################################################
# Generate shared file for all samples but the one #
########################################################

# The generated subsampled file will have all the samples except the left-out-one.

# Now let's remove that 1 sample from rest of the samples by removing from original count, fasta and taxa files.

mothur "#remove.groups(groups="${GROUPS}", fasta="${FASTA}", count="${COUNT}", taxonomy="${TAXONOMY}",  groups="${NUM}");
dist.seqs(fasta=current, cutoff=0.03);
cluster(column=current, count=current);
make.shared(list=current, count=current, label=0.03);
sub.sample(shared=current, label=0.03)"










###############
# Cleaning Up #
###############

echo PROGRESS: Cleaning up working directory.

# Making dir for storing intermediate files (can be deleted later)
mkdir -p "${OUTDIR}"/intermediate/

# Deleting unneccessary files
rm "${OUTDIR}"/*filter.unique.precluster*fasta
rm "${OUTDIR}"/*filter.unique.precluster*map
rm "${OUTDIR}"/*filter.unique.precluster*count_table

# Moving all remaining intermediate files to the intermediate dir
mv "${OUTDIR}"/stability* "${OUTDIR}"/intermediate/
