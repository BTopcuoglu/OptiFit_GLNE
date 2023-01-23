#! /bin/bash
# COURTNEY ARMOUR
# OCTOBER 2022
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
dbFASTA=${1:?ERROR: Need to define dbFASTAV4.}
dbTAX=${2:?ERROR: Need to define dbTAX.}
NPROC=${3:?ERROR: Need to define NPROC.}

# Other variables
OUTDIR=data/process/precluster/
TMP="${OUTDIR}"/tmp/

# Making reference output directory
mkdir -p "${TMP}"/

#run mothur to precluster silva reference sequences
mothur "#set.current(outputdir="${TMP}"/, processors="${NPROC}");
        filter.seqs(fasta="${dbFASTA}"); 
        unique.seqs(fasta=current, format=count);
        pre.cluster(fasta=current, count=current, diffs=2);
        dist.seqs(fasta=current);
        list.seqs(fasta=current);
        get.seqs(accnos=current, taxonomy="${dbTAX}");
        get.seqs(accnos=current, count=current);
        summary.seqs(fasta=current, count=current);"

# move and rename files
mv "${TMP}"/gg.filter.unique.precluster.fasta "${OUTDIR}"/gg.precluster.fasta
mv "${TMP}"/gg.filter.unique.precluster.pick.count_table "${OUTDIR}"/gg.precluster.count_table
mv "${TMP}"/gg.filter.unique.precluster.dist "${OUTDIR}"/gg.precluster.dist
mv "${TMP}"/gg.pick.tax "${OUTDIR}"/gg.precluster.tax

# Cleaning up reference dir
# rm -rf "${TMP}"/