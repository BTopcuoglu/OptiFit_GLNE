#! /bin/bash
# mothurContigs.sh
# Begum Topcuoglu
# Schloss Lab
# University of Michigan


export WORKDIR=data/process
export NUM=2003650

# We are now ready to use OptiFit to fit left out sample to the rest.

mothur "#cluster.fit(fasta="${WORKDIR}"/sample."${NUM}".fasta, column="${WORKDIR}"/sample."${NUM}".dist, count="${WORKDIR}"/sample."${NUM}".count_table, reffasta="${WORKDIR}"/without."${NUM}".fasta, refcolumn="${WORKDIR}"/without."${NUM}".dist, reflist="${WORKDIR}"/without.opti_mcc."${NUM}".list, method=closed)"
