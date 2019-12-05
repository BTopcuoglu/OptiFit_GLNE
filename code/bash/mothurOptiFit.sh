#! /bin/bash
# mothurContigs.sh
# Begum Topcuoglu
# Schloss Lab
# University of Michigan



# We are now ready to use OptiFit to fit left out sample to the rest.

mothur cluster.fit(fasta=data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.sample.fasta, column=data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.sample.dist, count=data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.pick.sample.count_table, reffasta=data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.fasta, refcolumn=data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.dist, reflist=data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.pick.opti_mcc.list, method=closed)
