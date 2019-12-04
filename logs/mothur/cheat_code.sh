cluster(column=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table


classify.otu(list=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)

data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta


get.oturep(fasta=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, label=0.03, method=abundance)


get.groups(shared=data/mothur/process//final.shared, groups=Neg1-Neg2-Neg3-Neg4-Neg5-Neg6-Neg7-Neg8)


remove.groups(shared=data/mothur/process//final.shared, groups=Water5a-Water5b-Water6a-Water6b-Water7a-Water7b-bufferP1.1-bufferP1.2-bufferP2.1-bufferP2.1b-bufferP2.2-bufferP2.2b-bufferP2.3-bufferP2.4-bufferP2.5-Neg1-Neg2-Neg3-Neg4-Neg5-Neg6-Neg7-Neg8-Mock1-Mock2-Mock3-Mock4-Mock5-Mock6-Mock7-Mock8)
