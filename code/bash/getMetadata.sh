#! /bin/bash
# getMetadata.sh
# Begum Topcuoglu
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Other variables
export OUTDIR=data/metadata/



###########################
# Download Metadata Table #
###########################

# Making output dir
mkdir -p "${OUTDIR}"

# Downloading metadata from https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015
curl -o "${OUTDIR}"/metadata.tsv https://raw.githubusercontent.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015/master/data/metadata.tsv
