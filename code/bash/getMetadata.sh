#! /bin/bash
# getMetadata.sh
# William L. Close
# Begum Topcuoglu
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Other variables
OUTDIR=data/metadata/



###########################
# Download Metadata Table #
###########################

# Making output dir
mkdir -p "${OUTDIR}"

# Downloading metadata from https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015
curl -o "${OUTDIR}"/metadata.tsv https://raw.githubusercontent.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015/master/data/metadata.tsv
