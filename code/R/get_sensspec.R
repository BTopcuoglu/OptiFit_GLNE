#!/usr/bin/env Rscript
# get_sensspec.R
# Courtney Armour
# March 2021
# Pat Schloss Lab
# University of Michigan

#############################################################################
# This script pulls the stats from the sensspec files for each sample to output
# a concatenated fill of values for all samples and optifit/opticlust methods
#############################################################################
library(tidyverse)

#only cluster one time with opticlust so only one sensspec file
opticlust_file <- "data/process/opticlust/shared/glne.precluster.opti_mcc.sensspec"

# get directory paths for all samples from optifit
optifit_sample_dirs <- list.dirs(path="data/process/optifit",recursive=F)

# get file lists for all samples and optifit/opticlust
# IN refers to cluster.fit on the left out sample
# OUT is cluster on the other 489 samples
optifit_in_files <- list.files(path=paste0(optifit_sample_dirs,"/in"),pattern="*.sensspec",full.names=T)
optifit_out_files <- list.files(path=paste0(optifit_sample_dirs,"/out"),pattern="*.sensspec",full.names=T)

#function to read files and add sample and type column
read_optifit_sensspec_files <- function(filenames){
  for(file in filenames){
    # Read the sensspec files and add type(e.g. optifit_in) and sampleID
    data <- read_tsv(file,col_types = cols(.default=col_double())) %>%
      mutate(sample=unlist(strsplit(file,"/"))[4],
             type=paste0(unlist(strsplit(file,"/"))[3],"_",unlist(strsplit(file,"/"))[5]))
  }
  return(data)
}

#generate table for each type
optifit_in.df <- map_df(optifit_in_files, read_optifit_sensspec_files)
optifit_out.df <- map_df(optifit_out_files, read_optifit_sensspec_files)
opticlust.df <- read_tsv(opticlust_file,col_types = cols(.default=col_double())) %>%
  mutate(sample=NA,type="opticlust")

#merge together and write output
all.df <- bind_rows(optifit_in.df,optifit_out.df,opticlust.df)

write_csv(all.df,"results/tables/merged_sensspec.csv")
