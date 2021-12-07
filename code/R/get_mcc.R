# get_mcc.R
# Courtney Armour
# November 2021
# Pat Schloss Lab
# University of Michigan

#############################################################################
# This script pulls the stats from the mothur sensspec files for each sample 
# to output a concatenated file of values for all samples and methods
#############################################################################
library(tidyverse)


outDir <- "results/tables/"

# grab input files
input <- commandArgs(trailingOnly=TRUE)
# input <- c("data/process/opticlust/shared/glne.precluster.opti_mcc.sensspec",
#            "data/process/optifit/split_1/train/glne.precluster.pick.opti_mcc.sensspec",
#            "data/process/optifit/split_2/train/glne.precluster.pick.opti_mcc.sensspec",
#            "data/process/optifit/split_1/test/glne.precluster.pick.subsample.renamed.fit.optifit_mcc.sensspec",
#            "data/process/optifit/split_2/test/glne.precluster.pick.subsample.renamed.fit.optifit_mcc.sensspec")

#function to read in mcc files
read_mcc <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    algorithm <- file_parts[grepl(pattern="^opticlust$|^optifit$",file_parts)]
    if(algorithm == "opticlust"){
        split <- NA
        subset <- NA
    }else{
        split <- file_parts[grepl(pattern="split",file_parts)]
        subset <- file_parts[grepl(pattern="train|test",file_parts)]
    }
    
    #get just the mcc score and add the other info
    mcc <- read_tsv(file,col_types = cols(.default = col_double()))  %>% 
        select(mcc)  %>% 
        mutate(split=split,
               algorithm=algorithm,
               subset=subset)
    return(mcc)
}

# map files to datafram
mergedMCC <- map_dfr(input,read_mcc)
write_csv(mergedMCC,paste0(outDir,"mergedMCC.csv"))