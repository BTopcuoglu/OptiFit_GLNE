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
library(assertthat)

outDir <- "data/learning/summary/"
if(!dir.exists(outDir)){dir.create(outDir)}

# grab input files
input <- commandArgs(trailingOnly=TRUE)
# input <- c("data/process/opticlust/shared/glne.precluster.opti_mcc.sensspec",
#            "data/process/optifit/split_1/train/glne.precluster.pick.opti_mcc.sensspec",
#            "data/process/optifit/split_2/train/glne.precluster.pick.opti_mcc.sensspec",
#            "data/process/optifit/split_1/test/glne.precluster.pick.subsample.renamed.fit.optifit_mcc.steps",
#            "data/process/optifit/split_2/test/glne.precluster.pick.subsample.renamed.fit.optifit_mcc.steps")

#function to read in mcc files
read_mcc <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    algorithm <- file_parts[grepl(pattern="^opticlust$|^optifit$",file_parts)]
    if(algorithm == "opticlust"){

        #get just the mcc score and add the other info
        mcc <- read_tsv(file,col_types = cols(.default = col_double()))  %>% 
            select(mcc)  %>% 
            mutate(split=NA,
                   algorithm=algorithm,
                   subset=NA,
                   state=NA)
    }else{
        split <- file_parts[grepl(pattern="split",file_parts)]
        subset <- file_parts[grepl(pattern="train|test",file_parts)]
        
        if(subset == "test"){
            #select the last iteration of the steps file
            steps <- read_tsv(file)  %>% 
                select(state,iter,mcc) %>% 
                filter( iter == max(iter)) 
            #make sure we have both states
            assert_that(nrow(steps) == 2)
            #get mcc score
            mcc <- steps  %>% 
                select(state,mcc)  %>% 
                mutate(split=split,
                       algorithm=algorithm,
                       subset=subset)  %>% 
                select(mcc,split,algorithm,subset,state)
            
        }else{
            #get just the mcc score and add the other info
            mcc <- read_tsv(file,col_types = cols(.default = col_double()))  %>% 
                select(mcc)  %>% 
                mutate(split=split,
                       algorithm=algorithm,
                       subset=subset,
                       state=NA)
        }
    }
    
    return(mcc)
}

# map files to datafram
mergedMCC <- map_dfr(input,read_mcc)
write_csv(mergedMCC,paste0(outDir,"merged_MCC.csv"))