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

# grab input files
input <- unlist(snakemake@input)

# #for testing:
# #opticlust mcc
# in1 <- "data/process/opticlust_denovo/glne.precluster.opti_mcc.sensspec"
# #optifit reference mcc per split
# in2 <- paste0("data/process/optifit_self/split_",seq(1,100),"/train/glne.precluster.pick.opti_mcc.sensspec")
# #optifit fit mcc per split
# in3 <- paste0("data/process/optifit_self/split_",seq(1,100),"/test/glne.precluster.pick.renamed.fit.optifit_mcc.steps")
# #greengenes mcc
# in4 <- "data/process/optifit_gg/glne.precluster.fit.optifit_mcc.sensspec"
# #vsearch gg mcc
# in5 <- "data/process/vsearch_gg/glne.vsearch_gg.userLabel.pick.sensspec"
# #vsearch denovo mcc
# in6 <- "data/process/vsearch_denovo/glne.vsearch_denovo.userLabel.pick.sensspec"
# input <- c(in1,in2,in3,in4,in5,in6)

#function to read in mcc files
read_mcc <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    algorithm <- file_parts[grepl(pattern="^opticlust|^optifit|^vsearch",file_parts)][1]
    if(algorithm == "optifit_self"){
        split <- file_parts[grepl(pattern="split",file_parts)]
        subset <- file_parts[grepl(pattern="train|test",file_parts)]
        
        if(subset == "test"){
            #select the last iteration of the steps file
            steps <- read_tsv(file)  %>% 
                select(state,iter,mcc) %>% 
                filter( iter == max(iter)) 
            #make sure we have both states
            assert_that(nrow(steps) == 2, msg=paste("something wrong with processing",file))
            
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
    }else{
        #get just the mcc score and add the other info
        mcc <- read_tsv(file)  %>% 
            select(mcc)  %>% 
            mutate(split=NA,
                   algorithm=algorithm,
                   subset=NA,
                   state=NA)
    }    
    return(mcc)
}

# map files to datafram
merged_mcc <- map_dfr(input,read_mcc)
write_csv(merged_mcc,snakemake@output[["merged_mcc"]])

