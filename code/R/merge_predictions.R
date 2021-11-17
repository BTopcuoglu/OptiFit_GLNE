##################################################################
# Author: Courtney Armour 
# Date: November 2021
# Description: merge prediction results from run_model.R
##################################################################

library(tidyverse)

outDir <- "data/learning/summary/"

# get files from command line
input <- commandArgs(trailingOnly=TRUE)

#separate optifit and opticlust files
opticlust_files <- input[grepl(pattern="opticlust",input)]
optifit_files <- input[grepl(pattern="optifit",input)]

print("opticlust files: ")
print(opticlust_files)
print("optifit files: ")
print(optifit_files)

#################################
### FUNCTIONS -------------------
#################################
read_files <- function(file){
    #get split and algorithm info from filename
    file_parts <- unlist(str_split(file,"/|\\."))
    split <- file_parts[grepl(pattern="split",file_parts)] 
    split <- gsub(pattern="prediction_results_",replacement="",split)
    algorithm <- file_parts[grepl(pattern="opti",file_parts)]
    
    pred <- read_csv(file) %>% 
        mutate(split=split,algorithm=algorithm)
        
    return(pred)
}


#################################
### MAIN -------------------
#################################

# merge opticlust files
opticlust_pred <- map_dfr(opticlust_files, read_files) 

# merge optifit files
optifit_pred <- map_dfr(optifit_files, read_files) 

# merge all files together
all_pred <- bind_rows(opticlust_pred,optifit_pred)

write_csv(all_pred,paste0(outDir,"merged_predictions.csv"))