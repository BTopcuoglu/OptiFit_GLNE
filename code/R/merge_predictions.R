################################
# Courtney R Armour
# November 2021
# Schloss Lab
# University of Michigan
################################

### SETUP -------------------
library(tidyverse)

outDir <- "data/learning/summary/"
if(!dir.exists(outDir)){dir.create(outDir)}

input <- commandArgs(trailingOnly = TRUE)
# input1 <- list.files(path="data/learning/results/opticlust/",pattern="prediction*",full.names=T)
# input2 <- list.files(path="data/learning/results/optifit/",pattern="prediction*",full.names=T)
# input <- c(input1,input2)

### FUNCTIONS -------------------
read_pred <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    split <- file_parts[grepl(pattern="split",file_parts)]
    split <- gsub("prediction_results_","",split)
    
    algorithm <- file_parts[grepl(pattern="^opticlust$|^optifit$",file_parts)]
    
    pred <- read_csv(file,na = c("","NA","NaN")) %>% #some of the Pos/Neg Pred Values were NaN 
        mutate(split=split,algorithm=algorithm)
  
  return(pred)
}

### MAIN -------------------

mergedPred <- map_dfr(input,read_pred)

write_csv(mergedPred,paste0(outDir,"merged_predictions.csv"))
