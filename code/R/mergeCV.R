################################
# Courtney R Armour
# November 2021
# Schloss Lab
# University of Michigan
################################

### SETUP -------------------
library(tidyverse)

outDir <- "data/learning/summary/"

input <- commandArgs(trailingOnly = TRUE)

### FUNCTIONS -------------------
read_cv <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    split <- file_parts[grepl(pattern="split",file_parts)]
    split <- gsub("cv_results_","",split)
    
    algorithm <- file_parts[grepl(pattern="^opticlust$|^optifit$",file_parts)]
    
    cv <- read_csv(file,na = c("","NA","NaN")) %>% #some of the Pos/Neg Pred Values were NaN 
        mutate(split=split,algorithm=algorithm)
  
  return(cv)
}

### MAIN -------------------

mergedCV <- map_dfr(input,read_cv)

write.csv(mergedCV,paste0(outDir,"merged_CV.csv"))