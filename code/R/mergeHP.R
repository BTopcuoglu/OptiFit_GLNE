################################
# Courtney R Armour
# January 2022
# Schloss Lab
# University of Michigan
################################

### SETUP -------------------
library(tidyverse)

outDir <- "data/learning/summary/"
if(!dir.exists(outDir)){dir.create(outDir)}

input <- commandArgs(trailingOnly = TRUE)
input1 <- list.files(path="data/learning/results/opticlust/",pattern="hp*",full.names=T)
input2 <- list.files(path="data/learning/results/optifit/",pattern="hp*",full.names=T)
input <- c(input1,input2)

### FUNCTIONS -------------------
read_hp <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    split <- file_parts[grepl(pattern="split",file_parts)]
    split <- gsub("hp_","",split)
    
    algorithm <- file_parts[grepl(pattern="^opticlust$|^optifit$",file_parts)]
    
    hp <- read_csv(file,na = c("","NA","NaN")) %>% 
        mutate(split=split,algorithm=algorithm)
  
  return(hp)
}

### MAIN -------------------

mergedHP <- map_dfr(input,read_hp)

write_csv(mergedHP,paste0(outDir,"merged_HP.csv"))
