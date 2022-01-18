################################
# Courtney R Armour
# Janurary 2022
# Schloss Lab
# University of Michigan
################################

### SETUP -------------------
library(tidyverse)

outDir <- "data/learning/summary/"
if(!dir.exists(outDir)){dir.create(outDir)}

input <- commandArgs(trailingOnly = TRUE)
# input1 <- list.files(path="data/learning/results/opticlust/",pattern="performance*",full.names=T)
# input2 <- list.files(path="data/learning/results/optifit/",pattern="performance*",full.names=T)
# input <- c(input1,input2)

### FUNCTIONS -------------------
read_perf <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    split <- file_parts[grepl(pattern="split",file_parts)]
    split <- gsub("performance_","",split)
    
    algorithm <- file_parts[grepl(pattern="^opticlust$|^optifit$",file_parts)]
    
    perf <- read_csv(file,na = c("","NA","NaN")) %>% #some of the Pos/Neg Pred Values were NaN 
        mutate(split=split,algorithm=algorithm)
  
  return(perf)
}

### MAIN -------------------

mergedPerf <- map_dfr(input,read_perf)

write_csv(mergedPerf,paste0(outDir,"merged_performance.csv"))
