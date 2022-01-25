################################
# Courtney R Armour
# January 2022
# Schloss Lab
# University of Michigan
################################

### SETUP -------------------
library(tidyverse)

outDir <- "results/tables/"
if(!dir.exists(outDir)){dir.create(outDir)}

input <- snakemake@input[["sensspec"]]
#file = "data/process/opticlust/split_1/sub/glne.precluster.opti_mcc.0.03.pick.sensspec"

### FUNCTIONS -------------------
read_sensspec <- function(file){
    file_parts <- unlist(str_split(file,"/|\\."))
    
    split <- file_parts[grepl(pattern="split",file_parts)]
    split <- gsub("prediction_results_","",split)
    
    algorithm <- file_parts[grepl(pattern="^opticlust$|^optifit$",file_parts)]
    
    sensspec <- read_tsv(file,na = c("","NA","NaN")) %>% #some of the Pos/Neg Pred Values were NaN 
        select(mcc)  %>% 
        mutate(split=split,algorithm=algorithm)
  
  return(sensspec)
}

### MAIN -------------------

mergedMCC <- map_dfr(input,read_sensspec)

write_csv(mergedMCC,paste0(outDir,"opticlust_20_mcc.csv"))