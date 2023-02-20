################################
# Courtney R Armour
# 
# Schloss Lab
# University of Michigan
################################

library(tidyverse)

infiles = snakemake@input[["results"]]

process_file <- function(file){
    
    file_parts = unlist(str_split(file,"/|\\."))

    algorithm = file_parts[3]
    metric = gsub("_split_[0-9]+","",file_parts[4])
    split =  gsub(".*_(split_[0-9]+)","\\1",file_parts[4])
    
    dat <- read_csv(file,na = c("","NA","NaN")) %>% 
        mutate(algorithm=algorithm,result_metric=metric,split=split)
        
    return(dat)
}

merged <- map_dfr(infiles,process_file)

result_metric <- merged %>% 
    pull(result_metric) %>% unique()

if(length(result_metric) == 1){
    write_csv(merged %>% select(-result_metric),
              paste0("results/ml/summary/merged_",result_metric,".csv"))
} else{
    stop("something went wrong...")
}

