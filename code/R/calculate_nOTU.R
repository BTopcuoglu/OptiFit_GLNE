library(tidyverse)
library(data.table)

input <- snakemake@input %>% 
    unlist(., use.names=FALSE) %>% 
    unique()

# input <- c(paste0("results/ml/opticlust_denovo/preproc_train_split_",seq(1,5),".csv"),
#            paste0("results/ml/optifit_self/preproc_train_split_",seq(1,5),".csv"),
#            paste0("results/ml/optifit_gg/preproc_train_split_",seq(1,5),".csv"),
#            paste0("results/ml/vsearch_denovo/preproc_train_split_",seq(1,5),".csv"),
#            paste0("results/ml/vsearch_gg/preproc_train_split_",seq(1,5),".csv"))

calc_n_feats <- function(file){
    algorithm = gsub("results/ml/(.+)/preproc_train_split_([0-9]+).csv","\\1",file)
    split = gsub("results/ml/(.+)/preproc_train_split_([0-9]+).csv","\\2",file)
    
    n_feats <- as_tibble(fread(file,sep=",")) %>% 
        select(-dx,-Group) %>% ncol()
        
    return(tibble(algorithm=algorithm,
                  split=split,
                  n_feats=n_feats))
}

data <- map_dfr(input,calc_n_feats)

write_csv(data,snakemake@output[["outfile"]])