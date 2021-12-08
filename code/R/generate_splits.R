library(tidyverse)

### FUNCTIONS #############
create_split <- function(groups,train_frac,split_num){
    
    groups_shuffled <- sample(groups,length(groups))
    
    train_n <- as.integer(length(groups)*train_frac)
    test_n <- length(groups) - train_n
    
    data <- tibble(group=groups_shuffled,
                   train_test=c(rep("train",train_n),
                                rep("test",test_n)))  %>% 
        arrange(group)
    
    return(data)
}

### VARIABLES ################
outdir <- "data/process/splits/"
if(!dir.exists(outdir)){dir.create(outdir)}

# Parsing command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Variables defined by user
shared <- read_tsv(args[1],col_types = cols(Group=col_character(),
                                            .default=col_double()))

n_splits <- args[2]

### MAIN ########################
# get the list of sampels from the shared file
groups <- shared  %>% 
    pull(Group)

#generate n random 80/20 splits, write files to outdir
set.seed(2021)
for(i in 1:n_splits){
    create_split(groups,0.8,i)  %>% 
        write_csv(.,paste0(outdir,"split_",i,".csv"))
}
