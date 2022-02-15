# This script examines the split files to output a matrix where each row 
# is a sample and each column is a split. The 1s in the matrix represent which 
# splits that sample is in the test set and 0s are when that sample is in the 
# training set. For example: 
#         | split1 | split2 
#   -------------------------
#   samp1 |    1   |   0
#   samp2 |    0   |   1
# this shows that sample 1 is in the test set in split1 and the training set in split 2
# while samp2 is in the training set in split1 and the test set in split2.

library(tidyverse)
library(assertthat)

split_files <- snakemake@input[["splits"]]
#split_files <- list.files("./data/process/splits/",pattern="split_",full.names = T)

get_split_binary <- function(file){
    #get split number from file name
    split_num <- basename(file) %>%
        str_replace(".csv","")
    
    #read in the split file
    split_info <- read_csv(file) 
    
    #check file for all 490 samples and 2 columns
    assert_that(nrow(split_info) == 490)
    assert_that(ncol(split_info) == 2)
    
    binary <- split_info %>% 
        #mutate(!!split := case_when(train_test == "test" ~ 1,
        mutate(split = case_when(train_test == "test" ~ 1,
                                    train_test == "train" ~ 0,
                                    TRUE ~ as.numeric(NA)))  %>% 
        select(-train_test)
    
    # set samples as columns and split as row
    wide_binary <- binary  %>% 
        pivot_wider(names_from=group,values_from=split)  %>% 
        mutate(split = split_num)  %>% 
        select(split,everything())

    return(wide_binary)
}

#because map_dfc doesnt check rownames, I'm making the table where columns are 
#samples and rows are splits, merging all the splits, then flipping
split_mapping <- map_dfr(split_files,get_split_binary) 

split_mapping  %>% 
    pivot_longer(-split,'variable','value')  %>% 
    pivot_wider(variable,split)  %>% 
    rename(Group=variable)  %>% 
    write_csv(snakemake@output[["splitTogetherFreq"]])