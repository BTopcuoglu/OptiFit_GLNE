library(tidyverse)

get_training_feats <- function(method,split){
    training <- paste0("data/process/",method,"/split_",split,"/train/glne.",method,".shared")
    
    train_preproc <- paste0("results/ml/",method,"/preproc_train_split_",split,".csv")
        
    print(paste0("training file: ",training))
    
    #get training features after preprocessing
    train_preproc_feats <- names(read.table(file = train_preproc,header = T,nrows = 1,sep=","))[c(-1,-2)]
    
    #Reading in shared file from training data
    if (method == "optifit_gg"){
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numRefOtus)  %>% 
        # remove de novo OTUs - only keep reference OTUs
        select(Group,starts_with("Ref_"))
    } else if (method == "vsearch_denovo") {
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        label=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numASVs)
    } else if (method == "vsearch_gg") {
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        label=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numRefOtus)  %>% 
        # remove de novo OTUs - only keep reference OTUs
        select(Group,starts_with("Ref_"))
    } else if (method %in% c("opticlust_denovo", "optifit_self")){
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numOtus)
    } else {
    stop(paste("unrecognized method",method))
    }
    
   
    
    #get the names of the features in the train set
    train_shared_feats <- names(train_shared  %>%  select(-Group)) 
    #subset to just the features that make it through preprocessing
    final_train_feats <- train_shared_feats[train_shared_feats %in% train_preproc_feats]
    
    return(final_train_feats)
}

get_training_counts <- function(method,split){
    training <- paste0("data/process/",method,"/split_",split,"/train/glne.",method,".shared")
    
    train_preproc <- paste0("results/ml/",method,"/preproc_train_split_",split,".csv")
        
    print(paste0("training file: ",training))
    
    #get training features after preprocessing
    train_preproc_feats <- names(read.table(file = train_preproc,header = T,nrows = 1,sep=","))[c(-1,-2)]
    
    #Reading in shared file from training data
    if (method == "optifit_gg"){
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numRefOtus)  %>% 
        # remove de novo OTUs - only keep reference OTUs
        select(Group,starts_with("Ref_"))
    } else if (method == "vsearch_denovo") {
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        label=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numASVs)
    } else if (method == "vsearch_gg") {
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        label=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numRefOtus)  %>% 
        # remove de novo OTUs - only keep reference OTUs
        select(Group,starts_with("Ref_"))
    } else if (method %in% c("opticlust_denovo", "optifit_self")){
    train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                        .default = col_double())) %>%
        select(-label, -numOtus)
    } else {
    stop(paste("unrecognized method",method))
    }
    
   
    
    #get the names of the features in the train set
    train_shared_feats <- names(train_shared  %>%  select(-Group)) 
    #subset to just the features that make it through preprocessing
    final_train_feats <- train_shared_feats[train_shared_feats %in% train_preproc_feats]
    
    sub_train_shared <- train_shared %>% 
        select(Group,all_of(final_train_feats)) 
    
    counts <- sub_train_shared %>% 
        mutate(count=rowSums(across(where(is.numeric))),
               set="train") %>% 
        select(Group,set,count) 
        
    return(counts)
}

get_testing_counts <- function(method,split){
    testing <- paste0("data/process/",method,"/split_",split,"/test/glne.",method,".shared")
 
    print(paste0("testing file: ",testing))
    
    # Reading in shared file of the testing data
    if (method == "optifit_self") { 
    test_shared <- read_tsv(testing, col_types = cols(Group=col_character(),
                                                    .default = col_double())) %>%
        select(Group,starts_with("Ref_")) %>% #select only OTUs in reference
        rename_all(str_replace, "Ref_", "") #remove Ref_ label to match train data
    } else if (method == "optifit_gg"){
    test_shared <- read_tsv(testing, col_types = cols(Group=col_character(),
                                                    .default = col_double())) %>%
        select(Group,starts_with("Ref_")) #select only OTUs in reference
    } else if (method == "vsearch_denovo") {
    test_shared <- read_tsv(testing, col_types = cols(Group=col_character(),
                                                        label=col_character(),
                                                        .default = col_double())) %>%
        select(Group,starts_with("ASV"))

    } else if (method == "vsearch_gg") {
    test_shared <- read_tsv(testing, col_types = cols(Group=col_character(),
                                                        label=col_character(),
                                                        .default = col_double())) %>%
        select(Group,starts_with("Ref_")) #select only OTUs in reference
    } else if (method == "opticlust_denovo"){
    test_shared <- read_tsv(testing, col_types = cols(Group=col_character(),
                                                    .default = col_double())) %>%
        select(-label, -numOtus) 
    } else {
    stop(paste("unrecognized method",method))
    }
    
    #get names fo features in the test set
    test_feats <- names(test_shared  %>% select(-Group))
    
    train_feats <- get_training_feats(method,split)
    
    #subset to features that are in the training set
    sub_test_feats <- test_feats[test_feats %in% train_feats]
    
    counts <- test_shared %>% 
        select(Group,all_of(sub_test_feats)) %>% 
        mutate(count=rowSums(across(where(is.numeric))),
               set="test") %>% 
        select(Group,set,count) 
    
    return(counts)
}

get_counts <- function(method,split){
    
    train_counts <- get_training_counts(method,split)
    test_counts <- get_testing_counts(method,split)
    
    return(bind_rows(train_counts,
                     test_counts))    
}

methods <- c("opticlust_denovo","optifit_self","optifit_gg","vsearch_denovo","vsearch_gg")
splits <- seq(1,100)

data <- expand_grid(methods,splits)  %>% 
    group_by(methods,splits) %>% 
    summarize(get_counts(methods,splits),
              .groups = "drop")

write_csv(data,"results/tables/counts.csv")
