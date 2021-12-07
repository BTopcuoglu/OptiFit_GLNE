library(tidyverse)
library(pROC)

### FUNCTIONS ###########################

get_pred_file_list <- function(algorithm){
  list.files(path=paste0("data/learning/results/",algorithm),
             pattern="prediction_results_split_[0-9]*.csv",
             full.names=T)
}

get_auc <- function(filename){
  alg <- str_replace(filename,"data/learning/results/(.*)/prediction_results_split_(\\d*).csv", "\\1")
  split <- str_replace(filename,"data/learning/results/(.*)/prediction_results_split_(\\d*).csv", "\\2")
  
  probs <- read_csv(filename,col_types = cols(Group = col_character(),
                                     dx = col_character(),
                                     .default = col_double())) 
  
  response <- probs$dx
  predictor <- probs$cancer
  
  split_auc <- as.numeric(auc(response, predictor,
                              levels=c("normal","cancer"),direction="<"))
  
  tibble(algorithm=alg,
         split=split,
         auc=split_auc)
}

### MAIN #################
algorithms <- c("optifit","opticlust")

allAUCs <- get_pred_file_list(algorithms) %>% 
  map_dfr(get_auc)

write_csv(allAUCs,"data/learning/summary/all_test_AUC.csv")