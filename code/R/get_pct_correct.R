library(tidyverse)

outDir <- "results/tables/"

### FUNCTIONS ##########################################
get_pred_file_list <- function(algorithm){
  list.files(path=paste0("data/learning/results/",algorithm),
             pattern="prediction_results_split_[0-9]*.csv",
             full.names=T)
}

get_categories <- function(filename,threshold){
  read_csv(filename,col_types = cols(Group=col_character(),
                                     dx=col_character(),
                                     .default=col_double())) %>%
    mutate(algorithm = str_replace(filename,
                                   "data/learning/results/(.*)/prediction_results_split_(\\d*).csv", "\\1"),
           split = str_replace(filename,
                               "data/learning/results/(.*)/prediction_results_split_(\\d*).csv", "\\2"),
           pred_dx = case_when(cancer > threshold ~ "cancer",
                               TRUE ~ "normal")) %>%
    select(Group,dx,pred_dx,algorithm,split) %>%
    mutate(category = case_when(dx == "normal" & pred_dx == "normal" ~ "TN",
                                dx == "normal" & pred_dx == "cancer" ~ "FP",
                                dx == "cancer" & pred_dx == "cancer" ~ "TP",
                                dx == "cancer" & pred_dx == "normal" ~ "FN",
                                TRUE ~ "NA"))
}

get_pct_correct <- function(algorithm, threshold){
  get_pred_file_list(algorithm) %>%
    map_dfr(get_categories,threshold) %>%
    select(Group,category,dx) %>%
    group_by(Group,dx) %>%
    count(category) %>%
    pivot_wider(id_cols=c("Group","dx"),names_from=category,values_from=n,values_fill=0) %>%
    mutate(correct = TN + TP,
           incorrect = FN + FP,
           total = TN + TP + FN + FP) %>%
    mutate(pct_correct = correct/total) %>%
    select(Group,dx,pct_correct) %>%
    mutate(algorithm = algorithm)
}

### MAIN ##############################
algorithms=c("opticlust","optifit")

pct_correct <- map_dfr(algorithms, get_pct_correct, threshold=0.5) %>%
  mutate(Group=as.character(Group))
  
write_csv(pct_correct,paste0(outDir,"pct_class_correct.csv"))


### TEST ########
# threshold = 0.5

# ex1 <- tibble(Group=rep("A",5),
#                dx=rep("normal",5),
#                cancer=c(0.1,0.2,0.3,0.4,0.45),
#                normal=c(1-cancer),
#                algorithm=rep("opticlust",5),
#                split=seq(1,5))  %>% 
#     mutate(pred_dx = case_when(cancer > threshold ~ "cancer",
#                                TRUE ~ "normal"))
# ex2 <- tibble(Group=rep("B",5),
#                dx=rep("normal",5),
#                cancer=c(0.1,0.2,0.45,0.55,0.6),
#                normal=c(1-cancer),
#                algorithm=rep("optifit",5),
#                split=seq(1,5))  %>% 
#     mutate(pred_dx = case_when(cancer > threshold ~ "cancer",
#                                TRUE ~ "normal"))  
# ex3 <- tibble(Group=rep("C",5),
#                dx=rep("cancer",5),
#                cancer=c(0.55,0.6,0.7,0.8,0.9),
#                normal=c(1-cancer),
#                algorithm=rep("opticlust",5),
#                split=seq(1,5))  %>% 
#     mutate(pred_dx = case_when(cancer > threshold ~ "cancer",
#                                TRUE ~ "normal"))
# ex4 <- tibble(Group=rep("D",5),
#                dx=rep("cancer",5),
#                cancer=c(0.4,0.45,0.7,0.8,0.9),
#                normal=c(1-cancer),
#                algorithm=rep("optifit",5),
#                split=seq(1,5))  %>% 
#     mutate(pred_dx = case_when(cancer > threshold ~ "cancer",
#                                TRUE ~ "normal"))

# data <- bind_rows(ex1,ex2,ex3,ex4)

# pct <- data  %>% 
#     select(Group,dx,pred_dx,algorithm,split) %>%
#     mutate(category = case_when(dx == "normal" & pred_dx == "normal" ~ "TN",
#                                 dx == "normal" & pred_dx == "cancer" ~ "FP",
#                                 dx == "cancer" & pred_dx == "cancer" ~ "TP",
#                                 dx == "cancer" & pred_dx == "normal" ~ "FN",
#                                 TRUE ~ "NA"))  %>% 
#     select(Group,algorithm,category,dx) %>%
#     group_by(Group,dx,algorithm) %>%
#     count(category) %>%
#     pivot_wider(id_cols=c("Group","dx","algorithm"),names_from=category,values_from=n,values_fill=0) %>%
#     mutate(correct = TN + TP,
#            incorrect = FN + FP,
#            total = TN + TP + FN + FP) %>%
#     mutate(pct_correct = correct/total) %>%
#     select(Group,dx,algorithm,pct_correct)
