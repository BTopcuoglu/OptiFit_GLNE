##################################################################
# Author: Courtney Armour (adapted from Begum main.R script)
# Date: November 2021
# Description: Script to run mikropml pipeline
#
##################################################################

library(mikropml)
library(tidyverse)
library(caret)

# User defined variables
input <- commandArgs(trailingOnly=TRUE)

trainingProc <- input[1] # preprocessed data from training samples 
testingProc  <- input[2] # preprocessed data from testing samples 
model <- input[3] # Type of model to use
outcome <- input[4] # Classifaction to predict
nprocs <- as.numeric(input[5])

#for testing the pipeline
# method = "opticlust_denovo" # "opticlust_denovo","optifit_self","optifit_gg","vsearch_denovo","vsearch_gg"
# trainingProc <- paste0("results/ml/",method,"/preproc_train_split_1.csv")
# testingProc <- paste0("results/ml/",method,"/preproc_test_split_1.csv")
# model <- "glmnet"
# outcome <- "dx"
# nprocs <- 12

split <- unlist(str_split(trainingProc,"/|\\."))[4] %>% 
  str_replace(.,"preproc_t.+_(split_.+)","\\1")

print(paste0("training file: ",trainingProc))
print(paste0("testing file: ",testingProc))
print(paste0("model: ",model))
print(paste0("outcome: ",outcome))
print(paste0("split: ",split))
print(paste0("nproc: ",nprocs))

# Other variables
# if (str_detect(trainingProc, "optifit")) { # Setting output dir based on source of input
#   outDir <- "data/learning/results/optifit/"
# } else if(str_detect(trainingProc, "opticlust")){
#   outDir <- "data/learning/results/opticlust/"
# } else if(str_detect(trainingProc, "full")){
#   outDir <- "data/learning/results/gg_full/"
# } else if(str_detect(trainingProc, "subsample_8000")){
#   outDir <- "data/learning/results/gg_subsample_8000/"
# } else {
#   stop("Not sure where to output data")
# }
type <- unlist(str_split(trainingProc,"/|\\."))[3]
outDir <- paste0("results/ml/",type,"/")

print(paste0("writing output to: ",outDir))

#read in data
train <- read_csv(trainingProc)
test <- read_csv(testingProc)

# save ID and dx
trainIDS <- train %>% select(dx,Group) 
testIDS <- test  %>% select(dx,Group)

######################## RUN PIPELINE ###########################

#setup for parallelization
message(paste0("registering future plan: using ",nprocs," cores"))
doFuture::registerDoFuture()
future::plan(future::multicore, workers = nprocs)

# ----------------------------------------------------------------------->
message("PROGRESS: Running the model.")

# merge testing and training
# allData <- bind_rows(train %>% mutate(grps = rep("train",nrow(.))),
#                      test %>% mutate(grps = rep("test",nrow(.))))  %>% 
#     select(dx,grps,everything())
allData <- bind_rows(train,test)  %>% 
  select(-Group)
    
#grps <- c(rep("TR",nrow(train)),rep("TE",nrow(test)))
grps <- c(paste0("TR",seq(1,nrow(train))),
          paste0("TE",seq(1,nrow(test))))

#set seed to split number
sd <- as.numeric(str_replace(split,"split_",""))

#set hyperparameters
new_hp = list(mtry = c(30,68,137,274,350))

# Run the model
results <- mikropml::run_ml(dataset = allData,
                            method = model,
                            outcome_colname = outcome,
                            groups = grps,
                            #group_partitions = list(train = c("TR"),test = c("TE")),
                            group_partitions = list(train = c(paste0("TR",seq(1,nrow(train)))),
                                                    test = c(paste0("TE",seq(1,nrow(test))))),
                            seed = sd,
                            hyperparameters = new_hp)

# OUTPUTS ##########################################################

# write out performance results
performance <- results$performance
write_csv(performance,file = paste0(outDir, "performance_", split, ".csv"))

# write out prediction probabilities for sample
prediction <- predict(results$trained_model,test,type = "prob")  %>% 
  bind_cols(.,testIDS) %>% 
  select(Group,dx,cancer,normal)
write_csv(prediction,file = paste0(outDir, "prediction_", split, ".csv"))

#hyperparameter performance
hyperparameters <- get_hp_performance(results$trained_model)$dat
write_csv(hyperparameters,paste0(outDir,"hp_",split,".csv"))

# write out model
saveRDS(results$trained_model,file = paste0(outDir, "model_",split,".rds"))

###################################################################
