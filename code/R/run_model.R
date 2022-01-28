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
# trainingProc <- "data/learning/results/opticlust/preproc_train_split_1.csv"
# testingProc <- "data/learning/results/opticlust/preproc_test_split_1.csv"
# model <- "glmnet"
# outcome <- "dx"
# nprocs <- 12

split <- unlist(str_split(trainingProc,"/|\\."))[5]  %>% 
    str_replace(.,"preproc_t.+_(split_.+)","\\1")

print(paste0("training file: ",trainingProc))
print(paste0("testing file: ",testingProc))
print(paste0("model: ",model))
print(paste0("outcome: ",outcome))
print(paste0("split: ",split))
print(paste0("nproc: ",nprocs))

# Other variables
if (str_detect(trainingProc, "optifit")) { # Setting output dir based on source of input
  outDir <- "data/learning/results/optifit/"
} else {
  outDir <- "data/learning/results/opticlust/"
}

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
new_hp = list(mtry = c(68,137,274,350))

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

performance <- results$performance

prediction <- predict(results$trained_model,test,type = "prob")  %>% 
  bind_cols(.,testIDS) %>% 
  select(Group,dx,cancer,normal)

#hyperparameter performance
hyperparameters <- get_hp_performance(results$trained_model)$dat
write_csv(hyperparameters,paste0(outDir,"hp_",split,".csv"))

# write out model
saveRDS(results$trained_model,file = paste0(outDir, "model_",split,".rds"))

# write out performance results
write_csv(performance,file = paste0(outDir, "performance_", split, ".csv"))

# write out prediction probabilities for sample
write_csv(prediction,file = paste0(outDir, "prediction_results_", split, ".csv"))

###################################################################
