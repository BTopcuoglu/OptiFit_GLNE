##################################################################
# Author: Courtney Armour 
# Date: December 2021
# Description: Script to preprocess data to run in mikropml pipeline
#
##################################################################

library(mikropml)
library(tidyverse)
library(caret)

# User defined variables
input <- commandArgs(trailingOnly=TRUE)

metadata <- input[1] # Metadata containing classification to predict
training <- input[2] # Subsampled shared from training samples 
testing  <- input[3] # Subsampled shared from testing samples 
ncores   <- as.numeric(input[4])

#for testing the pipeline
# metadata <- "data/metadata/metadata.tsv"
# training <- "data/process/opticlust/split_1/train/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared"
# #training <- "data/process/optifit/split_5/train/glne.precluster.pick.opti_mcc.0.03.subsample.shared"
# testing  <- "data/process/opticlust/split_1/test/glne.precluster.opti_mcc.0.03.subsample.0.03.pick.shared"
# #testing <- "data/process/optifit/split_5/test/glne.precluster.pick.subsample.renamed.fit.optifit_mcc.shared"
# ncores=12

split <- unlist(str_split(training,"/"))[4]

print(paste0("metadata file: ",metadata))
print(paste0("training file: ",training))
print(paste0("testing file: ",testing))
print(paste0("split: ",split))

# Other variables
if (str_detect(training, "optifit")) { # Setting output dir based on source of input
  outDir <- "data/learning/results/optifit/"
} else {
  outDir <- "data/learning/results/opticlust/"
}

######################## DATA PREPARATION ########################
# Features: 16S rRNA gene sequences(OTUs) in the stool

message("PROGRESS: Preparing data.")

doFuture::registerDoFuture()
future::plan(future::multicore, workers = ncores)

# Reading in shared file from training data
train_shared <- read_tsv(training, col_types = cols(Group=col_character(),
                                                    .default = col_double())) %>%
  select(-label, -numOtus)

# Reading in shared file of the testing data
if (str_detect(testing, "optifit")) { 
  test_shared <- read_tsv(testing, col_types = cols(Group=col_character(),
                                                  .default = col_double())) %>%
    select(-label, -numRefOtus) %>%
    rename_all(str_replace, "Ref_", "")
} else {
  test_shared <- read_tsv(testing, col_types = cols(Group=col_character(),
                                                  .default = col_double())) %>%
    select(-label, -numOtus) 
}


# Create test set ---------------------------------------------------------

# Creating list of otus in loo shared
train_labels <- names(train_shared)

# Creating list of otus in opti shared
test_labels <- names(test_shared)

# Finding list of otus common to both shared files
common_labels <- intersect(train_labels,test_labels)

# Finding otus that are missing in sample shared
missing_cols <- setdiff(train_labels, test_labels)

# Only keeping otus that are in both files and assigning as test data for ML prediction
test <- test_shared %>%
  select(common_labels)

# Adding in missing otus all coded as 0 so test has the same cols as the loo shared
test[missing_cols] <- 0

test <- test[train_labels]

# Create training set -----------------------------------------------------

# Merge metadata and OTU table.
# Group advanced adenomas and cancers together as cancer and normal, 
# high risk normal and non-advanced adenomas as normal
# Then remove the sample ID column

# Read in metadata and select only sample Id, diagnosis, and fit columns
meta <- read_tsv(metadata, col_types = cols(sample=col_character())) %>%
  select(sample, Dx_Bin, fit_result)

# Combining metadata and training data
data <- inner_join(train_shared, meta, by=c("Group"="sample")) %>%
  mutate(dx = case_when(Dx_Bin == "Adenoma" ~ "normal", # Recoding diagnoses
                        Dx_Bin == "Normal" ~ "normal",
                        Dx_Bin == "High Risk Normal" ~ "normal",
                        Dx_Bin == "adv Adenoma" ~ "cancer",
                        Dx_Bin == "Cancer" ~ "cancer",
                        TRUE ~ NA_character_),
         dx = as.factor(dx)) %>% # Encoding dx as factor
  #select(-Group, -Dx_Bin, -fit_result) %>%
  select(-Dx_Bin, -fit_result) %>%
  drop_na() %>%
  select(dx, everything()) %>%
  as.data.frame()
  
test <- inner_join(test, meta, by=c("Group"="sample"))  %>% 
    mutate(dx = case_when(Dx_Bin == "Adenoma" ~ "normal", # Recoding diagnoses
                        Dx_Bin == "Normal" ~ "normal",
                        Dx_Bin == "High Risk Normal" ~ "normal",
                        Dx_Bin == "adv Adenoma" ~ "cancer",
                        Dx_Bin == "Cancer" ~ "cancer",
                        TRUE ~ NA_character_),
         dx = as.factor(dx)) %>% # Encoding dx as factor
  #select(-Group, -Dx_Bin, -fit_result) %>%
  select(-Dx_Bin, -fit_result) %>%
  drop_na() %>%
  select(dx, everything()) %>%
  as.data.frame()

# check that there is no overlap in samples between test and train 
if( length(intersect(data$Group,test$Group)) != 0 ){
  stop("something went wrong with the data split, there are samples in testing and training.")
}

#remove group from data/test
# data <- data %>% 
#   select(-Group)
  
# test <- test %>% 
#   select(-Group)

######################## RUN PIPELINE ###########################

message("PROGRESS: Pre-processing data.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

set.seed(1989)

# ------------------Pre-process the full Dataset------------------------->
# We are doing the pre-processing to the training dataset,
# then applying the preprocessing to the test dataset

preProcValues <- preProcess(data, method = c("center", "scale"))

dataTransformed <- predict(preProcValues, data)

testTransformed <- predict(preProcValues, test)

write_csv(dataTransformed,file = paste0(outDir, "preproc_train_", split,".csv"))
write_csv(testTransformed,file = paste0(outDir, "preproc_test_", split,".csv"))
