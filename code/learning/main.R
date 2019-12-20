######################################################################
# Author: Begum Topcuoglu
# Date: 2018-12-20
# Title: Main pipeline for 7 classifiers in R programming language
###################################################################

###################################################################
# Description:

# This script will read in data from Baxter et al. 2016
#     - 0.03 subsampled OTU dataset
#     - CRC metadata: SRN information


# It will run the following machine learning pipelines:
#     - L2 Logistic Regression
#     - L1 and L2 Linear SVM
#     - RBF SVM
#     - Decision Tree
#     - Random Forest
#     - XGBoost
###################################################################

###################################################################
# Dependencies and Outputs:

# Be in the project directory.

# The outputs are:
#   (1) AUC values for cross-validation and testing for each data-split
#   (2) meanAUC values for each hyper-parameter tested during each split.
###################################################################


################### IMPORT LIBRARIES and FUNCTIONS ###############
# The dependinces for this script are consolidated in the first part
deps = c("dplyr", "tictoc", "caret" ,"rpart", "xgboost", "randomForest", "kernlab","LiblineaR", "pROC", "tidyverse");
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE, repos = "http://cran.us.r-project.org", dependencies=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
# Load in needed functions and libraries
source('code/learning/model_selection.R')
source('code/learning/model_pipeline_deployed.R') # has pipeline function defined here
##################################################################


# We will run main.R from command line with arguments
#  - These arguments will be saved into variable "input"
#  - First argument is the seed number which is the array index
#  - Second argument is the model name (one of the list above)
#  - Third argument is the perform permutation importance or not
#  - Fourth argument is the outcome we want to predict:
#        "dx"
#  - Fifth argument is the sample number we are predicting

input <- commandArgs(trailingOnly=TRUE)

outSubShared <- input[1] # Subsampled shared from samples after leaving one out
optifitSubShared <- input[2] # Subsampled shared after OptiFit clustering of left out sample
metadata <- input[3] # Metadata containing classification to predict
model <- input[4] # Type of model to use
outcome <- input[5] # Classifaction to predict

# sample_num <- as.numeric(input[3])
# 'data/process/without.opti_mcc.0.03.2003650.subsample.shared'
# 'data/process/sample.2003650.optifit_mcc.0.03.subsample.shared'
# 'data/process/metadata.tsv'





######################## DATA PREPARATION ########################
# Features: 16S rRNA gene sequences(OTUs) in the stool
# Labels: - Colorectal lesions of 490 patients.
#         - Defined as cancer or not.
#           (Cancer here means: SRN)
#           (SRNs are advanced adenomas+carcinomas)

# Reading in shared file from LOO
loo_shared <- read_tsv(outSubShared) %>%
  select(-label, -numOtus)

# Creating list of otus in loo shared
loo_labels <- names(loo_shared)



# Reading in shared file from OptiFit
opti_shared <- read_tsv(optifitSubShared) %>%
  select(-label, -numOtus) %>% 
  rename(str_replace, "Ref_", "")

# Creating list of otus in opti shared
opti_labels <- names(shared_one_out)



# Finding list of otus common to both shared files
common_cols <- intersect(opti_labels, loo_labels)



# Only keeping otus that are in both files and assigning as test data for ML prediction
test <- opti_shared %>%
  select(common_cols)

# Not needed because missing_cols aren't in opti_shared after using select
# # Finding list of otus not in one of the shared files
# missing_cols <- setdiff(loo_labels, opti_labels)
# # Recoding any otus not in loo_shared as all 0s
# opti_shared[missing_cols] <- 0
# test <- opti_shared




# Merge metadata and OTU table.
# Group advanced adenomas and cancers together as cancer and normal, high risk normal and non-advanced adenomas as normal
# Then remove the sample ID column

# Read in metadata and select only sample Id and diagnosis columns
meta <- read.delim('data/process/metadata.tsv', header=T, sep='\t') %>%
  select(sample, Dx_Bin, fit_result) %>%
  filter(sample != "2003650")

data <- inner_join(meta, shared, by=c("sample"="Group")) %>%
  mutate(dx = case_when(
    Dx_Bin== "Adenoma" ~ "normal",
    Dx_Bin== "Normal" ~ "normal",
    Dx_Bin== "High Risk Normal" ~ "normal",
    Dx_Bin== "adv Adenoma" ~ "cancer",
    Dx_Bin== "Cancer" ~ "cancer"
  )) %>%
  select(-sample, -Dx_Bin, -fit_result) %>%
  drop_na() %>%
  select(dx, everything())
# We want the diagnosis column to be a factor
data$dx <- factor(data$dx)


######################## RUN PIPELINE ###########################
# Choose which classification methods we want to run on command line
#                "L2_Logistic_Regression",
#                "Random_Forest",

set.seed(1989)

# Run the model
  # User can define the outcome and to do permutation or not here:
  # example: get_results(data, model, seed, 0, "dx)
  # OR pass as NA
results <- pipeline(data, test, model, outcome)

cv_auc <- results[1]
prediction <- results[2]

# Create a matrix with cv_aucs and test_aucs from 1 data split
aucs <- matrix(results[[1]], ncol=1)
# Convert to dataframe and add a column noting the model name
aucs_dataframe <- data.frame(aucs) %>%
  rename_at(1, ~ "cv_auc") %>%
  write_csv(path = paste0("data/temp/cv_results_", sample_num, ".csv"))

# Convert to dataframe and add a column noting the model name
predictions_dataframe <- data.frame(prediction) %>%
    write_csv(path = paste0("data/temp/prediction_results_", sample_num, ".csv"))

###################################################################
