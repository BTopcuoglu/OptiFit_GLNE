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
deps <- c("tictoc", "caret" ,"rpart", "xgboost", "randomForest", "kernlab","LiblineaR", "pROC", "tidyverse")

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

# User defined variables
input <- commandArgs(trailingOnly=TRUE)

looSubShared <- input[1] # Subsampled shared from samples after leaving one out
optifitSubShared <- input[2] # Subsampled shared after OptiFit clustering of left out sample
metadata <- input[3] # Metadata containing classification to predict
model <- input[4] # Type of model to use
outcome <- input[5] # Classifaction to predict

# Other variables
outDir <- "data/learning/"
sampleNum <- str_extract(looSubShared, "\\d+") # Pulling sample number from file path



######################## DATA PREPARATION ########################
# Features: 16S rRNA gene sequences(OTUs) in the stool
# Labels: - Colorectal lesions of 490 patients.
#         - Defined as cancer or not.
#           (Cancer here means: SRN)
#           (SRNs are advanced adenomas+carcinomas)

message("PROGRESS: Preparing data.")

# Reading in shared file from LOO
loo_shared <- read_tsv(looSubShared, col_types = cols()) %>%
  select(-label, -numOtus)

# Reading in shared file from OptiFit
opti_shared <- read_tsv(optifitSubShared, col_types = cols()) %>%
  select(-label, -numOtus) %>% 
  rename_all(str_replace, "Ref_", "")



# Create test set ---------------------------------------------------------

# Creating list of otus in loo shared
loo_labels <- names(loo_shared)

# Creating list of otus in opti shared
opti_labels <- names(opti_shared)

# Finding list of otus common to both shared files
common_labels <- intersect(opti_labels, loo_labels)

# Finding otus that are missing in optifit shared
missing_cols <- setdiff(loo_labels, opti_labels)

# Only keeping otus that are in both files and assigning as test data for ML prediction
test <- opti_shared %>%
  select(common_labels)

# Adding in missing otus all coded as 0 so test has the same cols as the loo shared
test[missing_cols] <- 0



# Create training set -----------------------------------------------------

# Merge metadata and OTU table.
# Group advanced adenomas and cancers together as cancer and normal, high risk normal and non-advanced adenomas as normal
# Then remove the sample ID column

# Read in metadata and select only sample Id, diagnosis, and fit columns
meta <- read_tsv(metadata, col_types = cols()) %>%
  select(sample, Dx_Bin, fit_result)


# Combining metadata and loo shared file to use as training data
data <- inner_join(loo_shared, meta, by=c("Group"="sample")) %>%
  mutate(dx = case_when(Dx_Bin == "Adenoma" ~ "normal", # Recoding diagnoses
                        Dx_Bin == "Normal" ~ "normal",
                        Dx_Bin == "High Risk Normal" ~ "normal",
                        Dx_Bin == "adv Adenoma" ~ "cancer",
                        Dx_Bin == "Cancer" ~ "cancer",
                        TRUE ~ NA_character_),
         dx = as.factor(dx)) %>% # Encoding dx as factor
  select(-Group, -Dx_Bin, -fit_result) %>%
  drop_na() %>%
  select(dx, everything())



######################## RUN PIPELINE ###########################
# Choose which classification methods we want to run on command line
#                "L2_Logistic_Regression",
#                "Random_Forest",

message("PROGRESS: Predicting outcome.")

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

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
  write_csv(path = paste0("data/temp/cv_results_", sampleNum, ".csv"))

# Convert to dataframe and add a column noting the model name
predictions_dataframe <- data.frame(prediction) %>%
    write_csv(path = paste0("data/temp/prediction_results_", sampleNum, ".csv"))

###################################################################
