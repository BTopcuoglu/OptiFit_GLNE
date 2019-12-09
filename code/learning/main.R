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
model <- input[2]
outcome <- input[3]
sample_num <- as.numeric(input[4])



######################## DATA PREPARATION ########################
# Features: 16S rRNA gene sequences(OTUs) in the stool
# Labels: - Colorectal lesions of 490 patients.
#         - Defined as cancer or not.
#           (Cancer here means: SRN)
#           (SRNs are advanced adenomas+carcinomas)

shared <- read.delim('data/process/without.opti_mcc.0.03.2003650.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus)

# Read in OTU table and remove label and numOtus columns
shared_one_out <- read.delim('data/process/sample.2003650.optifit_mcc.0.03.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus)

colnames(shared_one_out) <- str_replace_all(colnames(shared_one_out), "Ref_", "")

labels_all <- colnames(shared)

labels_one <- colnames(shared_one_out)

common_cols <- intersect(labels_one, labels_all)

shared_one_out <- shared_one_out %>%
  select(common_cols)

missing_cols <- setdiff(labels_all, labels_one)

shared_one_out[missing_cols] <- 0

test <- shared_one_out

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


write.csv(walltime, file=paste0("data/temp/walltime_", model, "_", seed, ".csv"), row.names=F)
###################################################################
