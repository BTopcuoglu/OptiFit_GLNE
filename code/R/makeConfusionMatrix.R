#!/usr/bin/env Rscript
# makeConfusionMatrix.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Combines output from L2 logistic regression ML pipeline into a confusion matrix.
# Usage: Rscript makeConfusionMatrix.R CVFILE1 CVFILE2 CVFILE3 ... CVFILEN PREDFILE1 PREDFILE2 PREDFILE3 ... PREDFILEN DXDIFFTHRESH CLASSTHRESH

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
metadata <- args[1] # Sample metadata with diagnosis information
mlFiles <- sort(args[2:(length(args)-2)]) # Output files from ML pipeline (cvAUC and prediction scores)
dxDiffThresh <- args[(length(args)-1)] # Threshold for wanting to investigate health data because prediction scores are too close
classThresh <- args[length(args)] # Threshold for calling normal based on prediction values

# Other variables
outDir <- "data/learning/summary/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Analysis ----------------------------------------------------------------

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)



# Collating results -------------------------------------------------------

# Extracting sample names from file names
sample_names <- str_extract(mlFiles, "\\d+") %>% 
  unique()

# Reading in cvAUC files and adding col for sample name
cv <- str_subset(mlFiles, "cv") %>% 
  map_df(read_csv, col_types = cols()) %>% 
  add_column(sample = sample_names, .before = "cv_auc")

# Reading in prediction files and adding col for sample name
prediction <- str_subset(mlFiles, "prediction") %>% 
  map_df(read_csv, col_types = cols()) %>% 
  add_column(sample = sample_names, .before = "cancer")

# Importing metadata and recoding to either be normal = 0 or cancer = 1
meta <- read_tsv(metadata, col_types = cols()) %>% 
  mutate(sample = as.character(sample),
         meta_class = case_when(Dx_Bin == "Adenoma" ~ 0, # Recoding diagnoses
                                Dx_Bin == "Normal" ~ 0,
                                Dx_Bin == "High Risk Normal" ~ 0,
                                Dx_Bin == "adv Adenoma" ~ 1,
                                Dx_Bin == "Cancer" ~ 1)) %>% 
  select(sample, meta_class) # Removing unnecessary cols

# Creating results df
results <- left_join(cv, prediction, by = "sample") %>% 
  mutate(dx_diff = abs(cancer - normal), # Difference between prediction scores for dx
         dx_flag = case_when(dx_diff <= dxDiffThresh ~ 1, # Noting whether classifications need to verified by searching 
                          dx_diff > dxDiffThresh ~ 0),
         pred_class = case_when(normal >= classThresh ~ 0, # Determining classification using threshold of classThresh
                                    normal < classThresh ~ 1)) %>% 
  left_join(meta, by = "sample") # Adding metadata classification

# Writing out collated results
write_tsv(results, path = paste0(outDir, "model_results.tsv"))



# Determining confusion matrix --------------------------------------------

# Different categories for confusion matrix classification
# TN = true neg, TP = true pos, FN = false neg, FP = false pos
conf_labels <- tibble(conf = c("TN", "TP", "FN", "FP"))

# Calculating confusion matrix groups
conf <- results %>% 
  mutate(conf = case_when(meta_class == 0 & pred_class == 0 ~ "TN",
                          meta_class == 1 & pred_class == 1 ~ "TP",
                          meta_class > pred_class ~ "FN",
                          meta_class < pred_class ~ "FP")) %>% 
  group_by(conf) %>% 
  count() %>% # Adding up different categories
  right_join(conf_labels, by = "conf") %>% # Filling in any missing confusion matrix categories
  replace_na(list(n = 0)) # Replacing missing values with 0

# Writing out confusion matrix results
write_tsv(conf, path = paste0(outDir, "confusion_matrix.tsv"))
