#!/usr/bin/env Rscript
# makeConfusionMatrix.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Combines output from L2 logistic regression ML pipeline into a confusion matrix.
# Usage: Rscript makeConfusionMatrix.R CVFILE1 CVFILE2 CVFILE3 ... CVFILEN PREDFILE1 PREDFILE2 PREDFILE3 ... PREDFILEN

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
mlFiles <- args # SRA run metadata

# Other variables
outDir <- "results/tables/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Analysis ----------------------------------------------------------------

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

mlFiles <- sort(list.files("data/learning/", full.names = T))
metadata <- "data/test_process/metadata.tsv"
dxDiffThresh <- 0.05 # Threshold for wanting to investigate health data because prediction scores are too close
classThresh <- 0.5 # Threshold for calling normal based on prediction values


# Importing metadata and recoding to either be normal = 0 or cancer = 1
meta <- read_tsv(metadata, col_types = cols()) %>% 
  mutate(sample = as.character(sample),
         meta_class = case_when(Dx_Bin == "Adenoma" ~ 0, # Recoding diagnoses
                        Dx_Bin == "Normal" ~ 0,
                        Dx_Bin == "High Risk Normal" ~ 0,
                        Dx_Bin == "adv Adenoma" ~ 1,
                        Dx_Bin == "Cancer" ~ 1)) %>% 
  select(sample, meta_class)


# Extracting sample names from file names
sample_names <- str_extract(mlFiles, "\\d+") %>% 
  unique()

# Subsetting out cross validated AUC score files
cv_files <- str_subset(mlFiles, "cv")

# Subsetting out prediction score files
prediction_files <- str_subset(mlFiles, "prediction")



# Reading in cvAUC files and adding col for sample name
cv <- map_df(cv_files, read_csv, col_types = cols()) %>% 
  add_column(sample = sample_names, .before = "cv_auc")

# Reading in prediction files and adding col for sample name
prediction <- map_df(prediction_files, read_csv, col_types = cols()) %>% 
  add_column(sample = sample_names, .before = "cancer")




# Creating classification df
classification <- left_join(cv, prediction, by = "sample") %>% 
  mutate(dx_diff = abs(cancer - normal), # Difference between prediction scores for dx
         dx_flag = case_when(dx_diff <= 0.05 ~ 1, # Noting whether classifications need to verified by searching 
                          dx_diff > 0.05 ~ 0),
         pred_class = case_when(normal >= classThresh ~ 0, # Determining classification using threshold of classThresh
                                    normal < classThresh ~ 1)) %>% 
  left_join(meta, by = "sample")


