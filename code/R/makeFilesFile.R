#!/usr/bin/env Rscript
# makeFilesFile.R
# William L. Close
# Pat Schloss Lab
# University of Michigan

# Purpose: Creates files file from SRA info for input into mothur.
# Usage: Rscript makeFilesFile.R SRAINFO SRASEQS1 SRASEQS2 SRASEQS3 ... SRASEQSN

# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sraInfo <- args[1] # SRA run metadata
sraSeqs <- args[-1] # Sequence files

# Checking if inputs are set
if (!file.exists(sraInfo)) {
  stop(paste(sraInfo, "does not exist."))
} 


# Other variables
outDir <- "data/process/"



# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Analysis ----------------------------------------------------------------

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE)

# Reading in the SRA run metadata
sra_info <- read_tsv(sraInfo, col_types = cols())

# Creating list sample names with paired files only
paired_samples <- sraSeqs %>% 
  str_extract("SRR\\d+") %>% 
  unique()

# Removing rows for samples with only one read file
sra_paired <- sra_info %>% 
  filter(Run %in% paired_samples)

# Creating the mothur files file
glne_files <- tibble(group = sra_paired$Sample_Name, # Sample name
                     run = sra_paired$Run) %>% # Sample run identifier
  mutate(forward = paste0(run, "_1.fastq.gz"), # Reconstructing name of forward seq file
         reverse = paste0(run, "_2.fastq.gz")) %>% # Reconstructing name of reverse seq file
  select(-run) # Dropping unneeded col


# Write out the files file to the mothur working dir
write_tsv(glne_files, path = paste0(outDir, "glne.files"), col_names = FALSE)
