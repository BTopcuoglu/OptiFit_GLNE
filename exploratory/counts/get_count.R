#!/usr/bin/env Rscript
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

#file <- read_tsv("/Users/armourc/Desktop/opticlust/data/process/optifit/2003650/in/glne.precluster.pick.subsample.optifit_mcc.shared")
file <- read_tsv(args[1])

group <- file$Group

count <- file %>% 
  select(-label,-Group,-numOtus) %>% 
  rowSums()

results <- as.data.frame(t(c(group,count)))

write_csv(results,paste0("./",group,"_count.csv"),col_names = FALSE)
