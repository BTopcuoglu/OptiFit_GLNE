################################
# Courtney R Armour
# August 2022
# Schloss Lab
# University of Michigan
################################

### SETUP -------------------
library(tidyverse)

outDir <- "data/learning/summary/"
if(!dir.exists(outDir)){dir.create(outDir)}

input <- commandArgs(trailingOnly = TRUE)

# input <- c("data/process/optifit/split_1/test/glne.precluster.pick.renamed.fit.optifit_mcc.0.03.subsample.shared",
#            "data/process/optifit/split_2/test/glne.precluster.pick.renamed.fit.optifit_mcc.0.03.subsample.shared")

### FUNCTIONS -----------------

get_frac_nonmapped <- function(file){
    file_parts <- unlist(str_split(file,"/"))  
    split <- file_parts[grepl(pattern="split",file_parts)]
      
    read_tsv(file) %>% 
        select(-label,-numOtus) %>% 
        mutate(RefSum=rowSums(select(.,starts_with("Ref_"))),
               OtherSum=rowSums(select(.,starts_with("Otu")))) %>% 
        select(Group,RefSum,OtherSum) %>% 
        mutate(total = RefSum+OtherSum) %>% 
        mutate(frac_nonmapped = OtherSum/total,
               split = split)
    
}

table <- map_dfr(input,get_frac_nonmapped)
write_csv(table,"results/tables/fracNonMapped.csv")