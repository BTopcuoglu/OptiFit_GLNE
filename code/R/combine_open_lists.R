library(tidyverse)

# Input and Output filenames
#if exists('snakemake') {
    infilename1 <- snakemake@input[[1]]
    infilename2 <- snakemake@input[[2]]
    outfilename <- snakemake@output[[1]]
# } else {
#     infilename1 = "data/process/vsearch_gg/glne.closed.list"
#     infilename2 = "data/process/vsearch_gg/glne.de_novo.list"
#     outfilename = "data/process/vsearch_gg/glne.vsearch_gg.list"
#}

print(infilename1)
print(infilename1)
print(outfilename)

#function to process file, re-name OTUs based on if they are reference(Ref) or denovo(new)
process_file <- function(filename){
    dat <- read_tsv(filename,col_names=FALSE)  %>% 
        rename(label=X1,
               numOtus=X2)

    label <- dat %>% pull(label)
    num <- dat %>% pull(numOtus)
    
    if(num > 9999){
        stop("Need to adjust code for more than 9999 OTUs.")
    }

    if(label == "Ref"){
        names(dat) <- c("label","numOtus",paste0("Ref_Otu",str_pad(seq(1,num),4,pad="0")))
    }else if(label == "new"){
        names(dat) <- c("label","numOtus",paste0("Otu",str_pad(seq(1,num),4,pad="0")))
    }else{
        stop(paste0("Unrecognized label ",label))
    }    
    return(dat %>% mutate(label = "userLabel"))
}

#process both input files
file1 <- process_file(infilename1)
file2 <- process_file(infilename2)

#merge the processed files
merge <- bind_cols(file1 %>% select(-label,-numOtus),
          file2 %>% select(-label,-numOtus))  %>% 
    mutate(label="userLabel",
           numOtus= (file1 %>% pull(numOtus) + file2 %>% pull(numOtus))) %>% 
    select(label,numOtus,everything())

write_tsv(merge,file=outfilename)