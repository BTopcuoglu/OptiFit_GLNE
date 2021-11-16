library(tidyverse)

outdir <- "analysis/"
if(!dir.exists(outdir)){dir.create(outdir)}

files <- commandArgs(trailingOnly = TRUE)

split_files <- files  %>% 
    lapply(read_csv)  %>% 
    bind_rows()  %>% 
    mutate(group=as.character(group))

counts <- split_files  %>% 
    group_by(group)  %>% 
    count(train_test)  

counts  %>% 
    ggplot(aes(x=train_test,y=n,color=train_test)) +
        geom_jitter(size=2,alpha=0.5,width=0.4,height = 0) +
        theme_classic() +
        scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
        theme(panel.grid.major.y = element_line(color="grey"),
              axis.text = element_text(size=14),
              axis.title = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14)) + 
        xlab("Group") + ylab("Number of times sample is in group") 
ggsave("analysis/view_splits.png",width = 6, height= 6)
