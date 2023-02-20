library(tidyverse)

outdir <- snakemake@params[["outdir"]]
if(!dir.exists(outdir)){dir.create(outdir)}

files <- snakemake@input[["splits"]]

split_files <- files  %>% 
    lapply(read_csv,show_col_types=F)  %>% 
    bind_rows()  %>% 
    mutate(group=as.character(group))

counts <- split_files  %>% 
    group_by(group)  %>% 
    count(train_test)  

set.seed(1234)
plot <- counts  %>% 
    ggplot(aes(x=train_test,y=n,color=train_test)) +
        geom_jitter(size=2,alpha=0.3,width=0.4,height = 0) +
        theme_classic(base_size=18) +
        scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
        theme(panel.grid.major.y = element_line(color="grey"),
              axis.text = element_text(size=14),
              axis.title = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14)) + 
        xlab("Group") + ylab("Number of times sample is in group") +
        # stat_summary(fun.data = "mean_sdl",fun.args = list(mult = 1),
        #              geom="errorbar",color="black")
        stat_summary(fun = mean, geom = "pointrange",
                     fun.max = function(x) mean(x) + sd(x),
                     fun.min = function(x) mean(x) - sd(x),
                     color="black") 
ggsave(paste0(outdir,"split_8020.png"),width = 6, height= 6)
