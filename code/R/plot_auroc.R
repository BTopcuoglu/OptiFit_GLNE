library(tidyverse)
library(cowplot)

data <- read_csv(snakemake@input[["perf"]])
pvals <- read_csv(snakemake@input[["pvals"]])
#data <- read_csv("data/learning/summary/merged_performance.csv")
#pvals <- read_csv("results/tables/pvalues.csv")

#check for significance, this plot shows "NS", throw error if p < 0.05
if((pvals  %>% pull(p_value) %>% min()) < 0.05){
       stop("need to adjust significance indicator in plot")
}

means <- data %>% 
  select(cv_metric_AUC,AUC,algorithm) %>% 
  rename(Train=cv_metric_AUC,
         Test=AUC) %>% 
  pivot_longer(!algorithm,names_to="type",values_to = "AUC") %>% 
  group_by(algorithm,type) %>% 
  summarise(mean_AUC = round(mean(AUC),digits=3)) %>% 
  mutate(type = factor(type,levels=c("Train","Test")),
         algorithm = factor(algorithm,levels=c("opticlust","optifit"),labels=c("OptiClust","OptiFit")))

plot <- data %>% 
  select(cv_metric_AUC,AUC,algorithm) %>% 
  rename(Train=cv_metric_AUC,
         Test=AUC) %>% 
  pivot_longer(!algorithm,names_to="type",values_to = "AUC") %>% 
  mutate(type = factor(type,levels=c("Train","Test")),
         algorithm = factor(algorithm,levels=c("opticlust","optifit"),labels=c("OptiClust","OptiFit"))) %>% 
  ggplot(aes(x=algorithm,y=AUC,color=algorithm)) +
    geom_jitter(height=0,width=0.2,alpha = 0.6,size=2) + 
    stat_summary(fun = mean, geom="pointrange",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),
                 color="black") +
    theme_bw() +
    ylab("AUROC") +
    xlab("") +
    facet_grid(.~type) +
    theme(legend.position = "none",
          axis.text= element_text(size=14),
          axis.title = element_text(size=14),
          panel.spacing = unit(2, "lines"),
          strip.text.x = element_text(size = 14)) +
    geom_text(data=means,aes(x=algorithm,y=mean_AUC,label=mean_AUC),
              color="black",hjust=-1.1,size=4.2) +
    scale_color_manual(values=c("#48A3DC","#F79271")) +
    geom_line(data = tibble(x=c(1,2),y=c(0.85,0.85)),
              aes(x=x,y=y),
             inherit.aes = FALSE) +
    geom_line(data = tibble(x=c(1,1),y=c(0.85,0.845)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
    geom_line(data = tibble(x=c(2,2),y=c(0.85,0.845)),
            aes(x=x,y=y),
            inherit.aes = FALSE) +
    geom_text(data = tibble(x=c(1.5),y=c(.86)),
              aes(x=x,y=y),label="NS",size=3.8,
              inherit.aes = FALSE) 
    
cowplot::plot_grid(plot,NULL,labels = c("A","B"),label_x=c(0,-48),rel_widths = c(100,1))
ggsave("results/figures/avg_auroc.png",height=5,width=8,bg='#ffffff')
