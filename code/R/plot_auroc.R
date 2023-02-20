library(tidyverse)
library(cowplot)

data <- read_csv(snakemake@input[["perf"]])
pvals <- read_csv(snakemake@input[["pvals"]])
output <- snakemake@output[["fig_auc"]]
colors <- snakemake@params[["colors"]]
#data <- read_csv("results/ml/summary/merged_performance.csv")
#pvals <- read_csv("results/tables/pvalues.csv")
#colors <- c("#20639b","#3caea3","#f5ad5b","#ed553b","#a989ba")

means <- data %>% 
  select(cv_metric_AUC,AUC,algorithm) %>% 
  rename(Train=cv_metric_AUC,
         Test=AUC) %>% 
  pivot_longer(!algorithm,names_to="type",values_to = "AUC") %>% 
  group_by(algorithm,type) %>% 
  summarise(mean_AUC = round(mean(AUC),digits=3)) %>% 
  mutate(type = factor(type,levels=c("Train","Test")),
         algorithm = factor(algorithm,levels=c("opticlust_denovo","optifit_self","optifit_gg","vsearch_denovo","vsearch_gg"),
                            labels = c("OptiClust de novo","OptiFit Self","OptiFit GreenGenes","VSEARCH de novo", "VSEARCH GreenGenes")))

set.seed(1234)
plot <- data %>% 
  select(cv_metric_AUC,AUC,algorithm) %>% 
  rename(Train=cv_metric_AUC,
         Test=AUC) %>% 
  pivot_longer(!algorithm,names_to="type",values_to = "AUC") %>% 
  mutate(type = factor(type,levels=c("Train","Test"),
                       labels=),
         algorithm = factor(algorithm,levels=rev(c("opticlust_denovo","optifit_self","optifit_gg","vsearch_denovo","vsearch_gg")),
                            labels = rev(c("OptiClust de novo","OptiFit Self","OptiFit GreenGenes","VSEARCH de novo", "VSEARCH GreenGenes")))) %>% 
  ggplot(aes(y=algorithm,x=AUC,color=algorithm)) +
    geom_jitter(height=0.2,width=0,alpha = 0.4,size=2) + 
    stat_summary(fun = mean, geom="pointrange",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),
                 color="black") +
    theme_bw(base_size=18) +
    xlab("AUROC") +
    ylab("") +
    facet_grid(.~type) +
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines")) +
    geom_text(data=means,aes(y=algorithm,x=mean_AUC,label=mean_AUC),
              vjust=4,size=4.2,color="black") +
    scale_color_manual(values=colors) 

cowplot::plot_grid(plot,NULL,labels = c("A","B"),label_x=c(0.24,-38),
                   label_y = c(0.98),rel_widths = c(100,1),label_size=14)
ggsave(output,height=5,width=8,scale=1.2,bg='#ffffff')

#median and IQR  
# data %>% 
#   select(cv_metric_AUC,AUC,algorithm) %>% 
#   rename(Train=cv_metric_AUC,
#          Test=AUC) %>% 
#   pivot_longer(!algorithm,names_to="type",values_to = "AUC") %>% 
#   mutate(type = factor(type,levels=c("Train","Test"),
#                        labels=),
#          algorithm = factor(algorithm,levels=rev(c("opticlust_denovo","optifit_self","optifit_gg","vsearch_denovo","vsearch_gg")),
#                             labels = rev(c("OptiClust de novo","OptiFit Self","OptiFit GreenGenes","VSEARCH de novo", "VSEARCH GreenGenes")))) %>% 
#   ggplot(aes(y=algorithm,x=AUC,color=algorithm)) +
#     geom_jitter(height=0.2,width=0,alpha = 0.4,size=2) + 
#     stat_summary(fun.data = median_hilow,fun.args=(conf.int=0.5), show.legend=FALSE,
#                  geom="crossbar",color="black",lwd=0.4,width=0.6) +
#     theme_bw(base_size=18) +
#     xlab("AUROC") +
#     ylab("") +
#     facet_grid(.~type) +
#     theme(legend.position = "none",
#           panel.spacing = unit(2, "lines")) +
# #     geom_text(data=means,aes(y=algorithm,x=mean_AUC,label=mean_AUC),
# #               vjust=4,size=4.2,color="black") +
#     scale_color_manual(values=colors)  
    
