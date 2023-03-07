library(tidyverse)
library(cowplot)
library(glue)
library(ggtext)

data <- read_csv(snakemake@input[["perf"]]) 
pvals <- read_csv(snakemake@input[["pvals"]])
output <- snakemake@output[["fig_auc"]]
colors <- snakemake@params[["colors"]]
order <- snakemake@params[["order"]]

names(colors) <- order
data <- data %>% 
  mutate(algorithm=case_when(algorithm == "opticlust_denovo" ~ "OptiClust de novo",
                             algorithm == "optifit_gg" ~ "OptiFit Greengenes",
                             algorithm == "vsearch_denovo" ~ "VSEARCH de novo",
                             algorithm == "vsearch_gg" ~ "VSEARCH Greengenes",
                             algorithm == "optifit_self" ~ "OptiFit Self",
                             TRUE ~ NA)) %>% 
  mutate(algorithm = factor(algorithm,levels=rev(order)))

means <- data %>% 
  select(cv_metric_AUC,AUC,algorithm) %>% 
  rename(Train=cv_metric_AUC,
         Test=AUC) %>% 
  pivot_longer(!algorithm,names_to="type",values_to = "AUC") %>% 
  group_by(algorithm,type) %>% 
  summarise(mean_AUC = format(round(mean(AUC),digits=3), nsmall=3))  %>% 
  mutate(type = factor(type,levels=c("Train","Test")))

set.seed(1234)
plot <- data %>% 
  select(cv_metric_AUC,AUC,algorithm) %>% 
  rename(Train=cv_metric_AUC,
         Test=AUC) %>% 
  pivot_longer(!algorithm,names_to="type",values_to = "AUC") %>% 
  mutate(type = factor(type,levels=c("Train","Test"))) %>% 
  ggplot(aes(y=algorithm,x=AUC,color=algorithm)) +
    geom_jitter(height=0.2,width=0,alpha = 0.4,size=2) + 
    stat_summary(fun = mean, geom="pointrange",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),
                 color="black",size=1,linewidth=1) +
    theme_bw(base_size=18) +
    xlab("AUROC") +
    ylab("") +
    facet_grid(.~type) +
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"),
          #axis.text.y = element_markdown(margin = margin(r=15))) +
          axis.text.y = element_markdown(hjust=1,margin=margin(r=15))) +
    geom_text(data=means,aes(y=algorithm,x=as.numeric(mean_AUC),label=mean_AUC),
              vjust=4,size=4.2,color="black") +
    scale_color_manual(values=colors) +
    scale_y_discrete(breaks=rev(order),
                     labels=c(glue("VSEARCH<br>*de novo*"),
                              glue("VSEARCH<br>Greengenes"),
                              glue("OptiFit<br>Greengenes"),
                              glue("OptiClust<br>*de novo*"),
                              glue("OptiFit<br>Self")))

cowplot::plot_grid(plot,NULL,labels = c("A","B"),label_x=c(0.15,-43),
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
#                             labels = rev(c("OptiClust de novo","OptiFit Self","OptiFit Greengenes","VSEARCH de novo", "VSEARCH Greengenes")))) %>% 
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
    

