library(tidyverse)

outDir <- snakemake@params[["outdir"]]
if(!dir.exists(outDir)){dir.create(outDir)}

mergedHP <- read_csv(snakemake@input[["hp"]]) 
colors <- snakemake@params[["colors"]]
order <- snakemake@params[["order"]]
names(colors) <- order

mergedHP <- mergedHP %>% 
  mutate(algorithm=case_when(algorithm == "opticlust_denovo" ~ "OptiClust de novo",
                             algorithm == "optifit_gg" ~ "OptiFit Greengenes",
                             algorithm == "vsearch_denovo" ~ "VSEARCH de novo",
                             algorithm == "vsearch_gg" ~ "VSEARCH Greengenes",
                             algorithm == "optifit_self" ~ "OptiFit Self",
                             TRUE ~ NA)) %>% 
  mutate(algorithm = factor(algorithm,levels=order)) 
  
plot <- mergedHP  %>% 
  group_by(mtry,algorithm)  %>% 
  summarise("mean_AUC" = mean(AUC),
            "sd_AUC" = sd(AUC),
            ymin_metric = mean_AUC - sd_AUC,
            ymax_metric = mean_AUC + sd_AUC)  %>% 
  ggplot(aes(x=mtry,y=mean_AUC,color=algorithm)) +  
    geom_point(size=4,alpha=0.8) + geom_line(alpha=0.5) +
    theme_bw(base_size=18) +
    #facet_wrap(.~algorithm,nrow=2) +
    geom_errorbar(aes(
      ymin = .data$ymin_metric,
      ymax = .data$ymax_metric),
      width = .001) +
    ylab("Mean AUC") +
    ggtitle("Hyperparameter Performance") +
    scale_color_manual(values=colors,name="Algorithm")
ggsave(paste0(outDir,"hp_performance.png"),width=10)


plot2 <- mergedHP  %>% 
  group_by(mtry,algorithm)  %>% 
  summarise("mean_AUC" = mean(AUC),
            "sd_AUC" = sd(AUC),
            ymin_metric = mean_AUC - sd_AUC,
            ymax_metric = mean_AUC + sd_AUC)  %>% 
  ggplot(aes(x=mtry,y=mean_AUC,color=algorithm)) +  
    geom_point(size=4,alpha=0.8) + geom_line(alpha=0.5) +
    theme_bw(base_size=18) +
    facet_wrap(.~algorithm,nrow=1) +
    geom_errorbar(aes(
      ymin = .data$ymin_metric,
      ymax = .data$ymax_metric),
      width = .001) +
    ylab("Mean AUC") +
    ggtitle("Hyperparameter Performance") +
    scale_color_manual(values=colors) +
    theme(legend.position="none")
ggsave(paste0(outDir,"hp_performance_grid.png"),width=13)
