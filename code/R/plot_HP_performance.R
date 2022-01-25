library(tidyverse)

outDir <- "results/figures/"
if(!dir.exists(outDir)){dir.create(outDir)}

mergedHP <- read_csv(snakemake@input[["hp"]])
#mergedHP <- read_csv("data/learning/summary/merged_HP.csv")

mergedHP  %>% 
  group_by(mtry,algorithm)  %>% 
  summarise("mean_AUC" = mean(AUC),
            "sd_AUC" = sd(AUC),
            ymin_metric = mean_AUC - sd_AUC,
            ymax_metric = mean_AUC + sd_AUC)  %>% 
  ggplot(aes(x=mtry,y=mean_AUC,color=algorithm)) +  
    geom_point(size=2,alpha=0.5) + geom_line(alpha=0.5) +
    theme_bw() +
    #facet_wrap(.~algorithm,nrow=2) +
    geom_errorbar(aes(
      ymin = .data$ymin_metric,
      ymax = .data$ymax_metric),
      width = .001) +
    ylab("Mean AUC") +
    ggtitle("Hyperparameter Performance")
ggsave(paste0(outDir,"hp_performance.png"))
