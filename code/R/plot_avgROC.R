library(tidyverse)
library(cowplot)

data <- read_csv(snakemake@input[["senspec"]])
#data <- read_csv("data/learning/summary/all_sens_spec.csv")

avg_sens <- data %>%
  group_by(specificity,algorithm) %>%
  summarise(mean_sensitivity = mean(sensitivity),
            sd_sensitivity = sd(sensitivity),
            .groups = "drop") %>%
  mutate(upper_sens = mean_sensitivity + sd_sensitivity,
         lower_sens = mean_sensitivity - sd_sensitivity) %>%
  mutate(upper_sens = case_when(upper_sens > 1 ~ 1,
                                TRUE ~ upper_sens),
         lower_sens = case_when(upper_sens < 0 ~ 0,
                                TRUE ~ lower_sens)) %>%
  mutate(fpr = 1-specificity)

plot <- avg_sens %>%
  mutate(algorithm = factor(algorithm,levels=c("opticlust","optifit"),labels=c("OptiClust","OptiFit"))) %>% 
  ggplot(aes(x=fpr,y=mean_sensitivity,
             ymin=lower_sens,ymax=upper_sens)) +
  geom_line(aes(color=algorithm),lwd=1.25,alpha=0.7) +
  geom_ribbon(aes(fill=algorithm),alpha=0.15) +
  coord_equal() +
  geom_abline(intercept = 0,lty="dashed",color="grey50") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("False Positive Rate") +
  ylab("Average True Positive Rate") +
  theme_classic() +
  theme(legend.position = c(0.28,0.98),
        legend.box = "vertical",
        axis.text.x = element_text(vjust = -0.5,size=14),
        axis.title.x = element_text(vjust= -0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.text=element_text(size=12)) +
  guides(color = guide_legend(nrow = 1,title=""),
         fill = guide_legend(nrow = 1,title="")) +
  scale_color_manual(values=c("#48A3DC","#F79271")) +
  scale_fill_manual(values=c("#48A3DC","#F79271")) 

plot_grid(plot,labels=c("C"),label_x=c(-0.01))
ggsave("results/figures/avg_roc.png",height=5,width=5,bg='#ffffff')
