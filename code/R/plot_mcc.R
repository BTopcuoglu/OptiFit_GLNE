library(tidyverse)
library(ggplot2)

mcc_file <- snakemake@input[["mcc"]]
colors <- snakemake@params[["colors"]]
outfile <- snakemake@output[["plot"]]

# mcc_file <- "results/ml/summary/merged_mcc.csv"
# colors <- c("#20639b","#3caea3","#f5ad5b","#ed553b","#a989ba")
# outfile <- "results/figures/mcc_scores.png"

mcc <- read_csv(mcc_file) %>% 
    mutate(Group=case_when(algorithm == "opticlust_denovo" ~ "OptiClust de novo",
                           algorithm == "optifit_gg" ~ "OptiFit GreenGenes",
                           algorithm == "vsearch_denovo" ~ "VSEARCH de novo",
                           algorithm == "vsearch_gg" ~ "VSEARCH GreenGenes",
                           subset == "train" ~ "OptiFit Self - Reference",
                           state == "fit" ~ "OptiFit Self - Fit",
                           TRUE ~ "OptiFit Self"))   %>% 
  filter(!(Group %in% c("OptiFit Self - Reference","OptiFit Self - Fit"))) %>%
  mutate(Group = factor(Group,
                        levels=c("OptiClust de novo","OptiFit Self","OptiFit GreenGenes",
                                 "VSEARCH de novo","VSEARCH GreenGenes")))
                        # levels=c("VSEARCH GreenGenes","VSEARCH de novo","OptiFit GreenGenes",
                        #          "OptiFit Self","OptiClust de novo")))

plot <- mcc  %>% 
  ggplot(aes(y=Group,x=mcc,fill=Group)) +
    geom_jitter(height=0.2,width=0,alpha=0.7,size=5,shape=21,color="black") +
    theme_bw(base_size=18) +
    ylab("") + xlab("MCC") +
    scale_x_continuous(breaks=seq(0.6,0.9,0.1),
                       minor_breaks=seq(0.55,0.85,0.1)) +
    theme(panel.grid.major.x = element_line(color="grey85"),
          panel.grid.minor.x = element_line(color="grey92"),
          panel.grid.major.y = element_blank(),
          text = element_text(size = 13)) +
    scale_fill_manual(values=colors)

ggsave(outfile,height=7,width=10)
