library(tidyverse)
library(cowplot)
library(png)

p1 <- readPNG(snakemake@input[["fig2ab"]])
p2 <- readPNG(snakemake@input[["fig2c"]])

ab <- ggplot2::ggplot() + ggplot2::annotation_custom(grid::rasterGrob(p1,
                                                width=ggplot2::unit(1,"npc"),
                                                height=ggplot2::unit(1,"npc")),
                               -Inf, Inf, -Inf, Inf)

c <- ggplot2::ggplot() + ggplot2::annotation_custom(grid::rasterGrob(p2,
                                                width=ggplot2::unit(0.95,"npc"),
                                                height=ggplot2::unit(0.95,"npc")),
                               -Inf, Inf, -Inf, Inf) +
    theme(panel.background = element_rect(fill = '#ffffff', colour = '#ffffff'))

plot_grid(ab,c,rel_widths=c(2,1.2))
#ggsave("submission/figures/fig2.tiff",width=7,height=3,units="in",bg='#ffffff')
ggsave("submission/figures/fig2.png",width=7,height=3,units="in",bg='#ffffff')