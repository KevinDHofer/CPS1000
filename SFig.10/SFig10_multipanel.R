library(patchwork)
library(cowplot)

r1c1 <- plot_grid(SFig10a, SFig10b, align = "h", axis = c("t", "b"), ncol = 2, rel_widths = c(0.4, 0.6), 
                  labels=c("A", "B"), label_size=14)

r2c1 <- plot_grid(SFig10d, SFig10e,
                  align = "v", axis = c("l", "r"), nrow = 2, rel_heights = c(0.5, 0.5), 
                  labels=c("D", "E"), label_size=14)

r3c1 <- plot_grid(SFig10g, SFig10h,
                      align = "h", axis = c("t", "b"), nrow = 1, rel_widths = c(0.6, 0.4), 
                      labels=c("G", "H"), label_size=14)

col1 <- plot_grid(r1c1, r2c1, r3c1, SFig10i, nrow = 4, rel_heights = c(0.1, 0.4, 0.15, 0.25),
                  labels=c("", "", "", "I"), label_size=14)

SFig10c <- (SFig10c1/SFig10c2 + plot_layout(ncol=1, guides = 'collect') & theme(legend.position='right'))

col2 <- plot_grid(SFig10c, SFig10f, NULL,
                  align = "v", axis = c("l", "r"), ncol = 1, rel_heights = c(0.5, 0.25, 0.25), 
                  labels=c("C", "F",""), label_size=14)

SFig10 <- plot_grid(col1, col2, ncol =2, rel_widths = c(0.6, 0.4)) & theme(plot.margin=margin(l=1,t=1,b=1,unit="cm"))

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.10/SFig10.pdf", plot=SFig10, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.10/SFig10.png", plot=SFig10, width=42, height=59.4, units = "cm", dpi = 300)
