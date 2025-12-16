library(patchwork)
library(cowplot)

ttnt <- plot_grid(SFig13e, SFig13f, align = "h", axis = c("t", "b"), ncol = 2,
                  labels=c("E", "F"), label_size=14)

row1 <- plot_grid(SFig13a, SFig13c, align = "h", axis = c("t", "b"), nrow = 1, rel_widths = c(0.625,0.375), 
                  labels=c("A", "C"), label_size=14)

row2 <- plot_grid(SFig13b, SFig13d, align = "h", axis = c("t", "b"), nrow = 1, rel_widths = c(0.625,0.375), 
                      labels=c("B", "D"), label_size=14)

row3 <- plot_grid(SFig13e, SFig13f, SFig13g, align = "h", axis = c("t", "b"), nrow = 1, rel_widths = c(0.2,0.2,0.6), 
                  labels=c("E", "F", "G"), label_size=14)

row4 <- plot_grid(SFig13h, SFig13i, SFig13j, align = "h", axis = c("t"), nrow = 1, rel_widths = c(0.15,0.25,0.4), 
                  labels=c("H", "I", "J"), label_size=14)

row5 <- plot_grid(SFig13k, SFig13l, NULL, align = "h", axis = c("t", "b"), nrow = 1, rel_widths = c(0.3,0.3,0.4), 
                  labels=c("K", "L", ""), label_size=14)

SFig13 <- plot_grid(row1, row2, row3, row4, row5, NULL, nrow=6, rel_heights = c(0.2, 0.2, 0.2, 0.15, 0.15, 0.10)) & theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))
                                    
ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.13/SFig13.pdf", plot=SFig13, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.13/SFig13.png", plot=SFig13, width=42, height=59.4, units = "cm", dpi = 300)
