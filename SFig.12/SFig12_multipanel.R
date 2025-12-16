library(patchwork)
library(cowplot)

r1c1 <- plot_grid(SFig12a, SFig12b, 
                  align = "h", 
                  axis = "tb", 
                  ncol = 2, 
                  rel_widths = c(1, 1), 
                  labels = c("A", "B"), 
                  label_size = 14)

row1<- plot_grid(r1c1, NULL, ncol = 2, rel_widths = c(0.5, 0.5))

row2 <- plot_grid(SFig12c, SFig12d, align = "h", axis = c("t", "b"), ncol = 2, rel_widths = c(0.65, 0.35), 
               labels=c("C", "D", ""), label_size=14)

SFig12 <- plot_grid(row1, row2, NULL, nrow =3, 
                    align = "v", axis = c("l"), rel_heights = c(0.2, 0.25, 0.55)) & theme(plot.margin=margin(l=1,t=1,b=1,r=1,unit="cm"))

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.12/SFig12.pdf", plot=SFig12, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.12/SFig12.png", plot=SFig12, width=42, height=59.4, units = "cm", dpi = 300)
