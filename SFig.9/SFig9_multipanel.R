library(patchwork)
library(cowplot)

row1 <- plot_grid(SFig9a, SFig9b, ncol = 2, rel_widths = c(0.5, 0.5), labels=c("A", "B"), label_size=14)

colde <- plot_grid(SFig9d, SFig9e, ncol = 1, rel_heights = c(0.6, 0.4), align = "v", axis = c("l", "r"), labels=c("D", "E"), label_size=14)

row2 <- plot_grid(SFig9c, colde, ncol = 2, rel_widths = c(0.6, 0.4), labels=c("C", ""), label_size=14)

row3 <- plot_grid(SFig9f, ncol = 1, labels=c("F"), label_size=14)

SFig.9 <- plot_grid(row1, row2, row3, NULL, rel_heights = c(0.15, 0.5, 0.25, 0.1), ncol = 1, 
                   align = "v", axis = c("l", "r")) & theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.9/SFig9.pdf", plot=SFig.9, width=42, height=59.4, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.9/SFig9.png", plot=SFig.9, width=42, height=59.4, units = "cm", dpi = 300)