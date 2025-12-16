library(patchwork)
library(cowplot)

row3 <- plot_grid(SFig8c, NULL, ncol = 2, rel_widths = c(0.4, 0.6))

SFig.8 <- plot_grid(SFig8a, SFig8b, row3, rel_heights = c(0.2, 0.65, 0.15), ncol = 1, 
                   align = "v", axis = c("l", "r"), 
                   labels=c("A", "B", "C"), label_size=14) & theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.8/SFig8.pdf", plot=SFig.8, width=42, height=59.4, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.8/SFig8.png", plot=SFig.8, width=42, height=59.4, units = "cm", dpi = 300)