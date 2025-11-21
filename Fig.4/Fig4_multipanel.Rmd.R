library(patchwork)
library(cowplot)

Fig4 <- plot_grid(Fig4a,Fig4b, NULL, nrow = 3, rel_heights = c(0.45, 0.1,0.45), 
                  labels=c("A", "B", ""), label_size=14) & theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))
Fig4

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.4/Fig4.pdf", plot=Fig4, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.4/Fig4.png", plot=Fig4, width=42, height=59.4, units = "cm", dpi = 300)