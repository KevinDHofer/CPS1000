library(patchwork)
library(cowplot)

library(png)

col1 <- plot_grid(Fig6a, Fig6b, 
                  align = "v", axis = c("r"), 
                  ncol = 1, rel_heights = c(0.4, 0.4), 
                  labels=c("A", "B"), label_size=14)

col2 <- plot_grid(Fig6c, ncol = 1, 
                  #rel_heights = c(0.8, 0.2),
                  labels=c("C", ""), label_size=14)

col3 <- plot_grid(Fig6d,Fig6e, ncol = 1, align = "v", axis = c("l", "r"),
                  rel_heights = c(0.5, 0.5), 
                  labels=c("D", "E"), label_size=14)

top <- plot_grid(col1, col2, col3, ncol=3, 
                 align = "h", axis = c("t", "b"),
                 rel_widths = c(0.2, 0.4, 0.4))

Fig6 <- plot_grid(top, NULL, nrow =2, rel_heights = c(0.25, 0.75)) & theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.6/Fig6.pdf", plot=Fig6, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.6/Fig6.png", plot=Fig6, width=42, height=59.4, units = "cm", dpi = 300)
