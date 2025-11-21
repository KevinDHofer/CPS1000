library(patchwork)
library(cowplot)

col1row3 <- plot_grid(Fig7c, Fig7d, NULL,
                            align = "h", axis = c("t", "b"), ncol = 3, rel_widths = c(0.4, 0.5, 0.1), 
                            labels=c("", "D", ""), label_size=14)

col1 <- plot_grid(Fig7a, Fig7b, col1row3, NULL,
                  align = "v", axis = c("l", "r"), ncol = 1, rel_heights = c(0.2, 0.2, 0.15, 0.45), 
                  labels=c("A", "B", "C"), label_size=14)

col2 <- plot_grid(NULL, Fig7e, Fig7f, NULL,
                  align = "v", axis = c("l", "r"), ncol = 1, rel_heights = c(0.025, 0.1, 0.35, 0.575), 
                  labels=c("E", "", "F", ""), label_size=14)

Fig7 <- plot_grid(col1, col2, ncol =2, rel_widths = c(0.65, 0.35)) & theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))
                                    
ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.7/Fig7.pdf", plot=Fig7, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.7/Fig7.png", plot=Fig7, width=42, height=59.4, units = "cm", dpi = 300)
