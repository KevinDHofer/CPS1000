library(patchwork)
library(cowplot)

Fig5b12 <- (Fig5b/Fig5b2 + plot_layout(ncol=1, guides = 'collect') & theme(legend.position='right'))

col1 <- plot_grid(Fig5a, Fig5c, NULL, Fig5e, NULL,
                  align = "v", axis = c("l", "r"), ncol = 1, rel_heights = c(0.2, 0.12, 0.25, 0.18, 0.05), 
                   labels=c("A", "C", "D", "E", ""
                            #, "F"
                            ), label_size=14)

col2 <- plot_grid(Fig5b12, Fig5f, Fig5g, NULL, 
                  align = "v", axis = c("r"),
                  ncol = 1, rel_heights = c(0.55, 0.22, 0.28, 0.05), 
                  labels=c("B", "F", "G", ""), label_size=14)

Fig5 <- plot_grid(col1, col2, ncol=2, rel_widths = c(0.6, 0.4)) & theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.5/Fig5.pdf", plot=Fig5, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.5/Fig5.png", plot=Fig5, width=42, height=59.4, units = "cm", dpi = 300)