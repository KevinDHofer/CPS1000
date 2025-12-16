library(patchwork)
library(cowplot)
library(ggplotify)
library(magick)
library(pdftools)

pdf("FigS11c.pdf", width = 11, height = 6, bg = "transparent")  # adjust size as needed
par(mar = c(0, 0, 0, 0))
draw(
  cps,
  merge_legends = TRUE,
  annotation_legend_side = "right",
  legend_gap = unit(0.25, "cm"),
  padding = unit(c(1, 1, 1, 1), "mm")
)
dev.off()

figs11c <- image_read_pdf("FigS11c.pdf") %>%
  image_transparent("white")

pdf("FigS11d.pdf", width = 15, height = 6, bg = "transparent")  # adjust size as needed
par(mar = c(0, 0, 0, 0))
draw(
  knis,
  merge_legends = TRUE,
  annotation_legend_side = "right",
  legend_gap = unit(0.25, "cm"),
  padding = unit(c(1, 1, 1, 1), "mm")
)
dev.off()

figs11d <- image_read_pdf("FigS11d.pdf") %>%
  image_transparent("white")

figs11d <- image_read_pdf("FigS11d.pdf") %>%
  image_transparent("white")

SFig11c <- ggdraw() + draw_image(figs11c, scale = 1.5)
SFig11d <- ggdraw() + draw_image(figs11d, scale = 1.5)

r2c1 <- plot_grid(SFig11e, SFig11f, nrow = 2, rel_heights = c(0.5, 0.5), 
                  labels=c("E", "F"), label_size=14)

r2c2 <- plot_grid(r2c1, NULL, ncol =2, rel_widths = c(0.6, 0.4))

r2 <- plot_grid(SFig11c, r2c2,
                ncol=2, rel_widths = c(0.5, 0.5), 
                labels=c("C", ""), label_size=14)

r1c1 <- plot_grid(NULL,SFig11a, nrow=2, rel_heights = c(0.15, 0.85))

r1 <- plot_grid(NULL, r1c1, NULL, SFig11b, ncol = 4, rel_widths = c(0.025,0.45, 0.025, 0.5), 
                  labels=c("A", "", "", "B"), label_size=14)

r3 <- plot_grid(SFig11d, NULL,
                  ncol = 2, rel_widths = c(0.7, 0.3), 
                  labels=c("D", ""), label_size=14)

SFig11 <- plot_grid(r1, r2, r3, NULL, nrow =4, rel_heights = c(0.15, 0.25, 0.25, 0.35)) & theme(plot.margin=margin(l=1,t=1,b=1,unit="cm"))
SFig11

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.11/SFig11.pdf", plot=SFig11, width=42, height=59.4, units = "cm")

#ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.11/SFig11.png", plot=SFig11, width=42, height=59.4, units = "cm", dpi = 300)



