library(patchwork)
library(cowplot)
library(ggplotify)
library(magick)
library(pdftools)

pdf("FigS5a.pdf", width = 15, height = 15, bg = "transparent")  # adjust size as needed
draw(
  cll,
  merge_legends = FALSE,
  annotation_legend_side = "right",
  legend_gap = unit(0.25, "cm")
)
dev.off()

figs5a <- image_read_pdf("FigS5a.pdf") %>%
  image_transparent("white")

pdf("FigS5b.pdf", width = 7, height = 15, bg = "transparent")  # adjust size as needed
draw(
  tall,
  merge_legends = FALSE,
  annotation_legend_side = "right",
  legend_gap = unit(0.25, "cm")
)
dev.off()

figs5b <- image_read_pdf("FigS5b.pdf") %>%
  image_transparent("white")

pdf("FigS5c.pdf", width = 7, height = 15, bg = "transparent")  # adjust size as needed
draw(
  aml,
  merge_legends = FALSE,
  annotation_legend_side = "right",
  legend_gap = unit(0.25, "cm")
)
dev.off()

figs5c <- image_read_pdf("FigS5c.pdf") %>%
  image_transparent("white")

pdf("FigS5d.pdf", width = 7, height = 15, bg = "transparent")  # adjust size as needed
draw(
  ball,
  merge_legends = TRUE,
  annotation_legend_side = "right",
  legend_gap = unit(0.25, "cm")
)
dev.off()

figs5d <- image_read_pdf("FigS5d.pdf") %>%
  image_transparent("white")

FigS5a <- ggdraw() + draw_image(figs5a, scale = 1.4)
FigS5b <- ggdraw() + draw_image(figs5b, scale = 1.4)
FigS5c <- ggdraw() + draw_image(figs5c, scale = 1.4)
FigS5d <- ggdraw() + draw_image(figs5d, scale = 1.4)


row1 <- plot_grid(FigS5a, 
                  ncol=1, 
                  labels=c("A"), label_size=14)

row2 <- plot_grid(FigS5c, FigS5d, FigS5b,
                  ncol=3, rel_widths=c(0.3, 0.3, 0.3),
                  labels=c("B", "C", "D"), label_size=14)


SFig.5 <- plot_grid(row1, row2,
                    nrow=2, rel_heights=c(0.5, 0.5)) &
  theme(plot.margin=margin(r=1,l=1,t=1,b=1,unit="cm"))


ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.5/SFig.5.pdf", plot=SFig.5, width=42, height=59.4, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/SFig.5/SFig.5.png", plot=SFig.5, width=42, height=59.4, units = "cm")
