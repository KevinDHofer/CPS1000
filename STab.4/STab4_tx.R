library(flextable)
library(readxl)
library(ggplot2)

tx <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.4/STab.4_drugs_tx.xlsx")

ft <- flextable(tx) |> 
  theme_alafoli() |>
  autofit()

ft <- ft |> 
  fontsize(size = 8, part = "header") |> 
  fontsize(size = 8, part = "body") |> 
  
  # Adjust padding for cleaner look
  padding(padding.top = 3, padding.bottom = 3, part = "all") |> 
  
  # Alignment: left-align text columns, right-align numbers
  align(j = 1, align = "left", part = "body") |>  
  align(j = 2:ncol(drugs), align = "left", part = "body") |>  # Numeric columns centered, not right
  align(align = "left", part = "header") |>
  
  # Header styling
  bold(part = "header") |>
  
  # Line spacing
  line_spacing(space = 1, part = "all") |> 
  
  # Make all fonts black
  color(color = "black", part = "all")

ft_grob <- gen_grob(ft, fit = "width", scaling = "min", hjust = 0)

# Wrap in a ggplot with margins
STab4 <- ggplot() +
  annotation_custom(ft_grob) +
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
STab4

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.4/STab4.pdf", plot=STab4, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.4/STab4.svg", plot=STab4, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.4/STab4.png", plot=STab4, width=21, height=29.7, units = "cm")
