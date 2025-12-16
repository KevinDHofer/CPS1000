library(flextable)
library(readxl)
library(ggplot2)


drugs <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.3/stable3_drugs.xlsx") |> 
  dplyr::select(-Location)

drugs$Drug <- sort(drugs$Drug)

ft <- flextable(drugs) |> 
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
  
  # Set consistent width if needed
  width(j=1, width = 1.25) |>  # Adjust as needed
  width(j=2, width = 3.5) |>
  width(j=3, width = 0.75) |>
  
  # Line spacing
  line_spacing(space = 1, part = "all") |> 
  
# Make all fonts black
color(color = "black", part = "all")

ft_grob <- gen_grob(ft, fit = "width", scaling = "min", hjust = 0)

# Wrap in a ggplot with margins
STab3 <- ggplot() +
  annotation_custom(ft_grob) +
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
STab3

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.3/STab3.pdf", plot=STab3, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.3/STab3.svg", plot=STab3, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.3/STab3.png", plot=STab3, width=21, height=29.7, units = "cm")
