library(flextable)
library(readxl)
library(ggplot2)
library(officer)

stable1 <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.1/stable1.xlsx")

ft <- flextable(stable1) |> 
  theme_alafoli() |>
  autofit()

ft <- ft |> 
  fontsize(size = 8, part = "header") |> 
  fontsize(size = 8, part = "body") |> 
  
  # Adjust padding for cleaner look
  padding(padding.top = 3, padding.bottom = 3, part = "all") |> 
  padding(padding.left = 3, padding.right = 3, part = "all") |>
  
  # Alignment: left-align text columns, right-align numbers
  align(j = 1, align = "center", part = "body") |>  # First column centered as well
  align(j = 2:ncol(stable1), align = "center", part = "body") |>  # Numeric columns centered, not right
  align(align = "center", part = "header") |>
  
  # Header styling
  bold(part = "header") |>
  
  # Set consistent width if needed
  width(j=1, width = 0.5) |>  # Adjust as needed
  width(j=2:ncol(stable1), width = 1) |>  # Adjust as needed
  
  # Line spacing
  line_spacing(space = 1.15, part = "all") |>
  
  # Add header with Greek letter mu (μ)
  vline(j = "Group") |>
  add_header_row(
    values = c("", "Concentration index (\u03BCM)"),  # \u03BC is the Unicode for μ
    colwidths = c(1, 5)
  ) |>
  
  # Make all fonts black
  color(color = "black", part = "all")

ft_grob <- gen_grob(ft, fit = "width", scaling = "min", hjust = 0)

# Wrap in a ggplot with margins
STab1 <- ggplot() +
  annotation_custom(ft_grob) +
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
STab1

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.1/STab1.pdf", plot=STab1, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.1/STab1.svg", plot=STab1, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.1/STab1.png", plot=STab1, width=21, height=29.7, units = "cm")
