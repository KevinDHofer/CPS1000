library(flextable)
library(readxl)
library(cowplot)
library(patchwork)
library(gtable)
library(grid)


data <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/raw data for figures/DDX3X_mutations_20241017.xlsx", 
                                       sheet = "table")
data <- data[,-9]
data$`Allele frequency` <- round(data$`Allele frequency`,2)

ft <- flextable(data) |> 
  theme_alafoli()
  
ft <-  ft |> 
  fontsize(size = 8, part = "header") |> 
  fontsize(size = 8, part = "body") |> 
  
  # Adjust padding for cleaner look
  padding(padding.top = 3, padding.bottom = 3, part = "all") |> 
  
  # Alignment: left-align text columns, right-align numbers
  align(align = "center", part = "all") |> 
  
  # Header styling
  bold(part = "header") |>
  
  # Set consistent width if needed
  width(j=c(1,2,3,6,8), width = 0.1) |>  # Adjust as needed
  width(j=c(4,7), width = 0.25) |>  # Adjust as needed
  width(j=c(5), width = 0.5) |>  # Adjust as needed
  
  # Line spacing
  line_spacing(space = 1.15, part = "all") |>
  
  # Add header with Greek letter mu (μ)
  vline(j = "ID") |>
  
  # Make all fonts black
  color(color = "black", part = "all")

SFig11a <- gen_grob(ft, fit = "width", scaling = "min", hjust = 0)

