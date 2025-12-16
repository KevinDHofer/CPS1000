library(flextable)
library(readxl)
library(ggplot2)
library(dplyr)

aml <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.6/aml_clin_trials.xlsx")

aml$`ORR (%)`<- round(aml$`ORR (%)`,2)

aml <- aml |> 
  mutate(`CRR (%)` = ifelse(is.na(`CRR (%)`), NA, round(as.numeric(`CRR (%)`), 2)),
         `PRR (%)` = ifelse(is.na(`PRR (%)`), NA, round(as.numeric(`PRR (%)`), 2))) |> 
  mutate(across(everything(), as.character))

ft <- flextable(aml) |> 
  theme_alafoli() |>
  autofit()

ft <- ft |> 
  fontsize(size = 8, part = "header") |> 
  fontsize(size = 8, part = "body") |> 
  
  # Adjust padding for cleaner look
  padding(padding.top = 3, padding.bottom = 3, part = "all") |> 
  
  # Alignment: left-align text columns, right-align numbers
  align(align = "left") |>
  
  # Header styling
  bold(part = "header") |>
  
  # Line spacing
  line_spacing(space = 1, part = "all") |> 
  
  # Make all fonts black
  color(color = "black", part = "all")

ft_grob <- gen_grob(ft, fit = "width", scaling = "min", hjust = 0)

# Wrap in a ggplot with margins
STab6 <- ggplot() +
  annotation_custom(ft_grob) +
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
STab6

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.6/STab6.pdf", plot=STab6, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.6/STab6.svg", plot=STab6, width=21, height=29.7, units = "cm")

ggsave(file="~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/STab.6/STab6.png", plot=STab6, width=21, height=29.7, units = "cm")
