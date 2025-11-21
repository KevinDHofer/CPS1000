library(tidyverse)
library(dplyr)
library(plotrix)
library(drc)
library(scales)

#load data
usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

#load table with disease entities  
diseases <- usedSamples_all |> 
  arrange(patientID) |> 
  dplyr::select(patientID, diagnosis) 

# merge the tables
merged_table <- cbind(diseases, viability_table) |> 
  relocate(diagnosis, .before = X10058.F4_1)
merged_table

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:318), names_to = "drug", values_to = "viability")
merged_table_long

#select one disease: CLL
merged_filtered <- merged_table_long |> filter(str_detect(diagnosis, "CLL"))

#remove last two characters from each string in 'drug' column
merged_filtered<- merged_filtered |> 
  mutate(drug_sub = str_sub(merged_filtered$drug, end = -3)) |> 
  filter(diagnosis != "")

#select venetoclax and diseases
merged_filtered <- merged_filtered |> 
  filter(drug_sub == "Venetoclax") |> 
  filter(diagnosis == "CLL") |> 
  mutate(drug = fct_recode(drug, "0.00064" = "Venetoclax_5",
                           "0.0032" = "Venetoclax_4", "0.016" = "Venetoclax_3", "0.08" = "Venetoclax_2", "0.4" = "Venetoclax_1"))

#make concentrations numeric
merged_filtered$drug <- as.numeric(as.character(merged_filtered$drug))

#calculate medians for each disease
merged_median <- merged_filtered |> 
  group_by(diagnosis, drug) |> summarize(med=median(viability))

merged_median_CLL <- merged_median |> 
  filter(diagnosis == "CLL")

# Fit the dose-response model, fixed at 1

#CLL
model_CLL <- drm(med ~ drug, data = merged_median_CLL, 
                 fct = LL.4(fixed = c(NA, 1, NA, NA), 
                            names = c("hill", "min_value", "Einf", "ec_50")))

# Create a new data frame for prediction

#CLL
new_data_CLL <- data.frame(drug = seq(min(merged_median_CLL$drug), max(merged_median_CLL$drug), length.out = 500))
new_data_CLL$pred <- predict(model_CLL, newdata = new_data_CLL)

# Add diagnosis column to each prediction dataframe
new_data_CLL$diagnosis <- "CLL"

all_predictions <- new_data_CLL

#add ic50

#CLL
coefs_CLL <- setNames(c(model_CLL$coefficients, 1), c("hill", "Einf", "ec_50", "min_value"))
#ic_50_CLL <- with(as.list(coefs_CLL),exp(log(ec_50) + (1 / hill) * log(min_value / (min_value - 2 * min_value))))
ic_50_CLL <- with(as.list(coefs_CLL), exp(log(ec_50) + (1 / hill) * log(min_value / (min_value - 2 * 0))))


#set colors
#mycolors <- setNames(c("#1F78B4", "#B2DF8A", "#A6CEE3", "darkgrey", "#E31A1C", "#FB9A99", "#FF7F00", "#CAB2D6"), c("MCL", "CLL", "T-PLL", "Sezary", "AML", "T-ALL", "B-ALL", "B-PLL"))  


# Update the geom_line to use diag instead of diagnosis
geom_line(data = all_predictions, aes(x = drug, y = pred, color = diag))

# Do the same for all_ic50
#all_ic50$diag <- factor(all_ic50$diagnosis, levels = c("CLL", "B-ALL", "MCL", "AML", "T-ALL", "Sezary"))

#plot
Fig7c <- merged_filtered |> 
  group_by(drug) |> 
  mutate(viability = median(viability)) |> 
  ggplot(aes(x=drug, y=viability)) +
  geom_point(color="black", fill="white", alpha=0.8)+
  geom_line(data = all_predictions, aes(x = drug, y = pred, color = diagnosis))+
  geom_ribbon(data = all_predictions, aes(x = drug, ymin = 0, ymax = pred, fill = diagnosis), alpha = 0.3, inherit.aes = FALSE) +
  labs(title="", x=bquote("Concentration ("*mu*"M)"), y = "Viability") +
  theme_classic() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=5), labels = trans_format("log10", math_format(10^.x)), 
                limits = c(0.00064, 1))+
  scale_y_continuous(breaks=seq(0,1.0, 0.2), limits=c(0,1.1)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0, color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), 
        axis.title.y = element_text(size=8), 
        axis.title.x = element_text(size=8))+
  geom_hline(yintercept=coefs_CLL[[2]], linetype="dashed", color="grey80")+
  annotate("text", x=0.75, y=coefs_CLL[[2]]+0.05, label="E[inf]", size=3, parse=TRUE)+
  geom_vline(xintercept = ic_50_CLL, linetype = "dashed", color = "grey80")+
  annotate("text", x=ic_50_CLL+0.008, y=0.75, label="IC[50]", size=3, parse=TRUE)+
  annotate("text", x=0.0025, y=0.6, label="AUC", color="black", size=3)+
  annotate("text", x=0.01, y=1, label="R^2", size=3, parse=TRUE)+
  geom_segment(aes(x = ic_50_CLL-0.0025, y = 0.47, xend = ic_50_CLL+0.00325, yend = 0.47), color="black")+
  geom_segment(aes(x = ic_50_CLL-0.0025, y = 0.47, xend = ic_50_CLL-0.0025, yend = 0.60), color="black")+
  geom_segment(aes(x = ic_50_CLL-0.0025, y = 0.60, xend = ic_50_CLL+0.00325, yend = 0.47), color="black", linetype = "dotted")+
  annotate("text", x=ic_50_CLL+0.06, y=0.5, label="Hill coefficient\n(slope)", size=3)
Fig7c <- Fig7c+theme(plot.margin = margin(t=1, b=1, r=1, unit = "cm"))
Fig7c