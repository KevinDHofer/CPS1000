library(tidyverse)
library(dplyr)
library("ComplexHeatmap")
library(circlize)
library(ggpubr)
library(readxl)

#load raw data
usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

drugCategory_modified <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/drugCategory_modified.csv", sep=";")

CoreDataSet <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/heatmap_CoreDataSet_External_Summary_jlu.xlsx") |> 
  dplyr::select(-SampleID) |> 
  filter(Diagnosis == "AML")

#load table with disease entities  
diseases <- usedSamples_all |> 
  arrange(patientID) |> 
  dplyr::select (patientID, diagnosis) 
diseases

#merge the tables
merged_table <- cbind(diseases, viability_table) |> 
  relocate(diagnosis, .before = X10058.F4_1)
merged_table

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:318), names_to = "drug", values_to = "viability")
merged_table_long

#remove last two characters from each string in 'drug' column
merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
merged_table_long

#create "other" group for less prevalent diseases
merged_other <- merged_table_long |>
  filter(diagnosis == "AML")
merged_other  

#create means across concentration levels
df_means <- merged_other |> 
  group_by(patientID, diagnosis, drug) |> 
  summarize(mean(viability))
df_means

#create means across drugs for pathways
df_means <- df_means |> dplyr::rename(viab = "mean(viability)")

#pivot wider
wide <- pivot_wider(df_means, names_from = patientID, values_from = viab)
wide

#add diagnosis
diag <- usedSamples_all[, c("patientID", "diagnosis")] |> 
  filter(diagnosis == "AML")

#merge with patMeta
#select top5 aberrations
i <- c(4:length(CoreDataSet))
CoreDataSet[,i] <- apply(CoreDataSet[,i], 2, function(x) as.numeric(as.character(x)))

top5 <- names(tail(sort(colSums(CoreDataSet[,i], na.rm=TRUE)), n=5))
top5

diag_merged_aml <- merge(diag, CoreDataSet, by.x = "patientID", by.y= "PatientID") |> 
  group_by(patientID) |> 
  slice_head(n=1) |> 
  as.data.frame()

#adjust the wide df  
wide.use <- names(wide)[(names(wide) %in% diag_merged_aml$patientID)]
wide.subset <- wide[, wide.use]
wide.subset.drug <- wide[, c("drug", wide.use)]

#merge with pathway table
df_pathway <- drugCategory_modified[, c("drug", "Pathway_mod")]

df_comb <- merge(df_pathway, wide.subset.drug, by = "drug", all = TRUE)

#scale
x <- data.frame(wide.subset)
y <- data.frame(t(scale(t(x))))

#create matrix
rownames(y) = wide.subset.drug$drug
colnames(y) = colnames(wide.subset)
z <- as.matrix(y)

#set colors
library(RColorBrewer)
library(paletteer)

#annotation for FLT3-ITD
`FLT3-ITD` <- diag_merged_aml[,c("FLT3-ITD")]
`FLT3-ITD` <- factor(`FLT3-ITD`, levels = c(1, 0,"unknown"))
flt3_anno <- `FLT3-ITD` %>% replace(is.na(.), "unknown")
flt3_colors <- c("black", "grey90", "#999999")
annotation_colors_flt3 <- setNames(flt3_colors, sort(unique(flt3_anno)))

#annotation for NPM1
NPM1 <- diag_merged_aml[,c("NPM1")]
NPM1 <- factor(NPM1, levels = c(1, 0,"unknown"))
NPM1_anno <- NPM1 %>% replace(is.na(.), "unknown")
NPM1_colors <- c("black", "grey90", "#999999")
annotation_colors_NPM1 <- setNames(NPM1_colors, sort(unique(NPM1_anno)))

#annotation for complex
`Complex karyotype (≥3)` <- diag_merged_aml[,c("complex caryotype\r\n( ≥3 aberrations)")]
`Complex karyotype (≥3)` <- factor(`Complex karyotype (≥3)`, levels = c(1, 0,"unknown"))
complex_anno <- `Complex karyotype (≥3)` %>% replace(is.na(.), "unknown")
complex_colors <- c("black", "grey90", "#999999")
annotation_colors_complex <- setNames(complex_colors, sort(unique(complex_anno)))

#annotation for +8
`+8` <- diag_merged_aml[,c("+8")]
`+8` <- factor(`+8`, levels = c(1, 0,"unknown"))
plus8_anno <- `+8` %>% replace(is.na(.), "unknown")
plus8_colors <- c("black", "grey90", "#999999")
annotation_colors_plus8 <- setNames(plus8_colors, sort(unique(plus8_anno)))

#annotation for inv(16)
`inv(16)` <- diag_merged_aml[,c("inv(16)")]
`inv(16)` <- factor(`inv(16)`, levels = c(1, 0,"unknown"))
inv16_anno <- `inv(16)` %>% replace(is.na(.), "unknown")
inv16_colors <- c("black", "grey90", "#999999")
annotation_colors_inv16 <- setNames(inv16_colors, sort(unique(inv16_anno)))


#annotation for pathways    
df_comb$Pathway_mod <- factor(df_comb$Pathway_mod,
                              levels = c("AKT/mTOR", "B-cell receptor", "Bromodomain/PLK", "Chemotherapy", 
                                         "DNA damage response", "JAK/STAT", "MAPK", "PI3K", "Stress pathway", "other"))

annotation_colors_row <- setNames(c("lightgrey", "#BC80BD", "#FCCDE5", "#B3DE69", "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", "#FFFFB3", "#8DD3C7"), c("other", "AKT/mTOR", "DNA damage response", "Chemotherapy" ,  "MAPK" , "B-cell receptor", "PI3K", "Stress pathway", "Bromodomain/PLK", "JAK/STAT"))

#z-score color
color_fun_reversed <- colorRamp2(c(4, median(z), -4), c("red", "white", "blue"))

#annotations
ha = HeatmapAnnotation(
  `FLT3-ITD` = flt3_anno,
  NPM1 = NPM1_anno,
  `Complex karyotype (≥3)` = complex_anno,
  `+8` = plus8_anno,
  `inv(16)` = inv16_anno,
  annotation_name_gp = gpar(fontsize = 8), 
  col = list(
    `FLT3-ITD` = annotation_colors_flt3,
    NPM1 = annotation_colors_NPM1,
    `Complex karyotype (≥3)` = annotation_colors_complex,
    `+8` = annotation_colors_plus8,
    `inv(16)` = annotation_colors_inv16
  ),
  annotation_legend_param = list(
    `FLT3-ITD` = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    NPM1 = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    `Complex karyotype (≥3)` = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    `+8`= list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    `inv(16)` = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8))
  ))

pway <- as.data.frame(df_comb[,1:2])
pway_ordered = pway[match(rownames(z), pway$drug), ]

ha2 = rowAnnotation(
  Pathway = pway_ordered$Pathway_mod, 
  annotation_name_gp = gpar(fontsize = 8), 
  col = list(Pathway = annotation_colors_row),
  annotation_legend_param = list(
    Pathway = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      nrow = 2
    )
  ),
  show_legend = FALSE
)

#heatmap
aml <- Heatmap(z, name = "Z-score", column_title = "AML", column_title_gp = gpar(fontsize = 8, fontface = "bold"), show_row_names = TRUE, show_column_names=FALSE, 
             row_names_gp = gpar(fontsize = 8), width = ncol(z)*unit(0.75, "mm"), 
             height = nrow(z)*unit(2.75, "mm"), col = color_fun_reversed, 
             top_annotation = ha, left_annotation = ha2, show_column_dend = FALSE, row_title = NULL, 
             #row_km = 5,
             #row_title_gp = gpar(fontsize = 8), 
             show_parent_dend_line = FALSE, show_heatmap_legend = FALSE)
