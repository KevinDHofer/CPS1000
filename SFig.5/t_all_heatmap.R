library(tidyverse)
library(dplyr)
library("ComplexHeatmap")
library(circlize)

#load raw data
usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

drugCategory_modified <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/drugCategory_modified.csv", sep=";")

CoreDataSet <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/heatmap_CoreDataSet_External_Summary_jlu.xlsx") |> 
  dplyr::select(-SampleID) |> 
  filter(Diagnosis == "T-ALL")

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
  filter(diagnosis == "T-ALL")
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
  filter(diagnosis == "T-ALL")

#merge with patMeta
#select top5 aberrations
i <- c(4:length(CoreDataSet))
CoreDataSet[,i] <- apply(CoreDataSet[,i], 2, function(x) as.numeric(as.character(x)))

top5 <- names(tail(sort(colSums(CoreDataSet[,i], na.rm=TRUE)), n=5))
top5

diag_merged_tall <- merge(diag, CoreDataSet, by.x = "patientID", by.y= "PatientID") |> 
  group_by(patientID) |> 
  slice_head(n=1) |> 
  as.data.frame()

#adjust the wide df  
wide.use <- names(wide)[(names(wide) %in% diag_merged_tall$patientID)]
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

#annotation for CDKN2A
CDKN2A <- diag_merged_tall[,c("CDKN2A")]
CDKN2A <- factor(CDKN2A, levels = c(1, 0,"unknown"))
CDKN2A_anno <- CDKN2A %>% replace(is.na(.), "unknown")
CDKN2A_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_CDKN2A <- CDKN2A_colors

#annotation for NOTCH1
NOTCH1 <- diag_merged_tall[,c("NOTCH1")]
NOTCH1 <- factor(NOTCH1, levels = c(1, 0,"unknown"))
NOTCH1_anno <- NOTCH1 %>% replace(is.na(.), "unknown")
NOTCH1_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_NOTCH1 <- NOTCH1_colors

#annotation for CDKN2B
CDKN2B <- diag_merged_tall[,c("CDKN2B")]
CDKN2B <- factor(CDKN2B, levels = c(1, 0,"unknown"))
CDKN2B_anno <- CDKN2B %>% replace(is.na(.), "unknown")
CDKN2B_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_CDKN2B <- CDKN2B_colors

#annotation for PTEN
PTEN <- diag_merged_tall[,c("PTEN")]
PTEN <- factor(PTEN, levels = c(1, 0,"unknown"))
PTEN_anno <- PTEN %>% replace(is.na(.), "unknown")
PTEN_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_PTEN <- PTEN_colors

#annotation for IL7R
IL7R <- diag_merged_tall[,c("IL7R")]
IL7R <- factor(IL7R, levels = c(1, 0,"unknown"))
IL7R_anno <- IL7R %>% replace(is.na(.), "unknown")
IL7R_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_IL7R <- IL7R_colors

#annotation for pathways    
df_comb$Pathway_mod <- factor(df_comb$Pathway_mod,
                              levels = c("AKT/mTOR", "B-cell receptor", "Bromodomain/PLK", "Chemotherapy", 
                                         "DNA damage response", "JAK/STAT", "MAPK", "PI3K", "Stress pathway", "other"))

annotation_colors_row <- setNames(c("lightgrey", "#BC80BD", "#FCCDE5", "#B3DE69", "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", "#FFFFB3", "#8DD3C7"), c("other", "AKT/mTOR", "DNA damage response", "Chemotherapy" ,  "MAPK" , "B-cell receptor", "PI3K", "Stress pathway", "Bromodomain/PLK", "JAK/STAT"))

#z-score color
color_fun_reversed <- colorRamp2(c(4, median(z), -4), c("red", "white", "blue"))

#annotations
ha = HeatmapAnnotation(
  CDKN2A = CDKN2A_anno,
  NOTCH1 = NOTCH1_anno,
  CDKN2B = CDKN2B_anno,
  PTEN = PTEN_anno,
  IL7R = IL7R_anno,
  annotation_name_gp = gpar(fontsize = 8), 
  col = list(
    CDKN2A = annotation_colors_CDKN2A,
    NOTCH1 = annotation_colors_NOTCH1,
    CDKN2B = annotation_colors_CDKN2B,
    PTEN = annotation_colors_PTEN,
    IL7R = annotation_colors_IL7R
  ),
  annotation_legend_param = list(
    CDKN2A = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    NOTCH1 = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    CDKN2B = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    PTEN = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    IL7R = list(
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
tall <- Heatmap(z, name = "Z-score", column_title = "T-ALL", column_title_gp = gpar(fontsize = 8, fontface = "bold"), show_row_names = TRUE, show_column_names=FALSE, 
             row_names_gp = gpar(fontsize = 8), width = ncol(z)*unit(0.75, "mm"), 
             height = nrow(z)*unit(2.75, "mm"), col = color_fun_reversed, 
             top_annotation = ha, left_annotation = ha2, show_column_dend = FALSE, row_title = NULL, 
             #row_km = 5,
             #row_title_gp = gpar(fontsize = 8), 
             show_parent_dend_line = FALSE, show_heatmap_legend = FALSE)

