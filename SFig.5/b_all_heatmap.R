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
  filter(Diagnosis == "B-ALL")

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
  filter(diagnosis == "B-ALL")
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
  filter(diagnosis == "B-ALL")

#merge with patMeta
#select top5 aberrations
i <- c(4:length(CoreDataSet))
CoreDataSet[,i] <- apply(CoreDataSet[,i], 2, function(x) as.numeric(as.character(x)))

top5 <- names(tail(sort(colSums(CoreDataSet[,i], na.rm=TRUE)), n=5))
top5
#plus BCR-ABL

diag_merged_ball <- merge(diag, CoreDataSet, by.x = "patientID", by.y= "PatientID") |> 
  group_by(patientID) |> 
  slice_head(n=1) |> 
  as.data.frame()

#adjust the wide df  
wide.use <- names(wide)[(names(wide) %in% diag_merged_ball$patientID)]
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

#annotation for KRAS
KRAS <- diag_merged_ball[,c("KRAS")]
KRAS <- factor(KRAS, levels = c(1, 0,"unknown"))
KRAS_anno <- KRAS %>% replace(is.na(.), "unknown")
KRAS_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_KRAS <- KRAS_colors

#annotation for SH2B3
SH2B3 <- diag_merged_ball[,c("SH2B3")]
SH2B3 <- factor(SH2B3, levels = c(1, 0,"unknown"))
SH2B3_anno <- SH2B3 %>% replace(is.na(.), "unknown")
SH2B3_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_SH2B3 <- SH2B3_colors

#annotation for TP53
TP53 <- diag_merged_ball[,c("TP53")]
TP53 <- factor(TP53, levels = c(1, 0,"unknown"))
TP53_anno <- TP53 %>% replace(is.na(.), "unknown")
TP53_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_TP53 <- TP53_colors

#annotation for IKFZ1
IKFZ1 <- diag_merged_ball[,c("IKFZ1")]
IKFZ1 <- factor(IKFZ1, levels = c(1, 0,"unknown"))
IKFZ1_anno <- IKFZ1 %>% replace(is.na(.), "unknown")
IKFZ1_colors <- c("1" = "black", "0" = "grey90", "unknown" = "#999999")
annotation_colors_IKFZ1 <- IKFZ1_colors
                                
#annotation for BCR-ABL
`BCR-ABL` <- diag_merged_ball[,c("BCR-ABL")]
`BCR-ABL` <- factor(`BCR-ABL`, levels = c(1, 0,"unknown"))
bcr_anno <- `BCR-ABL` %>% replace(is.na(.), "unknown")
bcr_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(bcr_anno)))
annotation_colors_bcr <- setNames(bcr_colors, sort(unique(bcr_anno)))

#annotation for pathways    
df_comb$Pathway_mod <- factor(df_comb$Pathway_mod,
                              levels = c("AKT/mTOR", "B-cell receptor", "Bromodomain/PLK", "Chemotherapy", 
                                         "DNA damage response", "JAK/STAT", "MAPK", "PI3K", "Stress pathway", "other"))

annotation_colors_row <- setNames(c("lightgrey", "#BC80BD", "#FCCDE5", "#B3DE69", "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", "#FFFFB3", "#8DD3C7"), c("other", "AKT/mTOR", "DNA damage response", "Chemotherapy" ,  "MAPK" , "B-cell receptor", "PI3K", "Stress pathway", "Bromodomain/PLK", "JAK/STAT"))

#z-score color
color_fun_reversed <- colorRamp2(c(4, median(z), -4), c("red", "white", "blue"))

#annotations
ha = HeatmapAnnotation(
  `BCR-ABL` = bcr_anno,
  KRAS = KRAS_anno,
  SH2B3 = SH2B3_anno,
  TP53 = TP53_anno,
  IKFZ1 = IKFZ1_anno,
  annotation_name_gp = gpar(fontsize = 8), 
  col = list(
    KRAS = annotation_colors_KRAS,
    SH2B3 = annotation_colors_SH2B3,
    TP53 = annotation_colors_TP53,
    IKFZ1 = annotation_colors_IKFZ1,
    `BCR-ABL` = annotation_colors_bcr
  ),
  annotation_legend_param = list(
    `BCR-ABL` = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    KRAS = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    SH2B3 = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    TP53 = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)),
    IKFZ1 = list(
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
      ncol = 1
    )
  ),
  show_legend = FALSE
)

#heatmap
ball <- Heatmap(z, name = "Z-score", column_title = "B-ALL", column_title_gp = gpar(fontsize = 8, fontface = "bold"), show_row_names = TRUE, show_column_names=FALSE, 
             row_names_gp = gpar(fontsize = 8), width = ncol(z)*unit(0.75, "mm"), 
             height = nrow(z)*unit(2.75, "mm"), col = color_fun_reversed, 
             top_annotation = ha, left_annotation = ha2, show_column_dend = FALSE, row_title = NULL, 
             #row_km = 5,
             #row_title_gp = gpar(fontsize = 8), 
             show_parent_dend_line = FALSE, show_heatmap_legend = FALSE)
