library(tidyverse)
library(readxl)
library(dplyr)
library(patchwork)
library(ggvenn)
library(ComplexHeatmap)
library(RColorBrewer)
library(factoextra)
library(ggpubr)
library(circlize)

#load data
samples <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSampleList_longitudinal.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viabilityTable_allSamples.csv")

treatments <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedPairs_analysis.csv")

#subtract suffixes from sampleID
treatments$sample1 <- gsub("-[1234]$", "", treatments$sample1)
treatments$sample2 <- gsub("-[1234]$", "", treatments$sample2)

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/patmeta_210324.RData")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/survival_190516.RData")

# merge tables
merged <- left_join(samples, viability_table[, !colnames(viability_table) %in% "patientID"], 
                    by = "sampleID") 

#viability table: 285 samples from 113 patients
viab_list_samples <- unique(merged$sampleID)
viab_list_pat <- unique(merged$patientID)

#add date of sampling
#subtract suffixes from sampleID
merged$sampleID <- gsub("-[1234]$", "", merged$sampleID)
#verify sampleID is unique
length(unique(merged$sampleID)) == nrow(merged)

merged_date <- left_join(merged, survT[, colnames(survT) %in% c("sampleID", "sampleDate")],
                         by = "sampleID") |> 
  relocate(sampleDate, .after = sampleID)

#check NA
sum(is.na(merged_date$sampleDate))

#nr of analyzed pairs: 285
nrow(merged_date)

#transfrom in long data frame
merged_table_long <- merged_date |>
  pivot_longer(cols = c(4:length(merged_date)), names_to = "drug", values_to = "viability")
merged_table_long

#add concentration level and remove last two characters from each string in 'drug' column
merged_table_long$concentration <- rep(seq(1:5), times=nrow(merged_table_long)/5)
merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
merged_table_long

#create means
means_viab <- merged_table_long |> 
  group_by(patientID, sampleID, sampleDate, drug) |> 
  summarize(viab=mean(viability))
means_viab

# Convert to wide format
means_viab_wide <- means_viab %>%
  pivot_wider(names_from = drug, values_from = viab)

#calculate pair-wise viability differences
create_long_differences <- function(treatments_df, viab_df) {
  # Get drug columns
  drug_cols <- setdiff(names(viab_df), c("patientID", "sampleID", "sampleDate"))
  
  # Get sample1 values
  sample1_long <- treatments_df %>%
    left_join(viab_df, by = c("sample1" = "sampleID")) %>%
    dplyr::select(patientID.x, pair, sample1, sample2, status, therapy, all_of(drug_cols)) %>%
    rename(patientID = patientID.x) %>%
    pivot_longer(cols = all_of(drug_cols), names_to = "drug", values_to = "sample1_value")
  
  # Get sample2 values
  sample2_long <- treatments_df %>%
    left_join(viab_df, by = c("sample2" = "sampleID")) %>%
    dplyr::select(patientID.x, pair, all_of(drug_cols)) %>%
    rename(patientID = patientID.x) %>%
    pivot_longer(cols = all_of(drug_cols), names_to = "drug", values_to = "sample2_value")
  
  # Combine and calculate differences
  differences_long <- sample1_long %>%
    left_join(sample2_long, by = c("patientID", "pair", "drug")) %>%
    mutate(
      difference = sample2_value - sample1_value,
      fold_change = sample2_value / sample1_value,
      percent_change = ((sample2_value - sample1_value) / sample1_value) * 100
    )
  
  return(differences_long)
}

differences_long <- create_long_differences(treatments, means_viab_wide)

viab_list_pat <- differences_long |> 
  filter(drug == "Ibrutinib") |> 
  dplyr::select(1:6)

#write.csv(viab_list_pat, "viab_list_pat.csv")

#summary
#nr of patients: 110
length(unique(viab_list_pat$patientID))

#nr of sample pairs: 144
nrow(viab_list_pat)

#nr of patients without therapy: 65
viab_list_pat |> filter(therapy == "none") |> 
  nrow()

#nr of patients with anytherapy: 79
viab_list_pat |> filter(therapy != "none") |> 
  nrow()

#nr of patients with chemo: 40
viab_list_pat |> filter(therapy == "chemo")|> 
  distinct(patientID) |> 
  nrow()

#nr of patients with targeted therapy: 19
viab_list_pat |> filter(therapy == "targeted") |> 
  distinct(patientID) |> 
  nrow()
  

#plot
p1 <- differences_long |> 
  ggplot(aes(difference))+
    geom_histogram(color="white", fill="grey40", alpha=0.5,bins = 50)+
    geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "black")+
  labs(x="Difference in mean viability", y="Count")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(size=8), 
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'))+
  scale_y_continuous(expand = c(0.02,0))+
  scale_x_continuous(expand = c(0.02,0), limits=c(-0.3, 0.3))
p1

Fig6a <- p1 + theme(plot.margin = margin(t=1, b=1, l=1, unit = "cm"))

#stats
#date intervals of samples (only available for 131 pairs)
df_date <- differences_long |> 
  filter(drug == "Ibrutinib") |> 
  dplyr::select("patientID", "pair", "sample1", "sample2") |> 
  merge(merged_date[, colnames(merged_date) %in% c("sampleID", "sampleDate")], by.x = "sample1", by.y ="sampleID") |> 
  rename(sampleDate1 = sampleDate)
df_date <- df_date |> 
  rename(sampleID = sample2) |> 
  left_join(merged_date[, colnames(merged_date) %in% c("sampleID", "sampleDate")], by ="sampleID") |> 
  rename(sampleDate2 = sampleDate, sample2 = sampleID) |> 
  mutate(diff_date = sampleDate2 - sampleDate1)

df_date$diff_date <- as.numeric(df_date$diff_date)
summary(df_date$diff_date)

#percent of sample pairs with more than 0.1 change in viability
differences_long |> 
  filter(abs(difference) > 0.1) |> 
  nrow()/nrow(differences_long)*100

#t-test with FDR within status and therapy
df_ttest <- differences_long |> 
  group_by(drug, status, therapy) |> 
  mutate(pval = t.test(sample1_value, sample2_value, paired = TRUE)$p.value) |> 
  ungroup() |>
  group_by(status, therapy) |>  # Group by status and therapy for adjustment
  mutate(padj = p.adjust(pval, method="BH"))

#compare treatment groups
df_ttest <- df_ttest |> 
  mutate(sign = ifelse(padj<0.1, "change", "no change"))

#set colors
colorList <- c(colorRampPalette(brewer.pal(3, 'Set2'))(3), 
               "grey80")
names(colorList) = c(unique(df_ttest$therapy), "Below 10% FDR")

color_sign <- c(colorRampPalette(brewer.pal(2, 'Set3'))(2))
names(color_sign) = c(unique(df_ttest$sign))

color_sign["no change"] <- "grey80"

#calculate chi square test
chi_test <- chisq.test(table(df_ttest$therapy, df_ttest$sign))
p_value <- chi_test$p.value
chi_stat <- chi_test$statistic

p22 <- df_ttest |> 
  ggplot(aes(x=therapy, fill=sign))+
  geom_bar(position="fill") +
  geom_text(
    aes(label=round(..count.. /tapply(..count.., ..x.., sum)[as.character(..x..)], 2)),
    stat="count",
    position=position_fill(vjust=0.5),
    size=3) +
  theme_classic()+
  scale_fill_manual(values = color_sign, labels = c("FDR<0.1", "No sign. change"))+
  theme(axis.text.y = element_text(color="black", size=8), 
        axis.text.x = element_text(color="black", size=8),
        axis.title = element_text(color="black", size=8), 
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8),
        plot.title = element_text(hjust=0.5, size=8, face="bold"))+
  labs(title="All treatments", fill="Drug response\nevolution",
       x="Treatment", y="Proportion of sample pair-drug interactions")+
  scale_x_discrete(labels=c("None", "Chemotherapy", "Targeted"))+
  annotate("text", x = Inf, y = Inf, 
           label = paste0("χ² = ", round(chi_stat, 3), 
                          ", p = ", format.pval(p_value, digits = 3)),
           hjust = 2, vjust = 1, size = 3)
p22

# Percentage of significant changes by therapy condition
sig_summary <- df_ttest %>%
  group_by(therapy) %>%
  summarise(
    total_tests = n(),
    significant_tests = sum(padj < 0.1, na.rm = TRUE),
    percent_significant = (significant_tests / total_tests) * 100,
    .groups = 'drop'
  )

p23 <- sig_summary |> 
  ggplot(aes(x=therapy, y=percent_significant, fill=therapy))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(color="black", size=8), 
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=8, face="bold"))+
  labs(title="All treatments", fill="Drug response\nevolution",
       x="", y="% significantly altered sampling pairs (FDR <0.1)")+
  scale_x_discrete(labels=c("Chemotherapy", "None", "Targeted"))+
  annotate("text", 
           x = length(unique(sig_summary$therapy)) / 2 + 0.5, # Center horizontally
           y = max(sig_summary$percent_significant) * 1.1,  # Position near top
           label = paste0("χ² = ", round(chi_stat, 3), ", p = ", format.pval(p_value, digits = 3)), 
           hjust = 0.5,  # Center horizontally
           vjust = 1,    # Align to top
           size = 3)+
  scale_fill_manual(values=colorList)
p23
SFig10a <-p23

#differences in targeted group, differentiated by status
df_ttest_t <-df_ttest |> 
  filter(therapy == "targeted")

# Percentage by status in targeted
sig_summary_targeted <- df_ttest_t %>%
  group_by(status) %>%
  summarise(
    total_tests = n(),
    significant_tests = sum(padj < 0.1, na.rm = TRUE),
    percent_significant = (significant_tests / total_tests) * 100,
    .groups = 'drop'
  ) |> 
  arrange(desc(percent_significant))

sig_summary_targeted$status <- factor(sig_summary_targeted$status, levels = sig_summary_targeted$status)

chi_test_t <- chisq.test(table(df_ttest_t$therapy, df_ttest_t$sign))
p_value_t <- chi_test_t$p.value
chi_stat_t <- chi_test_t$statistic

#set colors
colorListTargeted <- setNames(colorRampPalette(brewer.pal(3, 'Paired'))(3), unique(df_ttest_t$status))

p24 <- df_ttest |> 
  filter(therapy =="targeted") |> 
  ggplot(aes(x=status, fill=sign))+
  geom_bar(position="fill") +
  geom_text(
    aes(label=round(..count.. /tapply(..count.., ..x.., sum)[as.character(..x..)], 2)),
    stat="count",
    position=position_fill(vjust=0.5),
    size=3) +
  theme_classic()+
  scale_fill_manual(values = color_sign, 
                    labels = c("FDR<0.1", "No sign. change")
                    )+
  theme(axis.text.y = element_text(color="black", size=8), 
        axis.text.x = element_text(color="black", size=8),
        axis.title = element_text(color="black", size=8), 
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8),
        plot.title = element_text(hjust=0.5, size=8, face="bold"))+
  labs(title="Targeted therapy", fill="Drug response\nevolution",
       x="Timing", y="Proportion of sample pair-drug interactions")+
  scale_x_discrete(labels=c("Treatment start\nbetween samples", "Treatment start and stop\nbetween samples", "Sampling during\ntreatment"))+
  annotate("text", x = Inf, y = Inf, 
           label = paste0("χ² = ", round(chi_stat_t, 3), ", p = ", 
                          format.pval(p_value_t, digits = 3)),
           hjust = 2, vjust = 1, size = 3)
p24

p25 <- sig_summary_targeted |> 
  ggplot(aes(x=status, y=percent_significant, fill=status))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(color="black", size=8), 
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=8, face="bold"))+
  labs(title="Timing of targeted therapy", fill="Drug response\nevolution",
       x="", y="% significantly altered sample pairs (FDR <0.1)")+
  annotate("text", 
           x = length(unique(sig_summary_targeted$status)) / 2 + 0.5, # Center horizontally
           y = max(sig_summary_targeted$percent_significant) * 1.1,  # Position near top
           label = paste0("χ² = ", round(chi_stat_t, 3), ", p = ", format.pval(p_value_t, digits = 3)), 
           hjust = 0.5,  # Center horizontally
           vjust = 1,    # Align to top
           size = 3)+
  #scale_x_discrete(labels=c("Start between\nsamples", "Sampling during\ntreatment", "Start and stop\nbetween samples"))+
  scale_fill_manual(values=colorListTargeted)
p25
SFig10b <-p25

#nr of significant changes
df_ttest_therapy <- df_ttest |> 
  filter(padj <0.1) |> 
  group_by(patientID, pair, sample1, sample2, therapy) |> 
  summarize(nr = sum(padj<0.1), .groups = "drop")

#plot: number of sample pairs with significant effects
p12 <- df_ttest_therapy |> 
  group_by(therapy) |> 
  summarize(nr_max = max(nr)) |> 
  ggplot(aes(x=therapy, y=nr_max, fill=therapy))+
  geom_col()+
  theme_classic()+
#scale_fill_manual(values = color_treatments)+
  theme(axis.text.y = element_text(color="black", size=8), 
        axis.text.x = element_text(color="black", size=8),
        axis.title = element_text(color="black", size=8), 
        legend.position = "none")+
  labs(title="", 
       x="", y="Nr. of drugs with sign. change (FDR<0.1)")+
  scale_x_discrete(labels=c("None", "Chemotherapy", "Targeted"))
p12
  
#p-value scatterplot for drugs with FDR <0.1
shapeList <- setNames(c(0,1,2,3,4,5), unique(df_ttest$status))

#order drugs by max p
pval_list <- df_ttest |> 
  filter(padj <0.1) |> 
  group_by(drug) |> 
  arrange(pval, descending = FALSE)
drugOrder <- unique(pval_list$drug)

p3 <- df_ttest |> 
  filter(drug %in% drugOrder) |> 
  mutate(drug = factor(drug, levels = drugOrder),
  #status = ifelse(padj < 0.1, status, "Below 10% FDR"),
         therapy = ifelse(padj < 0.1, therapy, "Below 10% FDR")) |> 
  ggplot(aes(y=-log10(pval), x=drug, color=therapy, shape=status))+
  geom_point(size =2, alpha = 0.8)+
  theme_bw()+
  scale_color_manual(values = colorList, labels=c("Below 10% FDR",
                                                  "Chemotherapy",
                                                  "None",
                                                  "Targeted"))+
  scale_shape_manual(values = shapeList)+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust=0.5, 
                                   color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), 
        axis.title.y = element_text(size=8), 
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8),
        plot.title = element_text(size=8, face="bold", vjust=0.5),
        panel.grid.major.x = element_line(size=.1, color="grey50", 
                                          linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(title="", color="Treatment", shape = "Timing",
       x="", y="-log10(p)")
p3

#make volcano
library(ggrepel)

test <- df_ttest |> 
  mutate(therapy = ifelse(padj < 0.1, therapy, "Below 10% FDR")) |> 
  group_by(drug, therapy, status) |> 
  summarise(med_diff = median(difference), pval = median(pval), .groups = "drop")

p4 <- test |> 
  mutate(status=ifelse(status %in% c("treatment between", "start between"), status, "other")) |> 
  ggplot(aes(y=-log10(pval), x=med_diff, color=therapy, shape=status
             ))+
  geom_point(size =2, alpha = 0.8)+
  theme_classic()+
  geom_vline(xintercept = 0, linetype = "dashed", color="grey50")+
  geom_text_repel(aes(label= ifelse(therapy != "Below 10% FDR", drug, "")),
                  max.overlaps = 10,
                  color="black",
                  min.segment.length = 1,
                  point.padding = 0.5,
                  box.padding = 0.5,
                  size=3)+
  scale_color_manual(values = colorList, labels=c("Below 10% FDR",
                                                  "Chemotherapy",
                                                  "None",
                                                  "Targeted"))+
  #scale_shape_manual(values = shapeList)+
  theme(axis.text.x = element_text(color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), 
        axis.title = element_text(size=8), 
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8),
        plot.title = element_text(size=8, face="bold", vjust=0.5))+
  labs(title="", color="Treatment", shape = "Timing",
       x="Median difference in viability", y="-log10(p)")+
  xlim(-0.1,0.1)+
  scale_shape_discrete(labels = c("start between" = "Start between samples", 
                                "treatment between" = "Start and stop between samples",
                                "other" =  "Other"))
p4

#median difference in untreated patients with sign. changes
df_ttest |> 
  filter(therapy == "none", padj <0.1) |> 
  #group_by(drug) |> 
  summarize(med_diff = median(abs(difference)))

#heatmap for all patients
df_matrix <- differences_long |> 
  dplyr::select(patientID, drug, status, therapy, difference) |> 
  pivot_wider(names_from = c(patientID, status, therapy), values_from = difference)

#scale
x <- data.frame(df_matrix)
rownames(x) <- x$drug
x <- x |> 
  dplyr::select(-drug)
y <- data.frame(t(scale(t(x))))

#create matrix
z <- as.matrix(y)

#prepare annotation table
therapy_list <- as.data.frame(t(z)) |> 
  mutate(treatments = rownames(t(z)), .before=1) |> 
  mutate(therapy = ifelse(
    str_detect(treatments, "chemo"), "chemo", 
      ifelse(str_detect(treatments, "targeted"), "targeted", "none")), .before = 2) |> 
  mutate(status = ifelse(
    str_detect(treatments, "end.between"), "end between", 
      ifelse(str_detect(treatments, "once.treated"), "once treated",
          ifelse(str_detect(treatments, "start.between"), "start between", 
             ifelse(str_detect(treatments, "treatment.between"), "treatment between", 
                 ifelse(str_detect(treatments, "untreated"), "untreated", "within treatment"
                ))))), .before = 3) |> 
  dplyr::select(therapy, status)

#set colors
therapy_colors <- setNames(colorRampPalette(brewer.pal(3, "Set2"))(length(unique(therapy_list$therapy))),
                           unique(therapy_list$therapy))

status_colors <- setNames(colorRampPalette(brewer.pal(6, "Accent"))(length(unique(therapy_list$status))),
                          unique(therapy_list$status))
status_colors["untreated"] <- "grey80" 
status_colors["within treatment"] <-"#FFFF99"

color_fun_reversed <- colorRamp2(c(4, median(z), -4), c("red", "white", "blue"))

#annotations
ha = HeatmapAnnotation(
  Treatment = therapy_list$therapy, 
  Timing = therapy_list$status,
  annotation_name_gp = gpar(fontsize = 8), 
  col = list(Treatment = therapy_colors,
             Timing = status_colors),
  annotation_legend_param = list(
    Treatment = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      labels = c("Chemotherapy", "None", "Targeted"), # Custom labels
      ncol = 1
    ),
    Timing = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      labels = c("End between", "Once treated", "Start between", "Treatment between", "Untreated", "Within treatment"), # Custom labels
      ncol = 1
    ))
  )

#draw heatmap
set.seed(123)
h <- Heatmap(z, name = "Z-score", show_row_names = TRUE, show_column_names=FALSE, 
            row_names_gp = gpar(fontsize = 8), width = ncol(z)*unit(0.8, "mm"), height = nrow(z)*unit(2.5, "mm"), col = color_fun_reversed, 
            top_annotation = ha, 
            show_column_dend = TRUE, row_title = NULL, column_title = NULL, 
            column_km = 3,
            row_km =2,
            show_parent_dend_line = FALSE, heatmap_legend_param = list(
              direction = "vertical", title_gp = gpar(fontsize = 8, fontface = "bold"), 
              labels_gp = gpar(fontsize = 8)))

draw(h, merge_legends = TRUE,
     legend_gap = unit(0.5, "cm"))

#xyz <- grid.grabExpr(draw(h, merge_legends = TRUE, heatmap_legend_side="bottom", padding = unit(c(10, 10, 0, 0), "mm")))
#xyz

#heatmap for targeted and untreated only
df_matrix2 <- differences_long |> 
  filter(therapy != "chemo", status == "start between" | status == "untreated") |> 
  dplyr::select(patientID, drug, status, therapy, difference) |> 
  pivot_wider(names_from = c(patientID, status, therapy), values_from = difference)

#scale
x2 <- data.frame(df_matrix2)
rownames(x2) <- x2$drug
x2 <- x2 |> 
  dplyr::select(-drug)
y2 <- data.frame(t(scale(t(x2))))

#create matrix
z2 <- as.matrix(y2)

#prepare annotation table
therapy_list2 <- as.data.frame(t(z2)) |> 
  mutate(treatments = rownames(t(z2)), .before=1) |> 
  mutate(therapy = ifelse(
    str_detect(treatments, "chemo"), "chemo", 
    ifelse(str_detect(treatments, "targeted"), "targeted", "none")), .before = 2) |> 
  mutate(status = ifelse(
    str_detect(treatments, "end.between"), "end between", 
    ifelse(str_detect(treatments, "once.treated"), "once treated",
           ifelse(str_detect(treatments, "start.between"), "start between", 
                  ifelse(str_detect(treatments, "treatment.between"), "treatment between", 
                         ifelse(str_detect(treatments, "untreated"), "untreated", "within treatment"
                         ))))), .before = 3) |> 
  dplyr::select(therapy, status)

#set colors
therapy_colors2 <- setNames(c("#FC8D62", "grey"), unique(therapy_list2$therapy))

status_colors2 <- setNames(colorRampPalette(brewer.pal(3, "Set3"))(length(unique(therapy_list2$status))),
                          unique(therapy_list2$status))

status_colors2["start between"] <- "#A6CEE3"

color_fun_reversed <- colorRamp2(c(4, median(z2), -4), c("red", "white", "blue"))

#annotations
ha2 = HeatmapAnnotation(
  Treatment = therapy_list2$therapy, 
  Timing = therapy_list2$status,
  annotation_name_gp = gpar(fontsize = 8), 
  col = list(Treatment = therapy_colors2,
             Timing = status_colors2),
  annotation_legend_param = list(
    Treatment = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      labels = c("None", "Targeted"), # Custom labels
      ncol = 1
    ),
    Timing = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      labels = c("Start between", "Untreated"), # Custom labels
      ncol = 1
    ))
)

#draw heatmap
set.seed(123)
h2 <- Heatmap(z2, name = "Z-score", show_row_names = TRUE, show_column_names=FALSE, 
             row_names_gp = gpar(fontsize = 8), width = ncol(z)*unit(0.5, "mm"), height = nrow(z)*unit(2.6, "mm"), col = color_fun_reversed, 
             top_annotation = ha2, 
             show_column_dend = FALSE, cluster_columns = FALSE, row_title = NULL, column_title = NULL, 
             #column_km = 2,
             row_km =2,
             column_split = therapy_list2$therapy,
             show_parent_dend_line = FALSE, heatmap_legend_param = list(
               direction = "vertical", title_gp = gpar(fontsize = 8, fontface = "bold"), 
               labels_gp = gpar(fontsize = 8)))

draw(h2, merge_legends = TRUE,
     legend_gap = unit(0.5, "cm"))

#boxplot for direction and magnitude of effects
median_ranking_targeted <- df_ttest |> 
  filter(padj <0.1, therapy == "targeted") |> 
  group_by(drug, therapy) |> 
  mutate(median_diff = median(difference)) |> 
  arrange(median_diff)
median_ranking_order <- unique(median_ranking_targeted$drug)

p6 <- df_ttest |> 
  filter(padj <0.1, therapy == "targeted") |> 
  mutate(drug = factor(drug, levels = median_ranking_order)) |> 
ggplot(aes(y=difference, x=drug, color = therapy))+ 
         geom_boxplot(alpha = 0.8)+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  scale_color_manual(values = colorList, labels=c("Chemotherapy",
                                                  "None",
                                                  "Targeted"))+
  theme_bw()+
  theme(axis.text.y = element_text(color="black", size=8), 
        axis.text.x = element_text(angle=90, hjust = 1, vjust =0.5, color="black", size=8),
        axis.title = element_text(size=8), 
        plot.title = element_text(size=8, face="bold", hjust=0.5),
        #legend.title = element_text(size=8, face="bold"), 
        #legend.text=element_text(size=8),
        #legend.key.size = unit(0.5, 'cm'),
        legend.position = "none")+
  labs(title="Targeted, start between", color="Treatment",
       y="Difference in mean viability", x="")
p6 + coord_flip()+
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0))

### individual drug responses
#create means
viab <- merged_table_long |> 
  group_by(patientID, sampleID, sampleDate, drug) 
viab

# Convert to wide format
viab_wide <- viab %>%
  pivot_wider(names_from = drug, values_from = viability)

#calculate pair-wise viability differences
create_long_differences <- function(treatments_df, viab_df) {
  # Get drug columns
  drug_cols <- setdiff(names(viab_wide), c("patientID", "sampleID", "sampleDate", "concentration"))
  
  # Get sample1 values
  sample1_long <- treatments_df %>%
    left_join(viab_df, by = c("sample1" = "sampleID")) %>%
    dplyr::select(patientID.x, pair, sample1, sample2, status, therapy, concentration, all_of(drug_cols)) %>%
    rename(patientID = patientID.x) %>%
    pivot_longer(cols = all_of(drug_cols), names_to = "drug", values_to = "sample1_value")
  
  # Get sample2 values
  sample2_long <- treatments_df %>%
    left_join(viab_df, by = c("sample2" = "sampleID")) %>%
    dplyr::select(patientID.x, pair, concentration, all_of(drug_cols)) %>%
    rename(patientID = patientID.x) %>%
    pivot_longer(cols = all_of(drug_cols), names_to = "drug", values_to = "sample2_value")
  
  # Combine and calculate differences
  differences_long <- sample1_long %>%
    left_join(sample2_long, by = c("patientID", "pair", "drug", "concentration")) %>%
    mutate(
      difference = sample2_value - sample1_value,
      fold_change = sample2_value / sample1_value,
      percent_change = ((sample2_value - sample1_value) / sample1_value) * 100
    )
  
  return(differences_long)
}

differences_long_allconc <- create_long_differences(treatments, viab_wide)

#nr of drug/sample pair
differences_long_allconc |> 
  group_by(patientID, pair, drug) |> 
  nrow()

#t-test with FDR within status and therapy
df_ttest_allconc <- differences_long_allconc |> 
  group_by(drug, concentration, status, therapy) |> 
  mutate(pval = t.test(sample1_value, sample2_value, paired = TRUE)$p.value) |> 
  ungroup() |>
  group_by(status, therapy) |>  # Group by status and therapy for  adjustment
  mutate(padj = p.adjust(pval, method="BH"))

df_ttest_allconc_sig <- df_ttest_allconc

#Order drugs by their association similarity (hierarchical clustering)
pMat <- df_ttest_allconc_sig |> 
  filter(padj < 0.1) |>
  group_by(drug, concentration, therapy, status) |> 
  mutate(signed_log_pval = -log10(mean(pval)) * sign(mean(percent_change))) |>
  ungroup() |> 
  dplyr::select(drug, concentration, therapy, status, signed_log_pval) |>
  # Create a unique identifier for each condition
  unite("condition", therapy, status, concentration, sep = "_") |> 
  pivot_wider(names_from = condition, 
              values_from = signed_log_pval,
              values_fn = median) |>
  as.data.frame() |>
  column_to_rownames("drug") |> 
  as.matrix()

# Replace NAs with 0 (no significant effect)
pMat[is.na(pMat)] <- 0

# Perform hierarchical clustering
set.seed(123)
drug_clusters <- hclust(dist(pMat), method = "ward.D2")
drug_order <- rownames(pMat)[drug_clusters$order]

df_ttest_allconc_sig$drug <- factor(df_ttest_allconc_sig$drug, levels = drug_order)

#rename variables
df_ttest_allconc_sig$therapy <- recode(df_ttest_allconc_sig$therapy,
                                   "chemo" = "Chemotherapy",
                                   "targeted" = "Targeted", 
                                   "none" = "None")

#plot heatmap
Fig6c <- df_ttest_allconc_sig |> 
  filter(padj <0.1) |>
  group_by(drug, concentration, therapy, status) |> 
  mutate(signed_log_pval = -log10(mean(pval)) * sign(mean(percent_change))) |>  
  ggplot(aes(x=concentration, y=drug, fill=signed_log_pval))+
  geom_tile(size = 0.2, color = "white")+
  facet_wrap(therapy~status, ncol=4)+
  scale_fill_gradient2(low="blue", high="red", limits=c(-4, 4))+
  theme_bw()+
  labs(y="", x="Concentration index")+
  theme(strip.text.y = element_text(size = 0), 
                       strip.text.x = element_text(size=8),
                       axis.text = element_text(size =8, color="black"), axis.title = element_text(size =8, color="black"),
                       panel.spacing = unit(0.15, "lines"), 
                       legend.title = element_text(size=8, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(color = "black", linewidth = 1, fill = "grey90"),        
        plot.margin = margin(l = 0.5, r = 1, unit = "cm"))+
  labs(fill = "**-log<sub>10</sub>*P*<br>with direction**") +
  theme(legend.title = ggtext::element_markdown())
Fig6c <- Fig6c + theme(plot.margin = margin(l=1, t=1, b=1, unit = "cm"))
Fig6c

#create list of sign. drugs in targeted, start between
list <- df_ttest_allconc_sig |> 
  filter(padj <0.1, therapy == "Targeted", status == "start between") |> 
  group_by(drug) |> 
  mutate(median_diff = median(difference)) |> 
  arrange(median_diff)
list$drug <- factor(list$drug, levels = unique(list$drug))

#boxplot for direction and magnitude of effects
median_ranking_targeted_allconc <- df_ttest |> 
  filter(drug %in% unique(list$drug), therapy == "targeted", status == "start between") |> 
  group_by(drug) |> 
  mutate(median_diff = median(difference)) |> 
  arrange(median_diff)
median_ranking_order_allconc <- unique(median_ranking_targeted_allconc$drug)


# Calculate medians and add color variable
# load drug class
drugCat <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/drugCatagory.csv", sep=";", comment.char="#")
drugCat$Drug <- gsub("-", ".", drugCat$Drug)

#set colors
annotation_colors_row <- setNames(c("lightgrey", "#FCCDE5", "#B3DE69", "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", 
                                    "#FFFFB3", "#8DD3C7", "#A6CEE3", "#6DBD57", "#DE9E83", "#BC80BD"), 
                                  c("other", "DNA damage response", "Chemotherapy" ,  "MAPK" , "B-cell receptor", 
                                    "PI3K", "Stress pathway", "Bromodomain/PLK", "JAK/STAT", "Apoptosis", 
                                    "Cell cycle control", "Histone methylation", "PI3K/AKT/mTOR"))

df_with_median_color <- df_ttest |> 
  filter(drug %in% unique(median_ranking_targeted_allconc$drug), 
         therapy == "targeted", 
         status == "start between") |> 
  mutate(drug = factor(drug, levels = median_ranking_order_allconc)) |>
  group_by(drug) |>
  mutate(median_val = median(difference, na.rm = TRUE),
         median_color = ifelse(median_val > 0, "Above 0", "Below 0")) |>
  ungroup() |> 
  merge(drugCat[,c("Drug", "Pathway")], by.x="drug", by.y="Drug") |> 
  mutate(pathway_col = case_when(Pathway %in% c("EGFR", "HGF", "PKC", "NFkB", "Cytoskeleton", "Cytokine receptor") ~ "other",
                                 Pathway == "ABL (BCR)" ~ "B-cell receptor",
                                 Pathway == "DDR" ~ "DNA damage response",
                                 Pathway == "PI3K" ~ "PI3K/AKT/mTOR",
                                 Pathway %in% c("bromodomain", "PLK") ~ "Bromodomain/PLK",
                                 TRUE ~ Pathway),
         pathway_col = str_trim(pathway_col))

df_with_median_color$pathway_col <- factor(df_with_median_color$pathway_col, 
                                           levels = c("Apoptosis", "B-cell receptor", "Bromodomain/PLK", 
                                                      "Cell cycle control", "DNA damage response", "JAK/STAT", 
                                                      "MAPK", "PI3K/AKT/mTOR", "Stress pathway", "other"))

p6 <- df_with_median_color |>
  ggplot(aes(y = difference, x = drug, fill=pathway_col)) + 
  geom_boxplot(alpha = 0.75, outliers = TRUE) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth =0.5) +
  scale_fill_manual(values=annotation_colors_row)+
  theme_classic() +
  theme(axis.text.y = element_text(color = "black", size = 8), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 8),
        axis.title = element_text(size = 8), 
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'),) +
  labs(
    #title = "Targeted treatment (start between)", 
       fill = "Pathway",
       y = "Difference in mean viability", 
       x = "")
Fig6d <- p6 + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), 
        plot.margin = margin(t=1, unit = "cm"))
Fig6d

#plot OTX015 and TW-37
theme_set(theme_bw() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(color="black", angle = 45, hjust=1, vjust=1, size=8), 
  axis.text.y = element_text(color="black", size=8), 
  axis.title = element_text(size=8),
  legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), 
  legend.key.size = unit(0.5, 'cm'),
  #axis.line = element_line(size = 0.5)
))

Fig6e <- differences_long_allconc |> 
  mutate(drug = case_when(drug == "OTX015" ~ "OTX015 (BET)",
                          drug == "TW.37" ~ "TW-37 (BCL2/MCL1)",
                          TRUE ~ drug)) |> 
  pivot_longer(cols = c(sample1_value, sample2_value), names_to = "timepoint", values_to = "viability") |> 
  filter(drug %in% c("OTX015 (BET)","TW-37 (BCL2/MCL1)"), status == "start between", therapy == "targeted") |> 
  ggplot(aes(x=rev(factor(concentration)), y=viability, color=timepoint))+
  geom_boxplot(alpha = 0.3, width = 0.5, position=position_dodge(0.75), outlier.shape = NA)+
  geom_point(aes(color = timepoint), alpha = 0.25, 
           position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2))+ 
  labs(y = "Viability", color= "Sample")+ 
  facet_wrap(~drug, ncol=2) + 
  stat_compare_means(aes(group = timepoint), 
                     method = "t.test",
                     paired = TRUE,
                     label = "p.signif",
                     label.y = 1.2, 
                     size=3)+
  labs(x = bquote("Concentration ("*mu*"M)"))+
  scale_y_continuous(breaks=seq(0,1.0, 0.2), limits=c(0,1.38)) +
  scale_color_manual(labels = c("sample1_value" = "First", 
                                  "sample2_value" = "Second"),
                       values = c("grey30", "#6baed6"))+
  scale_x_discrete(labels = c("0.016", "0.08", "0.4", "2", "10"))+
  guides(color = guide_legend(override.aes = list(label = c(""))))+
  theme(strip.text.y = element_text(size = 0), 
        strip.text.x = element_text(size=8),
        axis.text = element_text(size =8, color="black"), axis.title = element_text(size =8, color="black"),
        panel.spacing = unit(0.15, "lines"), 
        legend.title = element_text(size=8, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(color = "black", linewidth = 1, fill = "grey90"),        
        plot.margin = margin(l = 0.5, r = 1, unit = "cm"))
Fig6e <- Fig6e + theme(plot.margin = margin(t=1, b=0.5, r=1, unit = "cm"))
Fig6e

#plot BCR inhibitors
plot_2 <- differences_long_allconc |> 
  pivot_longer(cols = c(sample1_value, sample2_value), names_to = "timepoint", values_to = "viability") |> 
  filter(drug %in% c("Dasatinib","Ganetespib", "Idelalisib", "AZD8055", "Cobimetinib", "Ibrutinib"), status == "start between", therapy == "targeted") |> 
  ggplot(aes(x=rev(factor(concentration)), y=viability, color=timepoint))+
  geom_boxplot(alpha = 0.3, width = 0.5, position=position_dodge(0.75), outlier.shape = NA)+
  geom_point(aes(color = timepoint), alpha = 0.25, 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2))+ 
  labs(y = "Viability", color= "Sample")+ 
  facet_wrap(~drug) + 
  stat_compare_means(aes(group = timepoint), 
                     method = "t.test",
                     paired = TRUE,
                     label = "p.signif",
                     label.y = 1.2, 
                     size=3)+
  labs(x = bquote("Concentration ("*mu*"M)"))+
  scale_y_continuous(breaks=seq(0,1.0, 0.2), limits=c(0,1.38)) +
  scale_color_discrete(labels = c("sample1_value" = "First", 
                                  "sample2_value" = "Second"))+
  scale_x_discrete(labels = c("0.016", "0.08", "0.4", "2", "10"))+
  guides(color = guide_legend(override.aes = list(label = c(""))))
plot_2

#max relative difference
test <- df_ttest_allconc |> 
  filter(therapy == "targeted", status == "start between") |> 
  group_by(drug, concentration) |> 
  summarize(med_s1 = median(sample1_value),
                            med_s2=median(sample2_value),
                            abs_diff=med_s2-med_s1,
                            rel_diff=(med_s2-med_s1)/med_s2*100) |> 
  arrange(desc(rel_diff))

#plot MCL1 and BCL2 inhibitors
plot_2 <- differences_long_allconc |> 
  pivot_longer(cols = c(sample1_value, sample2_value), names_to = "timepoint", values_to = "viability") |> 
  filter(drug %in% c("Venetoclax","A.1210477"), status == "start between", therapy == "targeted") |> 
  ggplot(aes(x=rev(factor(concentration)), y=viability, color=timepoint))+
  geom_boxplot(alpha = 0.3, width = 0.5, position=position_dodge(0.75), outlier.shape = NA)+
  geom_point(aes(color = timepoint), alpha = 0.25, 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2))+ 
  labs(y = "Viability", color= "Sample")+ 
  facet_wrap(~drug) + 
  stat_compare_means(aes(group = timepoint), 
                     method = "t.test",
                     paired = TRUE,
                     label = "p.signif",
                     label.y = 1.2, 
                     size=3)+
  labs(x = bquote("Concentration ("*mu*"M)"))+
  scale_y_continuous(breaks=seq(0,1.0, 0.2), limits=c(0,1.38)) +
  scale_color_discrete(labels = c("sample1_value" = "First", 
                                  "sample2_value" = "Second"))+
  scale_x_discrete(labels = c("0.016", "0.08", "0.4", "2", "10"))+
  guides(color = guide_legend(override.aes = list(label = c(""))))
plot_2
