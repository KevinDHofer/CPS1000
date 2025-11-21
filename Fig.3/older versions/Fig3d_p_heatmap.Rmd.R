library(tidyverse)
library(ggplot2)
library(dendextend)
library(DescTools)
library(ggh4x)
library(ggpattern)

tTest <- function(val, fac) {
  res <- t.test(val ~ fac, var.equal = TRUE)
  data.frame(p = res$p.value, 
             mean1 = res$estimate[[1]],
             mean2 = res$estimate[[2]],
             diff = res$estimate[[2]] - res$estimate[[1]])
}

#load data
usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/patmeta_210324.RData")

#load table with disease entities  
diseases <- usedSamples_all |> 
  arrange(patientID) |> 
  dplyr::select (patientID, diagnosis) 
diseases

# merge the tables
merged_table <- cbind(diseases, viability_table) |> 
  relocate(diagnosis, .before = X10058.F4_1)
merged_table

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:318), names_to = "drug", values_to = "viability")
merged_table_long

#separate M-CLL and U-CLL
merged_table_long <- mutate(merged_table_long, IGHV = patMeta[match(patientID, patMeta$Patient.ID),]$IGHV.status) %>%
  filter(diagnosis != "CLL" | !is.na(IGHV)) %>%
  mutate(diagnosis = ifelse(diagnosis == "CLL", paste0(IGHV, "-", diagnosis), diagnosis))

# only select diseases with at least five cases
diagSelected <- c("U-CLL", "M-CLL", "AML", "B-ALL", "T-ALL", "B-PLL", "T-PLL", "MCL")

merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
merged_table_long

#create means instead of of auc
viabTab.mean <- merged_table_long |> 
  filter(diagnosis %in% diagSelected) |> 
  group_by(patientID, diagnosis, drug) |> 
    summarize(viab = mean(viability))

#create p-values
allDiag <- unique(viabTab.mean$diagnosis)
pTabAll <- lapply(allDiag, function(diagBase) {
  #calculate p value matrix
  pTab <- lapply(allDiag[allDiag != diagBase], function(diagCompare) {
    testTab <- filter(viabTab.mean, diagnosis %in% c(diagBase, diagCompare)) %>%
      mutate(diagnosis = factor(diagnosis, levels = c(diagBase,diagCompare))) 
    res <- group_by(testTab, drug) %>% do(tTest(.$viab, .$diagnosis)) %>%
      mutate(diagCompare = diagCompare, diagBase = diagBase)
    res
  }) %>% bind_rows() %>% ungroup() %>% mutate(p.adj = p.adjust(p, method = "BH"))
}) %>% bind_rows()

#prepare data table for plot 
plotTab <- pTabAll %>% mutate(p = ifelse(p < 1e-12, 1e-12, p)) %>%
  mutate(pSign = -log10(p)*sign(diff)) %>%
  mutate(pSign = ifelse(p.adj < 0.1, pSign, 0)) %>%
  ungroup()

#Order drugs by their association similarity (hierarchical clustering)
pMat <- mutate(plotTab, diagPair = paste0(diagBase,"_",diagCompare)) %>%
  dplyr::select(drug, pSign, diagPair) %>%
  spread(key = diagPair, value = pSign) %>% data.frame() %>%
  column_to_rownames("drug") %>% as.matrix()

# Cut
nClust = 9
hcRes <- hclust(dist(pMat), method = "ward.D2")
drugOrder <- rownames(pMat)[hcRes$order]
#clustNum <- rev(c(15, 5, 9, 4, 3, 8, 6, 6, 7))
#clustIndex <- rep(seq(length(clustNum)), clustNum)
clustMap <- cutree(hcRes,nClust)

diagOrder <- c("U-CLL", "M-CLL","MCL", "B-PLL", "T-PLL", "B-ALL", "T-ALL","AML")
plotTab <- mutate(plotTab, drug = factor(drug, levels = drugOrder), 
                  diagBase = factor(diagBase, levels= diagOrder), 
                  diagCompare = factor(diagCompare, levels = diagOrder)) %>%
  arrange(drug) %>% mutate(drugClust = paste0("",nClust-clustMap[as.character(drug)]+1))

# Colour strips in y-direction
annotation_pathway <- c("lightgrey", "lightgrey", "#80B1D3", "#BEBADA", "lightgrey", "#FDB462", "#BC80BD", "lightgrey", "lightgrey")
strip <- strip_themed(background_y = elem_list_rect(fill = annotation_pathway))

#plot
Fig3d <- ggplot(plotTab, aes(x=diagCompare, y = drug, fill = pSign)) + 
  geom_tile(size = 0.2, color = "white") +
  facet_grid2(rows = vars(drugClust), cols = vars(diagBase), scales = "free", space = "free_y", strip = strip) + 
  scale_fill_gradient2(high = "blue", mid = "white", low = "red", midpoint = 0,
                       name = "-log10(p) \n with direction") +
  #ggtitle(sprintf("Relative drug sensitivity to %s", diagBase)) +
  theme_bw() +   theme(strip.text.y = element_text(size = 0), 
                       strip.text.x = element_text(size=8),
                       strip.background.y = element_rect(color=NA),
                       axis.text.y = element_text(size =8, color= "black"),
                       axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size=8, color = "black"),
                       panel.spacing = unit(0.15, "lines"), 
                       legend.title = element_text(size=8, face="bold")) +
  ylab("") + xlab("")

Fig3d

#ggsave("pvalue_heatmap_all_Vs_all_hierachical_clustering.pdf")
