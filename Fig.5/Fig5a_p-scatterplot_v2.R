library(tidyverse)
library(stringr)
library(ggplot2)

pTab <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/pTable_drug_VS_gene.csv", sep=";")

cols_to_convert <- (4:8)
pTab[cols_to_convert] <- lapply(pTab[cols_to_convert], function(x) {
  as.numeric(gsub(",", ".", x))
})

#Number of significant associations per gene (10% FDR)
plotTab1 <- filter(pTab, p.adj <= 0.1) %>% group_by(gene) %>%
  summarise(number = length(drug)) %>%
  arrange(desc(number)) %>% mutate(gene = factor(gene, levels = gene))

#P value scatter plot (IGHV is not shown)
# List a numuber of genes that should plot explicitly, based on plot above
top_genes <- as.vector(unlist(plotTab1[2:8,"gene"]))
showGenes <- sort(top_genes)
append(showGenes, "other")

#define colors for each gene
library(RColorBrewer)
colorList <- c(colorRampPalette(brewer.pal(7, 'Set2'))(7), 
               "grey80", "grey90")
  
names(colorList) = c(showGenes, "other", "Below 10% FDR")

colorList["del17p"] <- "#FC8D62" 
colorList["del11q"] <- "#E5C494"
colorList["trisomy12"] <- "#FFD92F"
colorList["DDX3X"] <- "#E78AC3"
colorList["TP53"] <- "#66C2A5"

#define shapes for viability
shapeList <- ifelse(names(colorList) %in% showGenes, 19,1)
names(shapeList) <- names(colorList)
shapeList["other"] <- 19



shapeList <- c(1, 25, 24)
names(shapeList) = c("Below 10% FDR", "Down", "Up")

#order drugs by their association similarity
pMat <- dplyr::select(pTab, drug, gene, p) %>% 
  filter(gene!="IGHV.status") %>%
  spread(key = gene, value = p) %>%
  column_to_rownames("drug") %>% as.matrix()

hc <- hclust(dist(pMat), method = "ward.D2")
drugOrder <- rownames(pMat)[hc$order]

#order drugs by max p
pval_list <- pTab |> 
  filter(gene != "IGHV.status") |> 
  group_by(drug) |> 
  summarize(pval = min(p)) |> 
  arrange(pval, descending = FALSE)

drugOrder2 <- pval_list$drug

#prepare plot table
plotTab <- pTab %>% filter(gene != "IGHV.status") %>%
  mutate(gene = ifelse(gene %in% showGenes, gene, "other")) %>%
  mutate(gene_FDR = ifelse(p.adj < 0.1, gene, "Below 10% FDR"),
         drug = factor(drug, levels = drugOrder2),
         gene = factor(gene, levels = names(colorList)),
         direction = ifelse(gene_FDR == "Below 10% FDR", "Below 10% FDR", ifelse(sign(diff)>0, "Down", "Up"))) 

Fig5a <- ggplot(plotTab, aes(x=drug, y = -log10(p), color = gene, shape = direction)) +
  geom_point(size =2, alpha = 0.8, stroke = 1) + scale_color_manual(values = colorList) +
  scale_shape_manual(values = shapeList, name = "Viability effect") + theme_bw() +
  labs(x="", y = expression(-log[10]*italic(P)), color = "Genetic aberration")+
  theme(axis.text.x = element_text(angle = 45, hjust =1, vjust=1, 
                                   color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), 
        axis.title.y = element_text(size=8), 
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8),
        panel.grid.major.x = element_line(size=.1, color="grey50", 
                                          linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = 'black', 
                                    fill = NA, 
                                    size = 1))+
  geom_hline(yintercept=-log10(0.004706046), linetype = "dashed", linewidth=0.3, color="black")+
  annotate("text", label="FDR 10%", size=3, x=60, y=3)
Fig5a <- Fig5a + theme(plot.margin = margin(t=1, l=1, unit = "cm"))
Fig5a


theme(legend.position = "none", 
      axis.text.x = element_text(angle = 0, color="black", size=8), 
      panel.grid = element_line(colour = "white"), 
      axis.text.y = element_text(color="black", size=8), 
      axis.title.y = element_text(size=8), 
      axis.title.x = element_text(size=8), 
      panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
      strip.background = element_rect(color = "black", linewidth = 1, fill = "grey90"),
      plot.title = element_text(hjust=0.5, size=8, face="bold"), 
      strip.text = element_text(size = 8), 
      plot.margin = margin(l = 0.5, r = 1, unit = "cm"))