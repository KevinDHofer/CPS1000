#create combined plots

#import hallmark and KEGG gene sets from CPS1000  and Knisbacher
fgsea_KEGG_df_cps <- readRDS("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.5/preparatory scripts for RNAseq/fgsea_KEGG_df_cps_v2.rds")
fgsea_hallmark_df_cps <- readRDS("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.5/preparatory scripts for RNAseq/fgsea_hallmark_df_cps_v2.rds")

fgsea_hallmark_df_knis <- readRDS("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.5/preparatory scripts for RNAseq/fgsea_hallmark_df_knis_v2.rds")
fgsea_kegg_df_knis <- readRDS("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/figures/Fig.5/preparatory scripts for RNAseq/fgsea_kegg_df_knis_v2.rds")

#hallmark: common gene sets
fgsea_hallmark_df_cps_filtered <- fgsea_hallmark_df_cps |> 
  filter(pathway %in% fgsea_hallmark_df_knis[1:10,]$pathway)

hallmarks_combined <- cbind(rbind(fgsea_hallmark_df_knis[1:10,], fgsea_hallmark_df_cps_filtered), 
                            Source=c(rep("Knisbacher", 10), rep("CPS1000",10)))

hallmarks_combined$Source <- factor(hallmarks_combined$Source, levels = c("CPS1000", "Knisbacher"))

#common top 10
top10_hallmark <- intersect(fgsea_hallmark_df_cps$pathway[1:10], 
                   fgsea_hallmark_df_knis$pathway[1:10])

#ggplot2
library(ggtext)
fgsea_hallmark_df_cps_plot <- fgsea_hallmark_df_cps |> 
  arrange(padj) |> 
  slice_head(n=10) |> 
  mutate(direction = as.character(ifelse(sign(NES) <0, "Down", "Up")),
         is_top10 = pathway %in% top10_hallmark,
         pathway = as.character(pathway),
         pathway_formatted = as.character(ifelse(is_top10, paste0("**", pathway, "**"), pathway)))

fgsea_hallmark_df_cps_plot$pathway_formatted <- factor(fgsea_hallmark_df_cps_plot$pathway_formatted, 
                                                  levels = sort_by(unique(fgsea_hallmark_df_cps_plot$pathway_formatted), 
                                                                   rev(fgsea_hallmark_df_cps_plot$pval)))


plot_fgsea_h <- fgsea_hallmark_df_cps_plot |> 
  ggplot(aes(x=-log10(pval), y=pathway_formatted, fill=direction)) +
  geom_col(alpha=0.75) +
  labs(
    y = "",
    x = expression(-log[10]*italic(P)),
    fill = "Direction",
    title = "Hallmark"
  ) +
  theme_bw() +
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(size=8), 
        plot.title = element_text(size=10, face="bold", hjust = 0.5),
        panel.border = element_rect(color = 'black', 
                                    fill = NA, 
                                    size = 1),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_discrete(expand = c(0.075,0.01))+
  scale_x_continuous(expand = c(0.025,0))+
  scale_fill_manual(values="lightcoral")
plot_fgsea_h

plot_fgsea_h <- plot_fgsea_h +
  theme(axis.text.y = element_markdown(color="black", size=8))
plot_fgsea_h

#plot KEGG pathways
#import KEGG gene sets from CPS1000 and Knisbacher
cps_kegg <- as.data.frame(fgsea_KEGG_df_cps) |> 
  arrange(pval)
cps_kegg <- cbind(cps_kegg, Source=rep("CPS1000",nrow(cps_kegg)))
cps_kegg$pathway <- gsub("REACTOME_", "", cps_kegg$pathway)

knis_kegg <- as.data.frame(fgsea_kegg_df_knis)
knis_kegg <- cbind(knis_kegg, Source=rep("Knisbacher",nrow(knis_kegg)))
knis_kegg$pathway <- gsub("REACTOME_", "", knis_kegg$pathway)

#kegg/reactome: common gene sets
cps_kegg_filtered <- cps_kegg |> 
  filter(pathway %in% knis_kegg[1:10,]$pathway)

fgsea_kegg_combined <- rbind(cps_kegg_filtered, knis_kegg[1:10,])

fgsea_kegg_combined <- fgsea_kegg_combined %>% 
  arrange(pval)
head(fgsea_kegg_combined)

fgsea_kegg_combined$pathway <- factor(fgsea_kegg_combined$pathway, levels = sort_by(unique(fgsea_kegg_combined$pathway), fgsea_kegg_combined$NES))

fgsea_kegg_combined <- fgsea_kegg_combined |> 
  arrange(desc(Source), NES) |> 
    mutate(order = row_number())
fgsea_kegg_combined$pathway <- factor(fgsea_kegg_combined$pathway, levels = unique(fgsea_kegg_combined$pathway))

# common top10 for KEGG
top10_kegg <- intersect(cps_kegg$pathway[1:10], 
                        knis_kegg$pathway[1:10])

# Define which pathways should be bold
cps_kegg_plot <- cps_kegg |> 
  arrange(padj) |> 
  slice_head(n=10) |> 
  mutate(direction = as.character(ifelse(sign(NES) <0, "Down", "Up")),
         is_top10 = pathway %in% top10_kegg,
         pathway = as.character(pathway),
         pathway_formatted = as.character(ifelse(is_top10, paste0("**", pathway, "**"), pathway)))

cps_kegg_plot$pathway_formatted <- factor(cps_kegg_plot$pathway_formatted, 
                                          levels = sort_by(unique(cps_kegg_plot$pathway_formatted), rev(cps_kegg_plot$pval)))

#plot
plot_fgsea_reactome <- cps_kegg_plot |> 
  ggplot(aes(x=-log10(pval), y=pathway_formatted, fill=direction)) +
  geom_col(alpha=0.75) +
  labs(
    y = "",
    x = "-log10(p value)",
    fill = "Direction",
    title = "Reactome"
  ) +
  theme_bw() +
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(size=8), 
        plot.title = element_text(size=10, face="bold", hjust = 0.5),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_discrete(expand = c(0.075,0.01))+
  scale_x_continuous(expand = c(0.025,0))
plot_fgsea_reactome

plot_fgsea_reactome <- plot_fgsea_reactome +
  theme(axis.text.y = element_markdown(color="black", size=8))
plot_fgsea_reactome

