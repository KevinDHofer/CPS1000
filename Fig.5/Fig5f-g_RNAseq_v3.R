#create combined plots
library(cowplot)

#run RNAseq_CPS_and_Knisbacher_v2

Fig5f <- plot_fgsea_h +
  theme(axis.text.y = element_markdown(color="black", size=8),
        plot.margin = margin(t=2, l=1, b=1, r=2, unit = "cm"))
Fig5f


# Connect to the Ensembl database
library(stringr)

names(DEres_knis)[1] <- "gene_id"
names(DEres_knis)[2] <- "gene_symbol"

DEres_knis$gene_id <- str_remove(DEres_knis$gene_id, "\\.[^.]*$")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get gene symbols for ENSEMBL IDs
gene_info_knis <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
                   filters = "ensembl_gene_id",
                   values = DEres_knis$gene_id,
                   mart = ensembl)

gene_info_cps <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
                        filters = "ensembl_gene_id",
                        values = DEres_cps$gene_id,
                        mart = ensembl)

# Merge with original data
DEres_knis_anno <- merge(DEres_knis, gene_info_knis, 
                             by.x = "gene_id", 
                             by.y = "ensembl_gene_id")

DEres_cps_anno <- merge(DEres_cps, gene_info_cps, 
                         by.x = "gene_id", 
                         by.y = "ensembl_gene_id")

DEres_combined <- rbind(cbind(DEres_knis_anno, source=rep("Knisbacher", nrow(DEres_knis_anno))), 
                        cbind(DEres_cps_anno, source=rep("CPS1000", nrow(DEres_cps_anno))))

#remove genes correlated with CD3E (>0.8) & remove hemoglobins
strong_positive_CD3D <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/R notebooks/strong_positive_CD3D.csv")
remove_gene <- c(strong_positive_CD3D$Gene, "HBB", "HBA2")

DEres_combined <- DEres_combined[!(DEres_combined$gene_symbol %in% remove_gene),]

#annotate genesets
gene_sets_k <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets_h <- msigdbr(species = "Homo sapiens", collection = "H")

# Create a simple mapping between gene symbols and hallmark pathways
gene_to_pathway_k <- gene_sets_k %>%
  dplyr::select(gene_symbol, gs_name)

gene_to_pathway_h <- gene_sets_h %>%
  dplyr::select(gene_symbol, gs_name)

# Add  annotations to df
annotated_genes_k <- DEres_combined %>%
  merge(gene_to_pathway_k, by = "gene_symbol")%>%
  group_by(gene_id) %>%
  summarize(pathways = paste(gs_name, collapse = "; "))

annotated_genes_h <- DEres_combined %>%
  merge(gene_to_pathway_h, by = "gene_symbol")%>%
  group_by(gene_id) %>%
  summarize(pathways = paste(gs_name, collapse = "; "))

# Merge back to original dataframe
DEres_combined_geneset <- DEres_combined %>%
  left_join(annotated_genes_k, by = "gene_id")

annotated_genes_h <- annotated_genes_h |> 
  mutate(pathways_h = pathways)

DEres_combined_geneset <- merge(DEres_combined_geneset, annotated_genes_h[,-2], by = "gene_id", all = TRUE)

#write.csv(DEres_combined_geneset, "DEres_combined_geneset_Tcell_corrected.csv")
#DEres_combined_geneset <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DEres_combined_geneset_Tcell_corrected.csv")

#create common volcano
DEres_combined_geneset$sig <- ifelse((DEres_combined_geneset$padj < 0.1 & abs(DEres_combined_geneset$log2FoldChange) > 2), TRUE, FALSE)
DEres_combined_geneset$source_sig <- ifelse(DEres_combined_geneset$sig == TRUE, DEres_combined_geneset$source, "non_sig")

mycolors <- setNames(c("lightgrey", "lightblue", "red"), unique(DEres_combined_geneset$source_sig))
mycolors

set.seed(123)
volc <- DEres_combined_geneset |> 
  filter(gene_biotype == "protein_coding") |> 
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color=source_sig), size=0.75, alpha = 0.7) +
  #geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray") +
  #geom_vline(xintercept = c(-2,2), linetype = "dashed", color = "gray")+
  geom_label_repel(
    aes(label = ifelse(
      (str_detect(pathways, "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION|KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION|KEGG_HEMATOPOIETIC_CELL_LINEAGE") | 
         str_detect(pathways_h, "HALLMARK_TNFA_SIGNALING_VIA_NFKB|HALLMARK_INFLAMMATORY_RESONSE|HALLMARK_ALLOGRAFT_REJECTION")) & 
        source_sig != "non_sig", 
      gene_symbol,
      ""
    ), fill=source_sig, segment.color=source_sig),
    color="white",
    box.padding = 2, 
    max.overlaps = Inf,
    size=3
  )+
  theme_classic() +
  labs(
    title = "",
    x = expression(log[2]~FC),
    y = expression(-log[10]*italic(P)),
    color = "Source"
  ) +
  scale_color_manual(values = mycolors, 
                     labels=c("CPS1000", "CLL-map Portal", "not significant"))+
  scale_fill_manual(values = mycolors, 
                    labels=c("CPS1000", "CLL-map Portal", "not significant"),
                    guide = "none")+
  scale_discrete_manual("segment.color", values = mycolors, guide = "none")+
  annotate("text", x=-5, y=8, label = "Down in DDX3X\nmutants", size = 3)+
  annotate("text", x=5, y=8, label = "Up in DDX3X\nmutants", size = 3)+
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(size=8), 
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'))


Fig5g <- volc + theme(plot.margin = margin(t=2, l=1, b=1.5, r=1, unit = "cm"))
Fig5g
