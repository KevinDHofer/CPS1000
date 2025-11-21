library(data.table) 
library(tidyverse)
library(xml2)
library(readxl)
library(DESeq2)
library(limma)
library(apeglm)
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/patmeta_191004.RData")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/ddsrna_180717.RData")

#Extract count matrix from dds object
count_matrix <- counts(dds)
count_df <- as.data.frame(count_matrix)

#add gene information from the rowData
gene_info <- as.data.frame(rowData(dds))
count_df_with_info <- cbind(gene_info$symbol, count_df)
colnames(count_df_with_info)[1] <- "gene_symbol"

#select DDX3X patients: 8 pat
DDX3X <- patMeta |> filter(DDX3X == 1) |> 
  drop_na(date.of.diagnosis)
DDX3X_pat_list <- DDX3X$Patient.ID

#number of DDX3X mutated patients with RNAseq data: 4 pat
ncol(dplyr::select(count_df_with_info, any_of(DDX3X$Patient.ID)))
DDX3X <- sort(colnames(count_df_with_info)[colnames(count_df_with_info) %in% DDX3X_pat_list])

#select U-CLL only: 165 pat
df_ucll <- patMeta |> 
  filter(IGHV.status == "U")
ucll_pat_list <- df_ucll$Patient.ID

#number of U-CLL patients with RNAseq data: 98 pat
ncol(dplyr::select(count_df_with_info, any_of(ucll_pat_list)))
ucll <- sort(colnames(count_df_with_info)[colnames(count_df_with_info) %in% ucll_pat_list])

dt_ucll <- count_df_with_info |> 
  dplyr::select(gene_symbol, any_of(ucll))

#prepare matrix for DESeq2
countData <- as.matrix(dt_ucll[,-1])
#rownames(countData) <- dt_ucll$gene_symbol

X <- c()
for (i in colnames(countData)) {
  ifelse(i %in% DDX3X, X <- append(X, "mutated"), X <- append(X,"unmutated"))
}
condition <- factor(X)

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = DataFrame(condition), 
                              design= ~ condition)

dds$condition <- relevel(dds$condition, ref = "unmutated")

#filter low count genes
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#run DESeq2
dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$pvalue),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

#Log fold change shrinkage
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_mutated_vs_unmutated", type="apeglm")
resLFC

# quality check
# distribution of adj. p-values
hist(res$pvalue, col="lightblue", breaks = 50, main = "p-value distribution")

#MA plot -> not working(?)
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

#alternative plotting
resDF <- as.data.frame(res)
resDF$sig <- ifelse(resDF$padj < 0.1, "FDR<0.1", "Not Sig")
ggplot(resDF, aes(x=log10(baseMean), y=log2FoldChange, color=sig)) +
  geom_point(alpha=0.4, size=0.75) +
  scale_color_manual(values=c("blue", "black")) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()+
  ylim(c(-2,2))

#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

#plot for LFC
resLFC_DF <- as.data.frame(resLFC)
resLFC_DF$sig <- ifelse(resLFC_DF$padj < 0.1, "FDR<0.1", "Not Sig")
ggplot(resLFC_DF, aes(x=log10(baseMean), y=log2FoldChange, color=sig)) +
  geom_point(alpha=0.4, size=0.75) +
  scale_color_manual(values=c("blue", "black")) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()

#PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition", pcsToUse = c(1,2))
plotPCA(vsd, intgroup="condition", pcsToUse = c(3,4))

#drop genes with padj = NA
DEres <- results(dds, tidy=TRUE)
DEres <- DEres |> 
  drop_na(padj, row, stat)

#List of all differentially expressed genes (10% FDR cut-off, adj < 0.1)
sigDE <- DEres %>% dplyr::select(row, log2FoldChange, pvalue, padj) %>%
  filter(padj <= 0.1) %>% arrange(pvalue) 

#Positive log2FoldChange indicates that the gene is up-regulated in DDX3X mutated samples and vice versa
sigDE %>%
  mutate_if(is.numeric, formatC, digits=2, format= "e") %>%
  DT::datatable()

# Prepare ranked gene list from DESeq2 results
# Sort by the stat column (log2FoldChange divided by standard error)
ranked_genes <- DEres %>%
  dplyr::arrange(desc(stat)) %>%
  dplyr::select(row, stat)

#add gene names
dt_ucll2 <- mutate(dt_ucll, gene_id = rownames(dt_ucll), .before = "gene_symbol") |> 
  drop_na(gene_symbol)
DEres <- merge(dt_ucll2[,1:2], DEres, by.x="gene_id", by.y="row")

DEres_cps <- DEres

# Check for duplicates
duplicate_genes <- ranked_genes$gene_symbol[duplicated(ranked_genes$gene_symbol)]

# Remove duplicates by keeping the entry with the highest absolute stat value
ranked_genes_unique <- ranked_genes %>%
  group_by(row) %>%
  slice_max(abs(stat), n = 1) %>%
  ungroup()

ranked_genes_unique <- merge(DEres[1:2], ranked_genes_unique, by.x="gene_id", by.y = "row", .before = 1)
ranked_genes_unique <- ranked_genes_unique |> 
  drop_na(stat, gene_symbol) |> 
  filter(gene_symbol != "")
ranked_genes_unique <- ranked_genes_unique |> 
  group_by(gene_symbol) %>%
  slice_max(abs(stat), n = 1) %>%
  ungroup()

#remove genes with high correlation (>0.8) with T cell genes (such as CD3E) & remove hemoglobins
strong_positive_CD3E <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/R notebooks/strong_positive_CD3D.csv")
gene_remove <- strong_positive_CD3E$Gene
ranked_genes_unique <- ranked_genes_unique[!(ranked_genes_unique$gene_symbol %in% strong_positive_CD3E$Gene),] |> 
  filter(gene_symbol != "HBB|HBA2")

#filter for protein-coding genes
library(stringr)
library(biomaRt)

names(ranked_genes_unique)[1] <- "gene_id"
names(ranked_genes_unique)[2] <- "gene_symbol"

ranked_genes_unique$gene_id <- str_remove(ranked_genes_unique$gene_id, "\\.[^.]*$")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get gene symbols for ENSEMBL IDs
gene_info_cps <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
                        filters = "ensembl_gene_id",
                        values = ranked_genes_unique$gene_id,
                        mart = ensembl)

# Merge with original data
ranked_genes_filtered <- merge(ranked_genes_unique, gene_info_cps, 
                               by.x = "gene_id", 
                               by.y = "ensembl_gene_id") |> 
  filter(gene_biotype == "protein_coding") |> 
  arrange(abs(stat))

# Create a named vector without duplicates
gene_ranks <- setNames(ranked_genes_filtered$stat, ranked_genes_filtered$gene_symbol)

gene_ranks_CPS1000 <- gene_ranks
#write.csv(gene_ranks_CPS1000, "~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/gene_ranks_CPS1000_v2.csv")

#### prepare for plotting

# Get gene sets
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run fgsea with the deduplicated data
fgsea_hallmark <- fgsea(
  pathways = hallmark_gene_sets,
  stats = gene_ranks_CPS1000,
  minSize = 15,
  maxSize = 500
)

# View results
fgsea_hallmark <- fgsea_hallmark %>% 
  arrange(pval)
head(fgsea_hallmark)

#setwd("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/R notebooks")
#write_csv(fgsea_hallmark, "fgsea_hallmark_baseline_v2.csv")

# Plot top enriched pathways
topPathways_hallmark <- unlist(fgsea_hallmark[head(order(fgsea_hallmark$pval), n=10), "pathway"])
plotGseaTable(hallmark_gene_sets[topPathways_hallmark], gene_ranks, fgsea_hallmark, 
              gseaParam = 0.5)

fgsea_hallmark_df_cps <- as.data.frame(fgsea_hallmark)
fgsea_hallmark_df_cps$pathway <- gsub("HALLMARK_", "", fgsea_hallmark_df_cps$pathway)
#saveRDS(fgsea_hallmark_df_cps, "fgsea_hallmark_df_cps_v2.rds")

topPathways_hallmark_clean <- fgsea_hallmark_df_cps$pathway[1:10]

#ggplot2
plot_fgsea_h <- fgsea_hallmark_df_cps |> filter(padj <0.1) |> 
  ggplot(aes(x=reorder(pathway, NES), y=NES, color= (-log10(padj)), size=size)) +
  geom_point(alpha=0.75) +
  coord_flip() +
  labs(
    y = "Normalized Enrichment Score (NES)",
    x = "",
    color = "-log10(FDR)",
    size = "Leading edge \nsize",
    title = "Hallmark (CPS, FDR<0.1)"
  ) +
  theme_bw() +
  scale_color_continuous(low = "grey80", high = "red")+
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(size=8), 
        plot.title = element_text(size=10, face="bold", hjust = 0.5),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'), 
  )+
  scale_size_continuous(range = c(1, 4))
plot_fgsea_h

#same for KEGG Pathways
KEGG_gene_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run fgsea
fgsea_KEGG_cps <- fgsea(
  pathways = KEGG_gene_sets,
  stats = gene_ranks_CPS1000,
  minSize = 15,
  maxSize = 500
)

# View results
fgsea_KEGG_cps <- fgsea_KEGG_cps %>% 
  arrange(pval)
head(fgsea_KEGG_cps)

# Plot top enriched pathways
fgsea_KEGG_df_cps <- as.data.frame(fgsea_KEGG_cps)
fgsea_KEGG_df_cps$pathway <- gsub("REACTOME_", "", fgsea_KEGG_df_cps$pathway)
#saveRDS(fgsea_KEGG_df_cps, "fgsea_KEGG_df_cps_v2.rds")

topPathways_KEGG_clean_cps <- fgsea_KEGG_df_cps$pathway[1:20]

# plot top10 for GO
plot_fgsea_kegg <- fgsea_KEGG_df_cps |> filter(
  pathway %in% topPathways_KEGG_clean_cps,
  #padj<0.1
  ) |> 
  ggplot(aes(x=reorder(pathway, NES), y=NES, color= (-log10(padj)))) +
  geom_point(aes(size=size), alpha=0.75) +
  coord_flip() +
  labs(
    y = "Normalized Enrichment Score (NES)",
    x = "",
    color = "-log10(FDR)",
    size = "Leading edge \nsize",
    title = "REACTOME (CPS, top20)"
  ) +
  theme_bw() +
  scale_color_continuous(low = "grey80", high = "red")+
  theme(axis.text = element_text(color="black", size=8), 
        axis.title = element_text(size=8), 
        plot.title = element_text(size=10, face="bold", hjust = 0.5),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'), 
  )+
  scale_size_continuous(range = c(1, 4))
plot_fgsea_kegg

#plot BCR pathway
pathway_name <- "REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR"
plotEnrichment(KEGG_gene_sets[[pathway_name]], 
               gene_ranks_CPS1000) +
  theme_bw() +
  labs(title = paste("RNA:", pathway_name),
       subtitle = paste("NES =", round(fgsea_KEGG_cps[pathway == pathway_name]$NES, 3),
                        ", p-adj =", round(fgsea_KEGG_cps[pathway == pathway_name]$padj, 3),
                        ", pval =", round(fgsea_KEGG_cps[pathway == pathway_name]$pval, 3)))


