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

#read the TSV file 
dt <- fread("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/cllmap_rnaseq_counts_full.tsv") 

#import bed file
bed <- as.data.frame(read.table("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/knisbacher2022.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#import excel file
df_meta <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                                                sheet = "Supplementary Table 1")

df_WES <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                      sheet = "Supplementary Table 3a")

df_WGS <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                     sheet = "Supplementary Table 3b")

df_driver <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                     sheet = "Supplementary Table 4a")

df_GSEA <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                     sheet = "Supplementary Table 6")

df_IGH <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                      sheet = "Supplementary Table 8b")

df_epi <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                     sheet = "Supplementary Table 12a")

df_OS_FFS <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DDX3X_BTKi/literature/Knisbacher 2022_DDX3X_pts/41588_2022_1140_MOESM3_ESM.xlsx", 
                        sheet = "Supplementary Table 14b")

#filter DDX3X patients: 33 pat
my.names <- df_WGS[2,]
colnames(df_WGS) <- my.names
df_WGS <-df_WGS[3:nrow(df_WGS),1:22]

#select DDX3X
DDX3X <- df_WGS |> filter(Hugo_Symbol == "DDX3X")

#number of DDX3X mutated patients with RNAseq data: 18 pat
ncol(dplyr::select(dt, any_of(DDX3X$participant_id)))

#merge with patMeta
meta.names <- df_meta[2,]
colnames(df_meta) <- meta.names
df_meta2 <-df_meta[4:nrow(df_meta),1:9]

#filter DDX3X: 33 pat
DDX3X_merged <- merge(DDX3X, df_meta2, by.x="participant_id", by.y="PARTICIPANT ID")

#DDX3X_females <- DDX3X_merged |> filter(SEX == "F")
#DDX3X_list_females <- DDX3X_females$participant_id

DDX3X_list <- DDX3X_merged$participant_id

#select U-CLL only: 547 pat (males 386 pat)
df_ucll <- df_meta2 |> 
  filter(`IGHV MUTATION STATUS` == "unmutated", 
         #SEX == "M"
         ) 

ucll_list <- df_ucll$"PARTICIPANT ID"

dt_ucll <- dt |> 
  dplyr::select(Name, Description, any_of(ucll_list))

#number of U-CLL patients with RNAseq data: 343 pat, males 247 pat
ncol(dt_ucll[,-(1:2)])

#remove IGHV and IGLV
#removed_IGHV <- dt_ucll |> filter((str_detect(Description, "IGHV")))
#removed_IGLV <- dt_ucll |> filter((str_detect(Description, "IGLV")))
#dt_ucll <- dt_ucll |> filter(!(Description %in% c(removed_IGHV$Description, removed_IGLV$Description)))

#prepare matrix for DESeq2
countData <- as.matrix(dt_ucll[,-(1:2)])
rownames(countData) <- dt_ucll$Name

X <- c()
for (i in colnames(countData)) {
  ifelse(i %in% DDX3X_list, X <- append(X, "mutated"), X <- append(X,"unmutated"))
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
  scale_color_manual(values=c("blue", "lightgrey")) +
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
  scale_color_manual(values=c("blue", "lightgrey")) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()

#PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition", pcsToUse = c(1,2))
plotPCA(vsd, intgroup="condition", pcsToUse = c(3,4))

#add gene names (without LFC shrinkage)
#resLFC_tidy <- as.data.frame(resLFC) %>% rownames_to_column("gene_id") %>%
  #as_tibble() |> 
  #mutate(stat = log2FoldChange / lfcSE)

DEres <- results(dds, tidy=TRUE) 
DEres <- merge(dt_ucll[,1:2], DEres, by.x="Name", by.y="row")

#List of all differentially expressed genes (10% FDR cut-off, adj < 0.1)
sigDE <- DEres %>% dplyr::select(Description, Name, log2FoldChange, pvalue, padj) %>%
  filter(padj <= 0.1) %>% arrange(pvalue) 

#Positive log2FoldChange indicates that the gene is up-regulated in DDX3X mutated samples and vice versa
sigDE %>%
  mutate_if(is.numeric, formatC, digits=2, format= "e") %>%
  DT::datatable()

#drop genes with padj = NA
DEres <- DEres |> 
  drop_na(padj)

DEres_knis <- DEres

# Prepare ranked gene list from DESeq2 results
# Sort by the stat column (log2FoldChange divided by standard error)
ranked_genes <- DEres %>%
  dplyr::arrange(desc(stat)) %>%
  dplyr::select(Name, Description, stat)

# Check for duplicates
duplicate_genes <- ranked_genes$Description[duplicated(ranked_genes$Description)]

# Remove duplicates by keeping the entry with the highest absolute stat value
ranked_genes_unique <- ranked_genes %>%
  group_by(Description) %>%
  slice_max(abs(stat), n = 1) %>%
  ungroup()

#remove genes with high correlation (>0.8) with T cell genes, such as CD3E
strong_positive_CD3E <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/R notebooks/strong_positive_CD3D.csv")
gene_remove <- strong_positive_CD3E$Gene
ranked_genes_unique <- ranked_genes_unique[!(ranked_genes_unique$Description %in% strong_positive_CD3E$Gene),]

#filter for protein-coding genes
library(stringr)
library(biomaRt)

names(ranked_genes_unique)[1] <- "gene_id"
names(ranked_genes_unique)[2] <- "gene_symbol"

ranked_genes_unique$gene_id <- str_remove(ranked_genes_unique$gene_id, "\\.[^.]*$")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get gene symbols for ENSEMBL IDs
gene_info_knis <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
                        filters = "ensembl_gene_id",
                        values = ranked_genes_unique$gene_id,
                        mart = ensembl)

# Merge with original data
ranked_genes_filtered <- merge(ranked_genes_unique, gene_info_knis, 
                         by.x = "gene_id", 
                         by.y = "ensembl_gene_id") |> 
  filter(gene_biotype == "protein_coding") |> 
  arrange(abs(stat))

# Create a named vector without duplicates
gene_ranks <- setNames(ranked_genes_filtered$stat, ranked_genes_filtered$gene_symbol)

gene_ranks_knisbacher <- gene_ranks
#write.csv(gene_ranks_knisbacher, "~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/gene_ranks_knisbacher_v2.csv")

####prepare for plotting

# Get gene sets
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run fgsea with the deduplicated data
fgsea_hallmark <- fgsea(
  pathways = hallmark_gene_sets,
  stats = gene_ranks_knisbacher,
  minSize = 15,
  maxSize = 500
)

# View results
fgsea_hallmark <- fgsea_hallmark %>% 
  arrange(pval)
head(fgsea_hallmark)

# Plot top enriched pathways
topPathways_hallmark <- unlist(fgsea_hallmark[head(order(fgsea_hallmark$pval), n=10), "pathway"])
plotGseaTable(hallmark_gene_sets[topPathways_hallmark], gene_ranks, fgsea_hallmark, 
              gseaParam = 0.5)

fgsea_hallmark_df_knis <- as.data.frame(fgsea_hallmark)
fgsea_hallmark_df_knis$pathway <- gsub("HALLMARK_", "", fgsea_hallmark_df_knis$pathway)

#saveRDS(fgsea_hallmark_df_knis, "fgsea_hallmark_df_knis_v2.rds")

topPathways_hallmark_clean <- fgsea_hallmark_df_knis$pathway[1:10]

#ggplot2
plot_fgsea_h <- fgsea_hallmark_df_knis |> filter(
  #pathway %in% topPathways_hallmark_clean, 
  padj <0.01) |> 
  ggplot(aes(x=reorder(pathway, NES), y=NES, color= (-log10(padj)), size=size)) +
  geom_point(alpha=0.75) +
  coord_flip() +
  labs(
    y = "Normalized Enrichment Score (NES)",
    x = "",
    color = "-log10(FDR)",
    size = "Leading edge \nsize",
    title = "Hallmark (Knisbacher, FDR<0.1)"
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
fgsea_kegg_knis <- fgsea(
  pathways = KEGG_gene_sets,
  stats = gene_ranks_knisbacher,
  minSize = 15,
  maxSize = 500
)

# View results
fgsea_kegg_knis <- fgsea_kegg_knis %>% 
  arrange(pval)
head(fgsea_kegg_knis)

# Plot top enriched pathways
fgsea_kegg_df_knis <- as.data.frame(fgsea_kegg_knis)
fgsea_kegg_df_knis$pathway <- gsub("KEGG_", "", fgsea_kegg_df_knis$pathway)
#saveRDS(fgsea_kegg_df_knis, "fgsea_kegg_df_knis_v2.rds")

topPathways_kegg_clean_knis <- fgsea_kegg_df_knis$pathway[1:10]

# plot top10 for GO
plot_fgsea_kegg <- fgsea_kegg_df_knis |> filter(
  #pathway %in% topPathways_kegg_clean_knis,
  padj <0.01) |> 
  ggplot(aes(x=reorder(pathway, NES), y=NES, color= (-log10(padj)))) +
  geom_point(aes(size=size), alpha=0.75) +
  coord_flip() +
  labs(
    y = "Normalized Enrichment Score (NES)",
    x = "",
    color = "-log10(FDR)",
    size = "Leading edge \nsize",
    title = "Reactome (Knisbacher, FDR<0.1)"
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

#same for Reactome Pathways
KEGG_gene_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run fgsea
fgsea_kegg_knis <- fgsea(
  pathways = KEGG_gene_sets,
  stats = gene_ranks_knisbacher,
  minSize = 15,
  maxSize = 500
)

# View results
fgsea_kegg_knis <- fgsea_kegg_knis %>% 
  arrange(pval)
head(fgsea_kegg_knis)
