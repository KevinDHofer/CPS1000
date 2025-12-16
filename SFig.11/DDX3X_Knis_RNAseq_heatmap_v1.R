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
library(ComplexHeatmap)
library(circlize)


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

#drop genes with padj = NA
DEres <- DEres |> 
  drop_na(padj)

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

#load common DE geneset table
DEres_combined_geneset <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DEres_combined_geneset.csv")

#remove genes with high correlation (>0.8) with T cell genes (such as CD3E) & remove hemoglobins
strong_positive_CD3E <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/R notebooks/strong_positive_CD3D.csv")
gene_remove <- c(strong_positive_CD3E$Gene, "HBB", "HBA2")
DEres_combined_geneset <- DEres_combined_geneset |> 
  filter(!gene_symbol %in% gene_remove)

#select gene sets
#TNF
list_knis_tnf <- DEres_combined_geneset |> 
  filter(gene_biotype == "protein_coding", 
         grepl("HALLMARK_TNFA_SIGNALING_VIA_NFKB", pathways_h),
         source == "Knisbacher", 
         abs(log2FoldChange) > 1
  ) |> 
  arrange(pvalue)

#mTORC
list_knis_mtor <- DEres_combined_geneset |> 
  filter(gene_biotype == "protein_coding", 
         grepl("HALLMARK_MTORC1_SIGNALING", pathways_h),
         source == "Knisbacher", 
         abs(log2FoldChange) > 1
  ) |> 
  arrange(pvalue)

#INFLAMMATORY
list_knis_infl <- DEres_combined_geneset |> 
  filter(gene_biotype == "protein_coding", 
         grepl("HALLMARK_INFLAMMATORY_RESPONSE", pathways_h),
         source == "Knisbacher", 
         abs(log2FoldChange) > 1
  ) |> 
  arrange(pvalue)

#select top 50 genes
list_top50_knis <- DEres_combined_geneset |> 
  filter(gene_symbol %in% c(list_knis_tnf$gene_symbol, list_knis_mtor$gene_symbol, list_knis_infl$gene_symbol))  |> 
  #filter(abs(log2FoldChange) > 1, gene_biotype == "protein_coding", source == "Knisbacher") |> 
  arrange(pvalue) |> 
  slice_head(n=50)

#write.csv(list_top50_knis, "~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/list_top50_knis.csv")

# Convert gene symbols to ENSG IDs
symbol_to_ensg <- mapIds(org.Hs.eg.db,
                         keys = DEres_combined_geneset$gene_symbol,
                         column = "ENSEMBL",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Create mapping for labeling, remove NAs
valid_conversions <- symbol_to_ensg[!is.na(symbol_to_ensg)]
ensg_to_symbol <- setNames(names(valid_conversions), valid_conversions)


# Get normalized expression data for genes
vst_counts <- assay(vsd)

#remove ending
rownames(vst_counts) <-sub("\\..*", "", rownames(vst_counts))

# Filter for genes you want to display
heatmap_data <- vst_counts[rownames(vst_counts) %in% valid_conversions, ]

# Convert ENSEMBL IDs to gene symbols
rownames(heatmap_data) <- ensg_to_symbol[rownames(heatmap_data)]

#create list of mutated patients
pat_list <- colnames(heatmap_data) 
mut_status <- setNames(
  ifelse(pat_list %in% DDX3X_list, "mutated", "unmutated"), 
  pat_list
)

#filter genes
heatmap_data <- heatmap_data[rownames(heatmap_data) %in% list_top50_knis$gene_symbol,]

#scale
x <- data.frame(heatmap_data)
y <- data.frame(t(scale(t(x))))

#create matrix
rownames(y) = rownames(heatmap_data)
colnames(y) = colnames(heatmap_data)
z <- as.matrix(y)    

# Define colors for mutation status
mut_colors <- c("mutated" = "#E78AC3", "unmutated" = "grey60")

color_fun_reversed <- colorRamp2(c(2, 0, -2), c("red", "white", "blue"))

ha_top <- HeatmapAnnotation(
  `DDX3X` = mut_status, 
  col = list(`DDX3X` = mut_colors
  ),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    `DDX3X` = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    ))
  )

#make common genes red
list_top50_cps <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/list_top50_cps.csv") |> 
  dplyr::select(-1)

common_sign <- rbind(list_top50_cps, list_top50_knis) |> 
  filter(padj <0.1)

#genes that are in both datasets
common_list <- intersect(list_top50_cps$gene_symbol, list_top50_knis$gene_symbol)

#genes that are shared AND significant
shared_and_significant <- intersect(common_list, common_sign$gene_symbol)

# Create a color vector
color_vector <- ifelse(rownames(z) %in% shared_and_significant, "red", "black")

mut_status <- mut_status[colnames(z)]

#heatmap
knis <- Heatmap(z, name = "Z-score", show_row_names = TRUE, show_column_names=FALSE, 
             row_names_gp = gpar(fontsize = 8, 
                                 fontface = "italic",
                                 col = color_vector),
             column_names_gp = gpar(fontsize = 8),
             show_row_dend = FALSE,
             show_column_dend = FALSE, 
             column_title = "Knisbacher et al. 2022",
             width = ncol(z)*unit(0.5, "mm"), height = nrow(z)*unit(2.5, "mm"), 
             col = color_fun_reversed, 
             top_annotation = ha_top,
             column_split = mut_status,
             row_title_gp = gpar(fontsize = 8), 
             show_parent_dend_line = FALSE, 
             column_title_gp = gpar(fontsize = 8, fontface = "bold"),
             #column_order = column_order,
             heatmap_legend_param = list(
               title_gp = gpar(fontsize = 8, fontface = "bold"),
               labels_gp = gpar(fontsize = 8)
             ))






