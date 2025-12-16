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
  vsd <- vst(dds, blind=FALSE)
  
  #drop genes with padj = NA
  DEres <- results(dds, tidy=TRUE)
  DEres <- DEres |> 
    drop_na(padj, row, stat)
  
  # Prepare ranked gene list from DESeq2 results
  # Sort by the stat column (log2FoldChange divided by standard error)
  ranked_genes <- DEres %>%
    dplyr::arrange(desc(stat)) %>%
    dplyr::select(row, stat)
  
  #add gene names
  dt_ucll2 <- mutate(dt_ucll, gene_id = rownames(dt_ucll), .before = "gene_symbol") |> 
    drop_na(gene_symbol)
  DEres <- merge(dt_ucll2[,1:2], DEres, by.x="gene_id", by.y="row")
  
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
  
  # Create a named vector without duplicates
  gene_ranks <- setNames(ranked_genes_unique$stat, ranked_genes_unique$gene_symbol)
  
  gene_ranks_CPS1000 <- gene_ranks
  #write.csv(gene_ranks_CPS1000, "~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/gene_ranks_CPS1000.csv")
  
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

  
  #plot top genes
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ComplexHeatmap)
  library(circlize)
  
  #load common DE geneset table
  DEres_combined_geneset <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/DEres_combined_geneset.csv")
  
  #remove genes with high correlation (>0.8) with T cell genes (such as CD3E) & remove hemoglobins
  strong_positive_CD3E <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/R notebooks/strong_positive_CD3D.csv")
  gene_remove <- c(strong_positive_CD3E$Gene, "HBB", "HBA2")
  DEres_combined_geneset <- DEres_combined_geneset |> 
    filter(!gene_symbol %in% gene_remove)
  
  #select gene sets
  #TNF
  list_cps_tnf <- DEres_combined_geneset |> 
    filter(gene_biotype == "protein_coding", 
           grepl("HALLMARK_TNFA_SIGNALING_VIA_NFKB", pathways_h),
           source == "CPS1000", 
           abs(log2FoldChange) > 1
           ) |> 
    arrange(pvalue)
  
  #leading edge
  list_le_cps_tnf <- unlist(fgsea_hallmark$leadingEdge[2])[1:20] #HALLMARK_TNFA_SIGNALING_VIA_NFKB
  
  #mTORC
  list_cps_mtor <- DEres_combined_geneset |> 
    filter(gene_biotype == "protein_coding", 
           grepl("HALLMARK_MTORC1_SIGNALING", pathways_h),
           source == "CPS1000", 
           abs(log2FoldChange) > 1
           ) |> 
    arrange(pvalue)
  
  #select leading edge
  list_le_cps_mtor <- unlist(fgsea_hallmark$leadingEdge[7])[1:20] #	HALLMARK_MYC_TARGETS_V1
  
  #INFLAMMATORY
  list_cps_infl <- DEres_combined_geneset |> 
    filter(gene_biotype == "protein_coding", 
           grepl("HALLMARK_INFLAMMATORY_RESPONSE", pathways_h),
           source == "CPS1000", 
           abs(log2FoldChange) > 1
           ) |> 
    arrange(pvalue)
  
  #select leading edge
  list_le_cps_infl <- unlist(fgsea_hallmark$leadingEdge[8])[1:20] #HALLMARK_INFLAMMATORY_RESPONSE
  
  #select top 100 genes
  list_top50_cps <- DEres_combined_geneset |> 
    filter(gene_symbol %in% c(list_cps_tnf$gene_symbol, list_cps_mtor$gene_symbol, list_cps_infl$gene_symbol))  |> 
    #filter(abs(log2FoldChange) > 2, gene_biotype == "protein_coding", source == "CPS1000") |> 
    arrange(pvalue) |> 
    slice_head(n=50)
  
  #write.csv(list_top50_cps, "~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/list_top50_cps.csv")
  
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
  
  # Filter for genes you want to display
  heatmap_data <- vst_counts[rownames(vst_counts) %in% valid_conversions, ]
  
  # Convert ENSEMBL IDs to gene symbols
  rownames(heatmap_data) <- ensg_to_symbol[rownames(heatmap_data)]
  
  #create list of mutated patients
  pat_list <- colnames(heatmap_data) 
  mut_status <- setNames(
    ifelse(pat_list %in% DDX3X_pat_list, "mutated", "unmutated"), 
    pat_list
  )
  
  #filter genes
  heatmap_data <- heatmap_data[rownames(heatmap_data) %in% list_top50_cps$gene_symbol,]
  
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
  
  #load viability data
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
  
  #remove last two characters from each string in 'drug' column
  merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
  merged_table_long
  
  #change order according to mean viability  
  merged_table_conc <- merged_table_long |> 
    mutate(conc = rep(seq(1:5), times = nrow(merged_table_long)/5))
  
  #focus on main diseases
  merged_filtered <- merged_table_conc |>
    filter(str_detect(diagnosis, "CLL"))
  
  #merge with patMeta
  merged_filtered <- merged_filtered |> 
    group_by(drug, patientID)
  
  df <- merge(merged_filtered, patMeta, by.x="patientID", by.y = "Patient.ID") 
  
  #focus on BCR
  df_bcr <- df |>
    filter(drug == "Ibrutinib") |> 
    dplyr::select(patientID, viability) |> 
    group_by(patientID) |> 
    summarize(mean_viab = mean(viability)) |> 
    filter(patientID %in% colnames(z))
  
  ibr_sens <- setNames(df_bcr$mean_viab, df_bcr$patientID) 
  remaining <- rep(NA, times=length(setdiff(pat_list, df_bcr$patientID)))
  names(remaining) <- setdiff(pat_list, df_bcr$patientID)
  
  ibr_sens <- c(ibr_sens, remaining)
  
  # Create top annotation
  #annotation for Ibr sensitivity
  library(RColorBrewer)
  
  # alignment: create a vector with names = colnames(z)
  aligned_ibr <- rep(NA_real_, ncol(z))
  names(aligned_ibr) <- colnames(z)
  
  aligned_mut <- rep(NA_real_, ncol(z))
  names(aligned_mut) <- colnames(z)
  
  # fill values for samples we do have
  common <- intersect(names(ibr_sens), colnames(z))
  aligned_ibr[common] <- ibr_sens[common]
  
  common_mut <- intersect(names(mut_status), colnames(z))
  aligned_mut[common_mut] <- mut_status[common_mut]
  
  # use the aligned version
  ibr_sens <- aligned_ibr
  mut_status <- aligned_mut
  
  ibr_colors <- colorRamp2(c(min(ibr_sens, na.rm = TRUE), 
                             mean(ibr_sens, na.rm = TRUE),
                             max(ibr_sens, na.rm = TRUE)), 
                           c("royalblue", "white", "orange"))
  
  ha_top <- HeatmapAnnotation(
    `DDX3X` = mut_status, 
    `Ibrutinib (mean viability)` = ibr_sens,
    col = list(`DDX3X` = mut_colors, 
               `Ibrutinib (mean viability)` = ibr_colors
    ),
    na_col = "grey80",
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
      `DDX3X` = list(
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      )
      , `Ibrutinib (mean viability)` = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)
      )
    ))
  
  #column_order <- order(ibr_sens, na.last = TRUE)  
  
  #make common genes red
  list_top50_knis <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/list_top50_knis.csv") |> 
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
  cps <- Heatmap(z, name = "Z-score", show_row_names = TRUE, show_column_names=FALSE, 
               row_names_gp = gpar(fontsize = 8, 
                                   fontface = "italic", 
                                   col = color_vector),
               column_names_gp = gpar(fontsize = 8),
               column_title = "CPS1000",
               width = ncol(z)*unit(0.5, "mm"), height = nrow(z)*unit(2.5, "mm"), 
               col = color_fun_reversed, 
               show_column_dend = FALSE,
               show_row_dend = FALSE, 
               column_split = mut_status,
               top_annotation = ha_top,
               row_title_gp = gpar(fontsize = 8), 
               show_parent_dend_line = FALSE, 
               column_title_gp = gpar(fontsize = 8, fontface = "bold"),
               #column_order = column_order,
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 8, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8)
               ))
  
  
  #select NLRP3 gene
  # Get normalized expression data for genes
  normalized_counts <- counts(dds, normalized = TRUE)
  
  heatmap_data <- normalized_counts[rownames(normalized_counts) %in% valid_conversions, ]
  rownames(heatmap_data) <- ensg_to_symbol[rownames(heatmap_data)]
  
  #plot
  library(ggpubr)
  mycolors <- c("#E78AC3", "grey60")
  names(mycolors) <- c("mutated","unmutated")
  
  nlrp3 <- heatmap_data[rownames(heatmap_data) == "NLRP3"]
  nlrp3 <- as.data.frame(nlrp3) 
  nlrp3$patientID <- colnames(heatmap_data)
  nlrp3 <- nlrp3 |> 
    mutate(DDX3X = ifelse(patientID %in% DDX3X_pat_list, "mutated", "unmutated"))
  
SFig11e <- nlrp3 |> 
    ggplot(aes(x=DDX3X, y=log2(nlrp3+1)))+
    geom_boxplot(aes(color = DDX3X), outliers = FALSE)+
    geom_jitter(aes(color = DDX3X), width=0.05, alpha=0.5)+
    stat_compare_means(aes(group = DDX3X), method = "t.test", label = "p.signif", label.x= 1.5, label.y = 10, size=3)+
    labs(title="NLRP3", y="Normalized expression")+
    scale_color_manual(values=mycolors) +
    theme_classic() + theme(
    axis.text = element_text(color="black", size=8), 
    axis.title = element_text(size=8),
    legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), 
    legend.key.size = unit(0.5, 'cm'),
    axis.line = element_line(size = 0.5),
    plot.title = element_text(hjust = 0.5, size=8, face="bold"),
    legend.position = "none", plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "cm"))
SFig11e
  
  #add sex and select DDX3Y 
  theme_set(theme_bw() + theme(
    strip.text.y = element_text(size = 0), 
    strip.text.x = element_text(size=8),
    axis.text = element_text(size =8, color="black"), 
    axis.title = element_text(size =8, color="black"),
    panel.spacing = unit(0.15, "lines"), 
    legend.title = element_text(size=8, face="bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    strip.background = element_rect(color = "black", linewidth = 1, fill = "grey90"),        
    plot.title = element_text(hjust=0.5, size=8, face="bold")
  ))
  
  DDX3Y <- heatmap_data[rownames(heatmap_data) == "DDX3Y",]
  DDX3Y <- as.data.frame(DDX3Y) 
  DDX3Y$patientID <- rownames(DDX3Y) 
  DDX3Y <- DDX3Y |> 
    merge(patMeta[,c("Patient.ID", "gender")], by.x="patientID", by.y="Patient.ID") |> 
    mutate(DDX3X = ifelse(patientID %in% DDX3X_pat_list, "mutated", "unmutated"))
  
  SFig11f <- DDX3Y |>
    mutate(gender = case_when(gender == "f" ~ "Female",
                              TRUE ~ "Male")) |> 
    ggplot(aes(x=DDX3X, y=log2(DDX3Y+1)))+
    geom_boxplot(aes(color = DDX3X), outliers = FALSE)+
    geom_jitter(aes(color = DDX3X), width=0.05, alpha=0.5)+
    facet_wrap(~gender, ncol=2)+
    stat_compare_means(aes(group = DDX3X), method = "t.test", label = "p.signif", label.x= 1.5, label.y = 16, size=3)+
    labs(title="DDX3Y", y="Normalized expression")+
    scale_color_manual(values=mycolors)+
    theme(legend.position = "none", plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
  SFig11f 
  
  
  
 
