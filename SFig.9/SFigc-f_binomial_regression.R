library(ggplot2)
library(tidyverse)
library(cowplot)
library(glmnet)
library(genefilter)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(mltools)
library(pheatmap)
library(RColorBrewer)

set.seed(1118)

usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/patmeta_210324.RData")

#load table with disease entities  
diseases <- usedSamples_all |> 
  arrange(patientID) |> 
  dplyr::select (patientID, diagnosis) 

# merge the tables
merged_table <- cbind(diseases, viability_table) |> 
  relocate(diagnosis, .before = X10058.F4_1)

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:318), names_to = "drug", values_to = "viability") |> 
  mutate(concIndex = rep(seq(1:5), times=nrow(merged_table)*63))
merged_table_long

#Prepare drug screening data (individual concentrations), only use disease that have at least 10 cases  
viabMat <- filter(merged_table_long, diagnosis %in% c("CLL","AML","B-ALL","T-ALL","MCL","B-PLL","T-PLL")) %>% 
  group_by(patientID, drug, concIndex) %>% summarise(viab = mean(viability)) %>%  
  ungroup() %>% select(-concIndex) %>%
  spread(key = patientID, value = "viab") %>% data.frame() %>%
  column_to_rownames("drug") %>% as.matrix()

#censoring the viability matrix
viabMat[viabMat > 1.2] <- 1.2
viabMat[viabMat < 0] <- 0

#remove features that do not show any variance
sds <- apply(viabMat, 1, sd)
viabMat <- viabMat[sds > genefilter::shorth(sds),] #???

#Remove highly correlated features, highly correlated features may lead to instability of LASSO model**

#function to remove highly correlated values and keep track of the removing
removeCorrelated <- function(x, cutoff = 0.6, method = "pearson", keep = NULL) {
  if (!is.null(keep)) {
    #if specified some feature to keep, then reorder the matrix, to make sure they are not collapsed
    posNow <- grepl(paste(keep, collapse = "|"), colnames(x))
    posNew <- rev(c(colnames(x)[posNow],colnames(x)[!posNow]))
    x <- x[,posNew]
  }
  #input is a feature matrix, features are in columns
  if (method == "binary") {
    #use binary similarity if input is a binary matrix,
    #maybe also usefull is the input is a sparse matrix
    simiMat <- 1 - as.matrix(dist(t(x), method = "binary"))
  } else if (method == "pearson") {
    #otherwise, using pearson correlation
    simiMat <- cor(x)
  } else if (method == "euclidean") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "euclidean"))
  } else if (method == "cosine") {
    # cosine similarity maybe prefered for sparse matrix
    cosineSimi <- function(x){
      x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
    }
    simiMat <- cosineSimi(t(x))
  } else if (method == "canberra") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "canberra"))/nrow(x)
  }
  
  #generate reduced matrix
  simiMat.ori <- simiMat
  simiMat[upper.tri(simiMat)] <- 0
  diag(simiMat) <- 0
  x.re <- x[,!apply(simiMat, 2, function(n) any(abs(n) >= cutoff))]
  
  #a matrix keeping track of the removed features
  mapReduce <- simiMat.ori
  diag(mapReduce) <- 0
  mapList <- lapply(colnames(x.re), function(i) colnames(mapReduce)[mapReduce[i,]>=cutoff])
  names(mapList) <- colnames(x.re)
  
  return(list(reduced = x.re, 
              mapReduce = mapList))
}

keepList <- c("Ibrutinib","Cobimetinib","Venetoclax","OTX","Nutlin")
reRes <- removeCorrelated(t(viabMat), cutoff = 0.7, method = "pearson", keep = keepList)
viabMat.re <- t(reRes$reduced)
mapReduce <- reRes$mapReduce

#Network plot shows the collapsed drugs
edgeTab <- lapply(names(mapReduce), function(n1) {
  if(length(mapReduce[[n1]]) != 0) {
    tibble(from = n1, to = mapReduce[[n1]])
  }
}) %>% bind_rows()

nodeTab <- bind_rows(tibble(name = unique(edgeTab$from), type = "kept"),
                     tibble(name = unique(edgeTab$to), type = "collapsed")) %>%
  mutate(id = seq(nrow(.)))

edgeTab <- mutate(edgeTab, from = nodeTab[match(from, nodeTab$name),]$id,
                  to = nodeTab[match(to, nodeTab$name),]$id)

drugNet <- tbl_graph(nodes = nodeTab, edges = edgeTab, directed = FALSE)

fig1 <- ggraph(drugNet, layout = "igraph", algorithm ="nicely") +
  geom_edge_link(color = "grey90") +
  geom_node_point(aes(color = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size=3) +
  scale_color_manual(values = c(kept = "red", collapsed = "grey50")) +
  labs(color="Correlated drugs")+
  theme_void() +
  theme(legend.title = element_text(size=8, face="bold"), 
             legend.text = element_text(size=8)) 

SFig9c <- fig1 + theme(plot.margin = margin(t=1, l=1, r=1, unit = "cm"))


"Because including highly similar features will lead to the instability 
of multivariate model, I need to first remove some of the highly correlated 
features(drugs that have high similarity of activity profile among samples). 
In the network plot, drugs (plus concentration) that are connected are 
highly similar (Pearson's r > 0.7). I removed the drugs colored by grey 
and kept the drugs colored by red.  If one of the red drugs is identified 
as a signature drug later, the connected grey drugs can also be considered 
as signature drugs."  

#Prepare disease information, separate M-CLL and U-CLL
merged_table_long |> 
  group_by()

diseaseTab <- tibble(patID = colnames(viabMat), 
                     diag = patMeta[match(colnames(viabMat), patMeta$Patient.ID),]$diagnosis, 
                     IGHV = patMeta[match(colnames(viabMat), patMeta$Patient.ID),]$IGHV.status) %>%
  mutate_all(as.character) %>%
  filter(!(diag == "CLL" & is.na(IGHV))) %>% mutate(diag = ifelse(diag == "CLL", paste0(IGHV,"-",diag), diag)) %>%
  mutate(diag = factor(diag, levels = c("U-CLL","M-CLL", "AML","B-ALL","T-ALL","MCL","B-PLL","T-PLL")))

diasaesTab <- diseaseTab |> drop_na(diag)

viabMat.re <- viabMat.re[,diseaseTab$patID]

#Binominal regression (logistic regression), one disease compared to all others

## Fit model
#Function for multi-variant binomial regression
runGlm.bin <- function(X, y, method = "ridge", repeats=20, folds = 3) {
  modelList <- list()
  lambdaList <- c()
  aucList <- c()
  coefMat  <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)
  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }
  
  for (i in seq(repeats)) {
      #balanced sampling
      vecFold <- mltools::folds(y,nfolds = folds, stratified = TRUE)
      res <- cv.glmnet(X,y, type.measure = "auc",
                        foldid = vecFold, alpha = alpha, standardize = FALSE,
                       intercept = TRUE, family = "binomial")
      lambdaList <- c(lambdaList, res$lambda.1se)
      aucList <- c(aucList, res$cvm[res$lambda == res$lambda.1se])
      modelList[[i]] <- res
      
      coefMat[,i] <- coef(res, s = "lambda.1se")[-1]
  }
  list(modelList = modelList, lambdaList = lambdaList, aucList = aucList, coefMat = coefMat)
}

#dropNA
diseaseTab_NA <- diseaseTab[is.na(diseaseTab$diag),]
diseaseTab <- diseaseTab[!is.na(diseaseTab$diag),]

viabMat.re <- viabMat.re[,!(colnames(viabMat.re) %in% diseaseTab_NA$patID)]

diseaseTab$diag <- factor(diseaseTab$diag, levels = unique(diseaseTab$diag))

#Run regression
X <- t(viabMat.re)
binRes <- lapply(unique(diseaseTab$diag), function(disease) {
  y <- ifelse(diseaseTab$diag == disease, 1,0)
  res <- runGlm.bin(X,y, method = "lasso", repeats = 50, folds = 3)
})
names(binRes) <- unique(diseaseTab$diag)

#Cross validation accuracy (auc) for each disease (How accurate each disease type can be predicted by drug response?)
mycolors <- setNames(c("#1F78B4", "#B2DF8A", "#33A02C", "#E31A1C", "#FB9A99", 
           "#FF7F00", "#CAB2D6", "#A6CEE3"), c("MCL", "U-CLL", "M-CLL", 
                                               "AML", "T-ALL", "B-ALL", "B-PLL", "T-PLL"))

pAuc <- lapply(names(binRes), function(n){
  tibble(auc = binRes[[n]][["aucList"]],
         disease = n)
}) %>% bind_rows() %>% group_by(disease) %>%
  summarise(meanAuc = mean(auc), sdAuc = sd(auc)) %>%
  ggplot(aes(x=disease, y = meanAuc, fill = disease)) + geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = meanAuc-sdAuc, ymax=meanAuc+sdAuc ), width = 0.25) +
  scale_fill_manual(values=mycolors) + 
  labs(x="", y="AUC of ROC", fill="Disease") + 
  scale_y_continuous(expand = c(0.01,0))+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8), axis.title.y = element_text(size=8), 
        legend.position = "none", 
        axis.line = element_line(size = 0.5))
pAuc 

SFig9e <- pAuc + theme(plot.margin = margin(t=2, l=3, r=3, b=2, unit = "cm"))

###### Summarise and plot results
sumCoef <- function(coefMat, coefCut = 0, freqCut =1) {
  meanCoef <- rowMeans(abs(coefMat))
  freqCoef <- rowMeans(coefMat != 0)
  subMat <- coefMat[meanCoef > coefCut & freqCoef >= freqCut,,drop=FALSE]
  eachTab <- data.frame(subMat) %>%
     rownames_to_column("drug") %>% gather(key = "rep",value = "coef",-drug) %>%
     mutate(rep = gsub("X","",rep))
  return(eachTab)
} 

coefTab <- lapply(names(binRes), function(d) {
   coefMat <- binRes[[d]]$coefMat
   coefTab <- sumCoef(coefMat=coefMat, freqCut = 1, coefCut =3) %>%
     mutate(disease = d)
}) %>% bind_rows()

#average coefficient for bar plot
meanTab <- coefTab %>% group_by(drug, disease) %>% 
  summarise(meanCoef = mean(coef)) %>% 
  spread(key = disease, value = meanCoef) %>%
  mutate_at(vars(-drug),replace_na,0) %>% 
  gather(key = "disease", value = "coef", -drug)


### Plot coefficient (drug importance) for each disease
plotList <- lapply(unique(coefTab$disease), function(n) {
  plotTab <- filter(coefTab, disease == n) %>%
    group_by(drug) %>% summarise(meanCoef = mean(coef),
                                 sdCoef = sd(coef)) %>%
    arrange(meanCoef) %>% mutate(drug = factor(drug, levels = drug))
  ggplot(plotTab, aes(x=drug, y=meanCoef)) + geom_bar(stat= "identity") +
    geom_errorbar(aes(ymin = meanCoef - sdCoef, ymax = meanCoef +sdCoef)) +
    ggtitle(n) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust =0.5),
                       plot.title = element_text(hjust =0.5, face = "bold")) 
})
gridExtra::grid.arrange(grobs = plotList, ncol=1)

## Plot combined drug response heatmap and coefficient**
#Drug response heatmap
colAnno <- diseaseTab %>% select(patID, diag) %>%
  arrange(diag) %>% data.frame() %>%
  column_to_rownames("patID")
plotMat <- viabMat[unique(meanTab$drug), rownames(colAnno)]

#scale and censoring
plotMat <- t(scale(t(plotMat)))
plotMat[plotMat > 4] <- 4
plotMat[plotMat < -4] <- -4

#hierarchical clustering on rows
hr <- hclust(dist(plotMat), method = "ward.D2")
plotMat <- plotMat[hr$order, ]

#annoation colors
annoColor <- list(Diagnosis = structure(RColorBrewer::brewer.pal(length(unique(colAnno$diag)),"Accent")), 
                                   names = sort(unique(as.character(colAnno$diag))))

breakList <- seq(-4,4,length.out = 100)

ph <- pheatmap(plotMat, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = colAnno, annotation_colors = annoColor,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breakList)),
         breaks = breakList,
         show_colnames = FALSE, show_rownames = TRUE, silent = TRUE)

#Mean coefficient bar plot
plotTab <- meanTab %>% ungroup() %>%
  mutate(drug = factor(drug, levels = rev(rownames(plotMat))))
barPlot <- ggplot(plotTab, aes(x=drug, y = coef, fill = disease)) + 
  geom_bar(stat = "identity", position = "stack") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = mycolors) + theme_bw() +
  labs(x="", y="Coefficient", fill="Disease")+
  theme_bw()+
  theme(legend.position = "right", legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), 
        legend.key.size = unit(0.5, 'cm'), axis.text = element_text(color="black", size=8), 
        axis.title.x = element_text(size=8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', 
                                    fill = NA, 
                                    size = 1)) +
  ylim(-10,10)+
  coord_flip()

pout <- ggdraw() + 
  draw_plot(barPlot)
pout

SFig9d <- pout+ theme(plot.margin = margin(t=1, l=1, r=1, unit = "cm"))

#pdf("combined_binomial.pdf",height = 15, width = 10)
#pout
#dev.off()

## Summarise signature drugs in a network plot

#Function for network plot of associations
plotSigNet <- function(meanTab) {
  absmax <- function(x) { x[which.max(abs(x))]}
  sumTab <- filter(meanTab, coef != 0) %>% 
    separate(drug, into = c("name","conc"),sep = "_") %>% 
    group_by(name, disease) %>%
    summarise(coef = absmax(coef)) %>%
    ungroup()
  
  nodeTab <- bind_rows(select(sumTab, name, coef) %>% mutate(type = "drug"),
                       tibble(name = unique(sumTab$disease), 
                              coef =15,
                              type = "disease")) %>%
    mutate(id = seq(nrow(.))) %>%
    mutate(direction = ifelse(type == "drug", ifelse(coef >0,"resistant","sensitive"),
                              name))  # Use disease name as the direction for disease nodes
  
  # Your color mapping
  mycolors <- setNames(c("#1F78B4", "#B2DF8A", "#33A02C", "#E31A1C", "#FB9A99", 
                         "#FF7F00", "#CAB2D6", "#A6CEE3"), 
                       c("MCL", "U-CLL", "M-CLL", "AML", "T-ALL", "B-ALL", "B-PLL", "T-PLL"))
  
  # Add the drug colors to the color palette
  all_colors <- c(mycolors, resistant = "pink", sensitive = "lightblue")
  
  edgeTab <- sumTab %>% mutate(from = seq(nrow(.)), to = nodeTab[match(disease,nodeTab$name),]$id) %>%
    select(from, to)
  
  tg <- tbl_graph(nodes = nodeTab, edges = edgeTab, directed = FALSE)
  
  ggraph(tg, layout = "igraph", algorithm ="nicely") +  
    geom_edge_link(color = "grey90") +
    geom_node_point(aes(size=abs(coef), shape = type, color = direction)) +
    geom_node_text(aes(label = name), size=3) +
    scale_size_continuous(range = c(1,20), guide = FALSE) +
    scale_shape_manual(values = c(disease= 16, drug = 16), guide = FALSE) +
    scale_color_manual(values = all_colors,
                       breaks = c("resistant", "sensitive"),
                       name = "Direction") +
    theme_void()+
    theme(legend.title = element_text(face="bold", size=8), legend.text = element_text(size=8))
}
set.seed(125)
SFig9f <- plotSigNet(meanTab) + theme(plot.margin = margin(t=1, r=1, l=2, unit = "cm"))
SFig9f 
