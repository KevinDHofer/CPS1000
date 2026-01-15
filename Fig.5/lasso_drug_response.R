#Libraries
library(Biobase)
library(gridExtra)
library(reshape2)
library(pheatmap)
library(genefilter)
library(Rtsne)
library(smallvis)
library(ggrepel)
library(RColorBrewer)
library(colorspace)
#library(jyluMisc)
library(robustbase)
library(ggplot2)
library(ggbeeswarm)
library(glmnet)
#library(ggtern)
library(gtable)
library(readxl)
library(survival)
library(maxstat)
library(survminer)
library(tidyverse)
library(dplyr)  
library(BloodCancerMultiOmics2017)
library(DESeq2)

set.seed(1118)

usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

screenData <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/patmeta_210324.RData")

#load table with disease entities  
diseases <- usedSamples_all |> 
  arrange(patientID) |> 
  dplyr::select (patientID, diagnosis) 

# merge the tables
merged_table <- cbind(diseases, screenData) |> 
  relocate(diagnosis, .before = X10058.F4_1)

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:318), names_to = "drug", values_to = "viability")
merged_table_long$drug <- gsub('.{2}$', '', merged_table_long$drug)

#1. Prepare drug screening data (individual concentrations), only CLL
viabMat <- filter(merged_table_long, diagnosis %in% c("CLL")) %>% 
  group_by(patientID, drug) %>% summarise(viab = mean(viability)) %>%  
  ungroup() %>%
  spread(key = patientID, value = "viab") %>% data.frame() %>%
  column_to_rownames("drug") %>% as.matrix() |> t()

## Data pre-processing

#2. Transcriptomic data: top 20 PCs of 5000 most variant genes

# load the most up to date datasets
load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/ddsrna_180717.RData")

#only CLL samples
dds$diag <- patMeta[match(dds$PatID, patMeta$Patient.ID),]$diagnosis
dds <- dds[,dds$diag %in% "CLL"]
dds <- estimateSizeFactors(dds)

#filter out non-protein coding genes
dds<-dds[rowData(dds)$biotype %in% "protein_coding",]

#filter out low count genes
minrs <- 100
rs  <- rowSums(counts(dds))
dds<-dds[ rs >= minrs, ]

#variance stabilize the data
#(includes normalizing for library size and dispsersion estimation) 
dds<-vst(dds)
exprMat <- assay(dds)

#filter out low variable genes
ntop<-5000
sds <- genefilter::rowSds(exprMat)
exprMat<-exprMat[order(sds,decreasing = T)[1:ntop],]
exprMat <- t(exprMat)

#Prepare PCA
pcRes <- prcomp(exprMat, center = T, scale. = TRUE)
rnaPCA <- pcRes$x[,1:20]
pcLoad <- pcRes$rotation[,1:20]

#3. Methylation array data: data from BlooodCancerMultiOmics2017, calculate top 20 PCs
meth_array <- load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/encPatientID_160218.RData")

#use methylation data in PACE
data("methData")
methData = t(assay(methData))
methPCA <- prcomp(methData, center = T, scale. = TRUE)$x[,1:20]
rownames(methPCA) <- encPatientID[match(rownames(methPCA), encPatientID$PatientID2),]$PatientID

#4. Genetic data (CNV, SNV): remove features with >20% missing data, then set missing data as 0 (wt).
#genetics
genData <- filter(patMeta, Patient.ID %in% rownames(viabMat)) %>%
  dplyr::select(-gender, -diagnosis, -IGHV.status, -Methylation_Cluster, -treatment) %>%
  mutate_at(vars(-Patient.ID), as.character) %>%
  mutate_at(vars(-Patient.ID), as.integer) %>%
  data.frame() %>% column_to_rownames("Patient.ID")

#remove gene with higher than 20% missing values
genData <- genData[,colSums(is.na(genData))/nrow(genData) <= 0.2]

#fill the missing value with majority
genData <- apply(genData, 2, function(x) {
  xVec <- x
  avgVal <- mean(x,na.rm= TRUE)
  if (avgVal >= 0.5) {
    xVec[is.na(xVec)] <- 1
  } else xVec[is.na(xVec)] <- 0
  xVec
})

#5. IGHV and methylation cluster
#create patBack
patBack <- filter(patMeta, Patient.ID %in% rownames(viabMat)) %>%
  dplyr::select(-project, -date.of.diagnosis, 
                #-treatment, 
                -date.of.first.treatment, -HIPO.ID)

#IGHV
translation <- c(`U` = 0, `M` = 1)
stopifnot(all(patBack$IGHV.status %in% c("U","M", NA)))
IGHVData <- matrix(translation[as.character(patBack$IGHV.status)], 
                   dimnames = list(patBack$Patient.ID, "IGHV"), ncol = 1)

#remove patiente with NA IGHV status
IGHVData<-IGHVData[!is.na(IGHVData), ,drop=F]

# Methylation cluster
translation <- c(`HP` = 2, `IP` = 1, `LP` = 0)
Mcluster <- matrix(translation[as.character(patBack$Methylation_Cluster)],
                   dimnames = list(patBack$Patient.ID, "ConsCluster"), ncol = 1)
Mcluster <- Mcluster[!is.na(Mcluster), ,drop=F]


#6. Pre-treatment: whether the patient was treated before screen sample was taken.
#pretreated <- matrix(patBack$pretreat, dimnames = list(patBack$Patient.ID, "pretreat"), ncol = 1)
#pretreated <- pretreated[!is.na(pretreated),,drop=FALSE]


#Function to clean and combine multi-omics data
generateData <- function(inclSet, onlyCombine = FALSE, censor = NULL, robust = FALSE) {
  
  dataScale <- function(x, censor = NULL, robust = FALSE) {
    #function to scale different variables
    if (length(unique(na.omit(x))) == 2){
      #a binary variable, change to -0.5 and 0.5 for 1 and 2
      x - 0.5
    } else if (length(unique(na.omit(x))) == 3) {
      #catagorical varialbe with 3 levels, methylation_cluster, change to -0.5,0,0.5
      (x - 1)/2
    } else {
      if (robust) {
        #continuous variable, centered by median and divied by 2*mad
        mScore <- (x-median(x,na.rm=TRUE))/(1.4826*mad(x,na.rm=TRUE))
        if (!is.null(censor)) {
          mScore[mScore > censor] <- censor
          mScore[mScore < -censor] <- -censor
        }
        mScore/2
      } else {
        mScore <- (x-mean(x,na.rm=TRUE))/(sd(x,na.rm=TRUE))
        if (!is.null(censor)) {
          mScore[mScore > censor] <- censor
          mScore[mScore < -censor] <- -censor
        }
        mScore/2
      }
    }
  }
  
  
  allResponse <- list() #list for storing response vector (drug viability)
  allExplain <- list() #list for storing explain matrices (multi-omics)
  
  for (drug in colnames(inclSet$drugs)) {
    y <- inclSet$drugs[,drug]
    
    #get overlapped samples for each dataset 
    overSample <- names(y)
    
    for (eachSet in inclSet) {
      overSample <- intersect(overSample,rownames(eachSet))
    }
    
    y <- dataScale(y[overSample], censor = censor, robust = robust)
    allResponse[[drug]] <- y
  }
  
  #generate explanatory variable table for each seahorse measurement
  expTab <- list()
  
  if ("gen" %in% names(inclSet)) {
    geneTab <- inclSet$gen[overSample,]
    #at least 3 mutated sample
    geneTab <- geneTab[, colSums(geneTab) >= 3]
    vecName <- sprintf("genetic(%s)", ncol(geneTab))
    expTab[[vecName]] <- apply(geneTab,2,dataScale)
  }
  
  if ("RNA" %in% names(inclSet)){
    
    #for PCA
    rnaPCA <- inclSet$RNA[overSample, ]
    colnames(rnaPCA) <- paste0("con.expression",colnames(rnaPCA), sep = "")
    expTab[["expression(20)"]] <- apply(rnaPCA,2,dataScale, censor = censor, robust = robust)
    
  }
  
  if ("meth" %in% names(inclSet)){
    methPCA <- inclSet$meth[overSample,]
    colnames(methPCA) <- paste("con.methylation",colnames(methPCA),sep = "")
    expTab[["methylation(20)"]] <- apply(methPCA[,1:20],2,dataScale, censor = censor, robust = robust)
  }
  
  if ("IGHV" %in% names(inclSet)) {
    IGHVtab <- inclSet$IGHV[overSample,,drop=FALSE]
    expTab[["IGHV(1)"]] <- apply(IGHVtab,2,dataScale)
  }
  
  if ("Mcluster" %in% names(inclSet)) {
    methTab <- inclSet$Mcluster[overSample,,drop=FALSE]
    colnames(methTab) <- "methylation cluster"
    expTab[["methCluster(1)"]] <- apply(methTab,2,dataScale)
  }
  
  if ("pretreated" %in% names(inclSet)){
    preTab <- inclSet$pretreated[overSample,,drop=FALSE]
    expTab[["pretreated(1)"]] <- apply(preTab,2,dataScale)
  }
  
  comboTab <- c()
  
  for (eachSet in names(expTab)){
    comboTab <- cbind(comboTab, expTab[[eachSet]])
  }
  vecName <- sprintf("all(%s)", ncol(comboTab))
  expTab[[vecName]] <- comboTab
  
  allExplain <- expTab
  
  if (onlyCombine) {
    #only return combined results, for feature selection
    allExplain <-expTab[length(expTab)]
  }
  
  return(list(allResponse=allResponse, allExplain=allExplain))
  
}

#Function for multi-variant regression
runGlm <- function(X, y, method = "ridge", repeats=20, folds = 3) {
  modelList <- list()
  lambdaList <- c()
  varExplain <- c()
  coefMat <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)
  
  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }
  
  for (i in seq(repeats)) {
    if (ncol(X) > 2) {
      res <- cv.glmnet(X,y, type.measure = "mse", family="gaussian", 
                       nfolds = folds, alpha = alpha, standardize = FALSE)
      lambdaList <- c(lambdaList, res$lambda.min)
      modelList[[i]] <- res
      
      coefModel <- coef(res, s = "lambda.min")[-1] #remove intercept row
      coefMat[,i] <- coefModel
      
      #calculate variance explained
      y.pred <- predict(res, s = "lambda.min", newx = X)
      varExp <- cor(as.vector(y),as.vector(y.pred))^2
      varExplain[i] <- ifelse(is.na(varExp), 0, varExp) 
      
   
    } else {
      fitlm<-lm(y~., data.frame(X))
      varExp <- summary(fitlm)$r.squared
      varExplain <- c(varExplain, varExp)
      
    }
    
  }
  list(modelList = modelList, lambdaList = lambdaList, varExplain = varExplain, coefMat = coefMat)
}

#Prepare clean data for feature selection, only genetic, IGHV included
inclSet<-list(gen=genData, 
              IGHV=IGHVData, 
              #Mcluster = Mcluster, 
              #pretreated=pretreated,
              drugs=viabMat)
cleanData <- generateData(inclSet, censor = 4, onlyCombine = TRUE)

#Perform lasso regression for feature selection (3-fold repeated CV)
lassoResults <- list()
for (eachMeasure in names(cleanData$allResponse)) {
  dataResult <- list()
  for (eachDataset in names(cleanData$allExplain)) {
    y <- cleanData$allResponse[[eachMeasure]]
    X <- cleanData$allExplain[[eachDataset]]
    
    
    glmRes <- runGlm(X, y, method = "lasso", repeats = 20, folds = 3)
    dataResult[[eachDataset]] <- glmRes 
  }
  lassoResults[[eachMeasure]] <- dataResult
  
}


#Function for the heatmap plot
lassoPlot <- function(lassoOut, cleanData, freqCut = 1, coefCut = 0.01, setNumber = "last") {
  plotList <- list()
  if (setNumber == "last") {
    setNumber <- length(lassoOut[[1]])
  } else {
    setNumber <- setNumber
  }
  
  for (seaName in names(lassoOut)) {
    #for the barplot on the left of the heatmap
    barValue <- rowMeans(lassoOut[[seaName]][[setNumber]]$coefMat)
    freqValue <- rowMeans(abs(sign(lassoOut[[seaName]][[setNumber]]$coefMat)))
    barValue <- barValue[abs(barValue) >= coefCut & freqValue >= freqCut] # a certain threshold
    barValue <- barValue[order(barValue)]
    if(length(barValue) == 0) {
      plotList[[seaName]] <- NA
      next
    }
    
    #for the heatmap and scatter plot below the heatmap
    allData <- cleanData$allExplain[[setNumber]]
    seaValue <- cleanData$allResponse[[seaName]]*2 #back to Z-score
    
    tabValue <- allData[, names(barValue),drop=FALSE]
    ord <- order(seaValue)
    seaValue <- seaValue[ord]
    tabValue <- tabValue[ord, ,drop=FALSE]
    sampleIDs <- rownames(tabValue)
    tabValue <- as.tibble(tabValue)
    
    #change scaled binary back to catagorical
    for (eachCol in colnames(tabValue)) {
      if (strsplit(eachCol, split = "[.]")[[1]][1] != "con") {
        tabValue[[eachCol]] <- as.integer(as.factor(tabValue[[eachCol]]))
      }
      else {
        tabValue[[eachCol]] <- tabValue[[eachCol]]*2 #back to Z-score
      }
    }
    
    tabValue$Sample <- sampleIDs
    #Mark different rows for different scaling in heatmap
    matValue <- gather(tabValue, key = "Var",value = "Value", -Sample)
    matValue$Type <- "mut"
    
    #For continuious value
    matValue$Type[grep("con.",matValue$Var)] <- "con"
    
    #for methylation_cluster
    matValue$Type[grep("cluster",matValue$Var)] <- "meth"
    
    #change the scale of the value, let them do not overlap with each other
    matValue[matValue$Type == "mut",]$Value = matValue[matValue$Type == "mut",]$Value + 10
    matValue[matValue$Type == "meth",]$Value = matValue[matValue$Type == "meth",]$Value + 20
    
    
    #color scale for viability
    idx <- matValue$Type == "con"
    
    myCol <- colorRampPalette(c('dark red','white','dark blue'), 
                              space = "Lab")
    if (sum(idx) != 0) {
      matValue[idx,]$Value = round(matValue[idx,]$Value,digits = 2)
      minViab <- min(matValue[idx,]$Value)
      maxViab <- max(matValue[idx,]$Value)
      limViab <- max(c(abs(minViab), abs(maxViab)))
      scaleSeq1 <- round(seq(-limViab, limViab,0.01), digits=2)
      color4viab <- setNames(myCol(length(scaleSeq1+1)), nm=scaleSeq1)
    } else {
      scaleSeq1 <- round(seq(0,1,0.01), digits=2)
      color4viab <- setNames(myCol(length(scaleSeq1+1)), nm=scaleSeq1)
    }
    
    
    
    #change continues measurement to discrete measurement
    matValue$Value <- factor(matValue$Value,levels = sort(unique(matValue$Value)))
    
    #change order of heatmap
    names(barValue) <-  gsub("con.", "", names(barValue))
    matValue$Var <- gsub("con.","",matValue$Var)
    matValue$Var <- factor(matValue$Var, levels = names(barValue))
    matValue$Sample <- factor(matValue$Sample, levels = names(seaValue))
    #plot the heatmap
    p1 <- ggplot(matValue, aes(x=Sample, y=Var)) + geom_tile(aes(fill=Value), color = "gray") + 
      theme_bw() + scale_y_discrete(expand=c(0,0)) + 
      theme(axis.text.y=element_text(hjust=0, size=8, color="black"), axis.text.x=element_blank(), 
            axis.ticks=element_blank(), panel.border=element_rect(colour="gainsboro"),  
            plot.title=element_text(size=8, color="black", face="bold"), 
            panel.background=element_blank(), panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank()) + 
      xlab("patients") + ylab("") + scale_fill_manual(name="Mutated", values=c(color4viab, `11`="gray96", `12`='black', `21`='#DEEBF7', `22`='#9ECAE1',`23` = '#3182BD'),guide=FALSE) + ggtitle(seaName)
    
    #Plot the bar plot on the left of the heatmap
    barDF = data.frame(barValue, nm=factor(names(barValue),levels=names(barValue)))
    
    p2 <- ggplot(data=barDF, aes(x=nm, y=barValue)) + 
      geom_bar(stat="identity", fill="lightblue", position = "identity", width=.8, alpha=0.75) + 
      theme_bw() + geom_hline(yintercept=0, size=0.5) + scale_x_discrete(expand=c(0,0.5)) + 
      scale_y_continuous(expand=c(0,0)) + coord_flip(ylim=c(min(barValue),max(barValue))) + 
      theme(panel.grid.major=element_blank(), panel.background=element_blank(), axis.ticks.y = element_blank(),
            panel.grid.minor = element_blank(), axis.text=element_text(size=8, color="black"), panel.border=element_blank(), axis.line.x = element_line(color="black", size=0.5))+
      xlab("") + ylab("") + geom_vline(xintercept=c(0.5), color="black")
    
    #Plot the scatter plot under the heatmap
    
    # scatterplot below
    scatterDF = data.frame(X=factor(names(seaValue), levels=names(seaValue)), Y=seaValue)
    
    p3 <- ggplot(scatterDF, aes(x=X, y=Y)) + geom_point(shape=21, fill="grey80", colour="black", size=1.2) + theme_bw() + 
      theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=8, color="black"), panel.border=element_rect(colour="dimgrey", size=0.1), panel.background=element_rect(fill="gray96"))
    
    #Scale bar for continuous variable
    
    Vgg = ggplot(data=data.frame(x=1, y=as.numeric(names(color4viab))), aes(x=x, y=y, color=y)) + geom_point() + 
      scale_color_gradientn(name="Z-score", colours =color4viab) + theme(legend.title=element_text(size=12), legend.text=element_text(size=10))
    
    #Assemble all the plots togehter
    
    # construct the gtable
    wdths = c(1.5, 0.2, 1.3*ncol(matValue), 0.7, 1)
    hghts = c(0.1, 0.0010*nrow(matValue), 0.16, 0.5)
    gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
    
    ## make grobs
    gg1 = ggplotGrob(p1)
    gg2 = ggplotGrob(p2)
    gg3 = ggplotGrob(p3)
    gg4 = ggplotGrob(Vgg)
    
    ## fill in the gtable
    gt = gtable_add_grob(gt, gtable_filter(gg2, "panel"), 2, 1) # add barplot
    gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 2, 3) # add heatmap
    gt = gtable_add_grob(gt, gtable_filter(gg1, "title"), 1, 3) #add title to plot
    gt = gtable_add_grob(gt, gtable_filter(gg3, "panel"), 4, 3) # add scatterplot
    gt = gtable_add_grob(gt, gtable_filter(gg2, "axis-b"), 3, 1) # y axis for barplot
    gt = gtable_add_grob(gt, gtable_filter(gg3, "axis-l"), 4, 2) # y axis for scatter plot
    gt = gtable_add_grob(gt, gtable_filter(gg1, "axis-l"), 2, 4) # variable names
    #gt = gtable_add_grob(gt, gtable_filter(gg4, "guide-box"), 2, 5) # scale bar for continous variables
    
    
    #plot
    #grid.draw(gt)
    plotList[[seaName]] <- gt
  }
  return(plotList)
}

heatMaps <- lassoPlot(lassoResults, cleanData, freqCut = 0.8, coefCut = 0.05)

#Show some examples, based on figure 5B in the JCI paper.**
showDrugs <- c("Ibrutinib", "PRT062607", "ONO.4059")

plotList <- heatMaps[showDrugs]
plotList <- plotList[!is.na(plotList)]
Fig5d <- gridExtra::grid.arrange(grobs = plotList,
                        ncol =1)

#Show some examples, based on figure 5B in the JCI paper.**
showDrugs <- c("Fludarabine", "Nutlin.3a")

plotList <- heatMaps[showDrugs]
plotList <- plotList[!is.na(plotList)]
SFig10i <- gridExtra::grid.arrange(grobs = plotList,
                                 ncol =1)