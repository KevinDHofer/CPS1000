library(tidyverse)
library(dplyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(survival)
library(ggsurvfit)
library(maxstat)
library(survminer)
library(DescTools)
library(drc)
library(patchwork)
library(scales)

#import raw data
#load tables
usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")
viability_table_all <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viabilityTable_allSamples.csv")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/20250625_timeline.RData")

usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/survival_190516.RData")

#load SMART trial
load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/SMARTrial/SMARTrial_data.RData")
#only 10 patients with ibr monotherapy and response assessment
#no patients with fludarabine treatment


pat_list <- usedSamples_all |> 
  filter(diagnosis == "CLL") |> 
  arrange(patientID) |> 
  pull(patientID)

df_tx <- therapyTab[pat_list] |> 
  bind_rows(.id = "patientID")

#write.csv2(df_tx, "~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/treatments_cps1000.csv")

#nr of individual patients
length(unique(df_tx$patientID))

#treatments
#patients with ibr tx, select all ibr therapy lines
ibr <- df_tx |> 
  filter(regimen %in% c("Ibr")) |> 
  rename(ibr_start=start, ibr_end=end) |> 
  group_by(patientID) |> 
  arrange(patientID, ibr_start) |> 
  #slice_head(n=1) |> 
  mutate(ibr_line = seq(1,n()))

#nr of individual ibr patients
ibr_tx_patID <- unique(ibr$patientID)

#calculate TTNT
ibr_start <- df_tx |> 
  filter(regimen %in% c("Ibr")) |> 
  dplyr::select(patientID, regimen, ibr_start = start, ibr_end=end) |> 
  rename(ibr_regimen = regimen)

#filter for treatments after ibr and calculate ttnt
ibr_ttnt <- df_tx %>%
  filter(patientID %in% ibr_tx_patID) |> 
  merge(ibr_start, by = "patientID") |>
  group_by(patientID, ibr_start) |> 
  filter(start > ibr_start) |> #all treatments starting after ibr
  mutate(ttnt = as.numeric(difftime(as.POSIXct(start), as.POSIXct(ibr_start), units = "days"))) |> 
  arrange(start) |> 
  slice_head(n=1) |> #select next treatment after ibr
  ungroup() |> 
  dplyr::select(patientID, ibr_start, ttnt)

#combine dataframes
ibr <- merge(ibr, ibr_ttnt, by=c("patientID", "ibr_start"), all.x = TRUE)

df_surv <- survT |> 
  dplyr::select(patID, sampleDate, LKA, died, lastUpdate) |> 
  group_by(patID) |> 
  arrange(sampleDate) |> 
  slice_tail(n=1) |> 
  dplyr::select(-sampleDate)

ibr_surv <- merge(ibr, df_surv, by.x="patientID", by.y="patID", all.x = TRUE) |> 
  mutate(next_tx = ifelse(
    is.na(ttnt), FALSE, TRUE),
    ttnt = ttnt/365,
    ttnt = ifelse(
      next_tx==TRUE, ttnt, as.numeric(difftime(as.POSIXct(LKA), as.POSIXct(ibr_start), units = "days"))/365)
  ) |> 
  relocate(next_tx, .after = "ttnt")

#samples
viability_table_all$sampleID <- gsub("-[0-9]$", "", viability_table_all$sampleID)
viability_table_all$sampleID <- gsub("-[2]A$", "", viability_table_all$sampleID)

#select untreated samples before ibrutinib start, select sample closest to ibr start
samples_ibr_tx <- viability_table_all |> 
  filter(patientID %in% ibr_tx_patID) |> 
  arrange(patientID) |> 
  dplyr::select (patientID, sampleID) |> 
  left_join(survT[, colnames(survT) %in% c("sampleID", "sampleDate")],
            by = "sampleID") |> 
  merge(ibr_surv[,c("patientID", "ibr_line", "ibr_start")], by="patientID") |> 
  mutate(diff_d_ibr = as.numeric(difftime(as.POSIXct(sampleDate), as.POSIXct(ibr_start), units = "days")))

ibr_pat_treatments <- df_tx |> 
  filter(patientID %in% ibr_tx_patID) |> 
  merge(ibr_surv[,c("patientID", "ibr_start")], by="patientID")

samples_ibr_tx2 <- samples_ibr_tx %>%
  dplyr::select(-ibr_start) |> 
  left_join(ibr_pat_treatments, by = c("patientID")) %>%
  mutate(during_treatment = sampleDate >= start & sampleDate <= end) |> 
  arrange(patientID, diff_d_ibr)

#filter for samples taken before ibr, without treatment
samples_ibr_tx3 <- samples_ibr_tx2 |> 
  group_by(patientID, sampleID, sampleDate, diff_d_ibr) |>
  mutate(during_treatment_any = ifelse(sum(as.numeric(during_treatment) >0), TRUE, FALSE)) |> 
  filter(diff_d_ibr < 0,during_treatment_any == FALSE) |> #untreated samples before ibr
  dplyr::select(patientID, sampleID, sampleDate, ibr_line, diff_d_ibr, ibr_start, during_treatment) |> 
  distinct() |> 
  arrange(patientID, desc(diff_d_ibr)) |> 
  filter(ibr_line == 1) |>  #filter for first time ibr
  group_by(patientID) |> 
  slice_head(n=1)

ibr_tx_patID_pre <- unique(samples_ibr_tx3$patientID)

#filter for samples during ibr treatment
samples_ibr_tx4 <- samples_ibr_tx2 |> 
  group_by(patientID, sampleID, sampleDate, diff_d_ibr) |>
  mutate(during_treatment_any = ifelse(sum(as.numeric(during_treatment) >0), TRUE, FALSE)) |> 
  filter(diff_d_ibr >0, regimen %in% c("Ibr", "O-Ibr"), during_treatment_any == TRUE) |>
  dplyr::select(patientID, sampleID, sampleDate, diff_d_ibr, ibr_start, during_treatment) |> 
  distinct(sampleDate) |> 
  arrange(patientID, sampleID, desc(diff_d_ibr)) |> 
  group_by(patientID) |> 
  slice_head(n=1)

ibr_tx_patID_post <- unique(samples_ibr_tx4$patientID)

#combine with viability data
diseases <- samples_ibr_tx3[,c("patientID", "sampleID", "sampleDate")]

#combine with viability data: during ibr treatment
diseases <- samples_ibr_tx4[,c("patientID", "sampleID", "sampleDate")]

#merge the tables
merged_table <- merge(diseases, viability_table_all, by=c("patientID", "sampleID"))

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:length(merged_table)), names_to = "drug", values_to = "viability")
merged_table_long

#remove last two characters from each string in 'drug' column
merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
merged_table_long <- merged_table_long |> 
  mutate(conc_level=rep(seq(1,5), times=nrow(merged_table_long)/5)) 


#cap supravital viability results (cut-off <1.2)
merged_table_long$viability[merged_table_long$viability >1.2] <-1.2

#create means and AUC, add IGHV and TP53
load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/patmeta_210324.RData")

#combine treatments with viability:ibr
viabTab.ibr <- merged_table_long |> 
  dplyr::filter(drug == "Ibrutinib") |> 
  mutate(conc = rep(c(2, 0.4, 0.08, 0.016, 0.0032), times=nrow(merged_table_long)/315)) |> 
  group_by(patientID, sampleDate, drug) |> 
  arrange(patientID, drug, conc) |>
  summarise(viab.mean = mean(viability), 
            viab.auc = AUC(conc, viability)*0.5) |> 
  merge(patMeta[,c("Patient.ID", "IGHV.status", "del17p", "TP53")], by.x="patientID", by.y="Patient.ID")

ibr_surv_viab <- merge(ibr_surv, viabTab.ibr, by="patientID") |> 
  filter(!is.na(ttnt), ibr_line == 1) |> 
  mutate(date_diff_d = as.numeric(difftime(as.POSIXct(ibr_start), as.POSIXct(sampleDate), units = "days")))

#plot 
ibr_surv_viab |> 
  ggplot(aes(x=factor(next_tx), y=viab.mean))+
  geom_boxplot(outliers=FALSE)+
  geom_point()+
  stat_compare_means(method="t.test")+
  facet_wrap(~IGHV.status)

#combine treatments with viability:tw37
viabTab.tw37 <- merged_table_long |> 
  dplyr::filter(drug == "TW.37") |> 
  mutate(conc = rep(c(10, 2, 0.4, 0.08, 0.016), times=nrow(merged_table_long)/315)) |> 
  group_by(patientID, sampleDate, drug) |> 
  arrange(patientID, drug, conc) |>
  #summarise(viab.mean = mean(viability), viab.auc = AUC(conc, viability)*0.5) |> 
  merge(patMeta[,c("Patient.ID", "IGHV.status", "del17p", "TP53")], by.x="patientID", by.y="Patient.ID")

#combine treatments with viability
tw37_surv_viab <- merge(ibr_surv, viabTab.tw37, by="patientID") |> 
  filter(!is.na(ttnt)) |> 
  mutate(date_diff_d = as.numeric(difftime(as.POSIXct(ibr_start), as.POSIXct(sampleDate), units = "days"))) |> 
  filter(date_diff_d <0)

#plot 
SFig13g1 <- tw37_surv_viab |> 
  #mutate(tp53_del17p = as.character(ifelse(TP53==1|del17p==1, 1, 0))) |> 
  ggplot(aes(x=factor(conc), y=viability, color=next_tx, ))+
  geom_boxplot(aes(group=interaction(conc, next_tx)), width = 0.5, outliers=FALSE, position = position_dodge(width = 0.75))+
  geom_jitter(aes(group=interaction(conc, next_tx), shape=IGHV.status), size=1, alpha = 0.75, 
              position = position_jitterdodge(jitter.width=0.1, dodge.width = 0.75))+
  labs(title="", y="Viability", x=bquote("Concentration ("*mu*"M)"), color = "Next treatment\nduring follow-up")+
  theme_bw()+
  theme(
    strip.text.y = element_text(size = 0), 
    strip.text.x = element_text(size=8),
    axis.text.y = element_text(size =8, color="black"), 
    axis.text.x = element_text(size =8, color="black", angle=45,hjust=1, vjust=1), 
    axis.title = element_text(size =8, color="black"),
    panel.spacing = unit(0.15, "lines"), 
    legend.title = element_text(size=8, face="bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    strip.background = element_rect(color = "black", linewidth = 1, fill = "grey90"),        
    plot.title = element_text(hjust=0.5, size=8, face="bold"), 
    legend.key.size = unit(0.5, 'cm'))+
  scale_y_continuous(breaks=seq(0.2,1.0, 0.2), limits=c(0.2,1.2)) +
  stat_compare_means(method = "t.test", label="p.signif", size=3, label.y = 1.15)+
  #scale_alpha_manual(name="TP53 and/or del(17p)", values=c(0.5,1), labels=c("unmutated", "mutated"))+
  scale_shape_discrete(name="IGHV", labels=c("unmutated", "mutated"))+
  guides(color = guide_legend(override.aes = list(label = ""))) +
  facet_wrap(~drug)
SFig13g1 

  
#dose-response curves: ibr
viabTab.ibr.allconc <- merged_table_long |> 
  dplyr::filter(drug == "Ibrutinib") |> 
  mutate(conc = rep(c(2, 0.4, 0.08, 0.016, 0.0032), times=nrow(merged_table_long)/315)) |> 
  group_by(patientID, sampleDate, drug) |> 
  arrange(patientID, drug, conc)

viabTab.all.mean <- merged_table_long |> 
  group_by(patientID, sampleDate, drug) |> 
  summarize(mean.viab = mean(viability))
  
  # Fit the dose-response model
  # ibr
  conc_seq <- data.frame(conc = seq(min(viabTab.ibr.allconc$conc), max(viabTab.ibr.allconc$conc), length.out = 500))
  predictions_list <- list()
  
  for (pat in unique(viabTab.ibr.allconc$patientID)) {
  
  # Fit model    
  model_cll <- drm(viability ~ conc, data = viabTab.ibr.allconc[viabTab.ibr.allconc$patientID==pat,], 
                   fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50")))
  
  # Create a new data frame for prediction
  new_data_cll <- data.frame(conc = conc_seq$conc)  # Note: use conc_seq$conc since conc_seq is a data frame
  
  # Create predictions
  new_data_cll$viability <- predict(model_cll, newdata = new_data_cll)
  new_data_cll$patientID <- pat
  
  # Store in list
  predictions_list[[pat]] <- new_data_cll
  }

  # Combine all predictions
  all_predictions_ibr <- bind_rows(predictions_list)
  
  #plot drug response curves
  #ibr
  p0 <- all_predictions_ibr |> 
    merge(ibr_surv_viab, by="patientID", all.x = TRUE) |> 
    filter(!is.na(next_tx)) |> 
    ggplot(aes(x=conc, y=viability, color=next_tx, group=patientID)) +
    geom_line()+
    scale_x_log10()+
    theme_classic()+
    #scale_color_manual(values=mycolors)+
    labs(title = "Ibrutinib", y="Viability", x=bquote("Log Concentration ("*mu*"M)"), color = "Next treatment\nduring follow-up")+
    theme(legend.position = "right", legend.title=element_text(size=8, face="bold"), 
          legend.text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1, color="black", size=8), 
          axis.text.y = element_text(color="black", size=8), axis.line = element_line(size = 0.5), 
          axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
          plot.title = element_text(hjust=0.5, size=8, face="bold"), 
          legend.key.size = unit(0.5, 'cm'))
  p0
  
  #boxplots
  SFig13g2 <- viabTab.ibr.allconc |> 
    merge(ibr_surv_viab, by="patientID", all.x = TRUE) |> 
    filter(!is.na(next_tx)) |> 
    ggplot(aes(x=factor(conc), y=viability, color=next_tx))+
    geom_boxplot(aes(group=interaction(conc, next_tx)), width = 0.5, outliers=FALSE, position = position_dodge(width = 0.75))+
    geom_jitter(aes(group=interaction(conc, next_tx), shape=IGHV.status), size=1, alpha = 0.75, 
                position = position_jitterdodge(jitter.width=0.1, dodge.width = 0.75))+
    labs(title="", y="Viability", x=bquote("Concentration ("*mu*"M)"), color = "Next treatment\nduring follow-up")+
    theme_bw()+
    theme(
      strip.text.y = element_text(size = 0), 
      strip.text.x = element_text(size=8),
      axis.text.y = element_text(size =8, color="black"), 
      axis.text.x = element_text(size =8, color="black", angle=45,hjust=1, vjust=1), 
      axis.title = element_text(size =8, color="black"),
      panel.spacing = unit(0.15, "lines"), 
      legend.title = element_text(size=8, face="bold"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
      strip.background = element_rect(color = "black", linewidth = 1, fill = "grey90"),        
      plot.title = element_text(hjust=0.5, size=8, face="bold"), 
      legend.key.size = unit(0.5, 'cm'))+
    scale_y_continuous(breaks=seq(0.4,1.0, 0.2), limits=c(0.4,1.2)) +
    stat_compare_means(method = "t.test", label="p.signif", size=3, label.y = 1.15)+
    #scale_alpha_manual(name="TP53 and/or del(17p)", values=c(0.5,1), labels=c("unmutated", "mutated"))+
    scale_shape_discrete(name="IGHV", labels=c("unmutated", "mutated"))+
    guides(color = guide_legend(override.aes = list(label = ""))) +
    facet_wrap(~drug.x)
  SFig13g2
    
  SFig13g <- SFig13g2+SFig13g1 + theme(
    plot.margin = margin(t=1, b=1, l=1, unit = "cm")) +
    plot_layout(guides = "collect")
  
  #follow-up 
  ibr_fu <- viabTab.ibr.allconc |> 
    merge(ibr_surv_viab, by="patientID", all.x = TRUE) |> 
    filter(!is.na(next_tx)) |> 
    mutate(fu = as.numeric(difftime(as.POSIXct(LKA), as.POSIXct(ibr_start), units = "days"))/365) |> 
    group_by(patientID) |> 
    slice_head(n=1) |> 
    dplyr::select(patientID, fu)
  summary(ibr_fu)
  
  #drug response parameters: run junyans drug curve script on new samples?
  load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/ic50Tab.RData")
  #check overlap: 6 samples
  check <- viability_table |> 
    filter(sampleID %in% viabTab.ibr.allconc$sampleID)
  
  #Fit sigmoid curve for each drug and each sample
  sumData <- viabTab.ibr.allconc |> 
    dplyr::rename(Drug = drug, viab = viability) |> 
    dplyr::select(patientID, sampleID, Drug, conc, viab)
  
  testF <- function(m0, m1) {
    Fstat <- ((sum(residuals(m0)^2) - sum(residuals(m1)^2))/(m0$df.residual-m1$df.residual))/(sum(residuals(m1)^2)/m1$df.residual)
    1 - pf(Fstat, m0$df.residual-m1$df.residual, m1$df.residual)
  }
  
  calcAUC <- function(m1, minConc =0, maxConc=15, n=100) {
    concSeq <- data.frame(conc = seq(minConc, maxConc, length.out = n))
    valueConc <- data.frame(viab = predict(m1, concSeq),
                            conc = concSeq[,1])
    valueConc <- valueConc[order(valueConc$conc),]
    areaTotal <- 0
    for (i in seq(nrow(valueConc)-1)) {
      areaEach <- (valueConc$viab[i] + valueConc$viab[i+1])*(valueConc$conc[i+1]-valueConc$conc[i])*0.5
      areaTotal <- areaTotal + areaEach
    }
    areaTotal/(valueConc[nrow(valueConc),]$conc-valueConc[1,]$conc)
  }
  
  calcParm <- function(df_input) {
    df <- as.data.frame(df_input)
    df$conc <- as.numeric(df$conc)
    df$viab <- as.numeric(df$viab)
    
    formula_obj <- viab ~ conc
    
    out <- tryCatch({
      m1 <- eval(substitute(
        drm(formula_obj,
            data = df,
            fct = LL2.3u(),
            robust = "tukey",
            lowerl = c(0, 0, NA),
            upperl = c(4, 1, NA)),
        list(formula_obj = formula_obj, df = df)
      ))
      
      fitRes <- m1$coefficients
      m0 <- lm(viab ~ 1, data = df)
      pred_vals <- predict(m1)
      r2 <- cor(pred_vals, df$viab, use = "pairwise.complete.obs")^2
      auc <- calcAUC(m1)
      
      tibble(
        HS = fitRes[[1]],
        Einf = fitRes[[2]],
        IC50 = fitRes[[3]],
        pF = testF(m0, m1),
        R2 = r2,
        AUC = auc
      )
    }, error = function(e) {
      message("Fit failed for a group: ", e$message)
      tibble(
        HS = NA_real_,
        Einf = NA_real_,
        IC50 = NA_real_,
        pF = NA_real_,
        R2 = NA_real_,
        AUC = NA_real_
      )
    })
    
    return(out)
  }
  
  ic50Tab_ibr <- sumData %>%
    group_by(patientID, Drug) %>% nest() %>%
    mutate(mm = map(data, calcParm)) %>%
    unnest(mm) %>%
    dplyr::select(-data) %>% ungroup()
  
  #poor fitting?
  
  #inspect
  annotation_data <- ic50Tab_ibr |> 
    distinct(patientID, HS)
  
  all_predictions_ibr |> 
    merge(ic50Tab_ibr, by="patientID", all.x = TRUE) |> 
    merge(ibr_surv_viab, by="patientID") |> 
    ggplot(aes(x=conc, y=viability, group=patientID, color=next_tx)) +
    geom_line()+
    geom_text(data=annotation_data, aes(label=paste0("HS=", round(HS, 2)), x=0.1, y=0.6), 
              size=3, inherit.aes = FALSE)+
    geom_point(data=viabTab.ibr.allconc, aes(x=conc, y=viability, group=patientID), inherit.aes = FALSE)+
    scale_x_log10()+
    theme_bw()+
    facet_wrap(~patientID, ncol=6)+
    #scale_color_manual(values=mycolors)+
    labs(title = "Ibrutinib drug response", y="Viability", x=bquote("Log Concentration ("*mu*"M)"), color = "Next treatment")+
    theme(legend.position = "right", legend.title=element_text(size=8, face="bold"), 
          legend.text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1, color="black", size=8), 
          axis.text.y = element_text(color="black", size=8),  
          axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
          plot.title = element_text(hjust=0.5, size=8, face="bold"), 
          legend.key.size = unit(0.5, 'cm'))
  
  #poor fit
  
  #combine
  ibr_cox_surv_viab <- merge(ibr_surv_viab, ic50Tab_ibr, by="patientID")
  
  #cox regression
  #univariate
  com_uni <- function(response, time, endpoint, scale =FALSE) {  
    
    if (scale) {
      #calculate z-score
      response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
    }
    
    tryCatch({
      surv <- coxph(Surv(time, endpoint) ~ response) 
      resTab <- tibble(p = summary(surv)[[7]][,5], 
                       HR = summary(surv)[[7]][,2], 
                       lower = summary(surv)[[8]][,3], 
                       higher = summary(surv)[[8]][,4],
                       concordance = summary(surv)$concordance[1],
                       concordance_se = summary(surv)$concordance[2])
    }, error = function(err) {
      resTab <- tibble(p = NA, 
                       HR = NA, 
                       lower = NA, 
                       higher = NA,
                       concordance = NA,
                       concordance_se = NA)
    })
    
  }
  
  
  #univariate analysis for ibr response
  parameters <- colnames(ic50Tab_ibr[,c(3:8)])
  parameters <- append(parameters, c("viab.mean", "viab.auc"))
  drugs <- c("Ibrutinib")
  
  resTTNT_uni <- data.frame()
  
  for (par in parameters) {
    for (drug in drugs) {
      res <- filter(ibr_cox_surv_viab, !is.na(HS)) %>%
        do(com_uni(.[[par]], .$ttnt, .$next_tx, TRUE))
      
      resTTNT_uni_row <- data.frame(
        parameter = par,
        cox = res,
        Drug = drug
      )
      resTTNT_uni <- rbind(resTTNT_uni, resTTNT_uni_row)
      
    }
  }
resTTNT_uni <- resTTNT_uni |> 
  arrange(Drug, cox.p)


#univariate analysis for all drugs
all_drugs_surv <- merge(
  viabTab.all.mean, ibr_surv_viab[,-which(names(ibr_surv_viab) %in% c("drug", "sampleDate"))], by="patientID") |> 
  pivot_wider(names_from = "drug", values_from = "mean.viab")

resTTNT_uni_all <- data.frame()

drugs <- unique(viabTab.all.mean$drug)

for (drug in drugs) {
    res <- filter(all_drugs_surv, !is.na(ttnt)) %>%
      do(com_uni(.[[drug]], .$ttnt, .$next_tx, TRUE))
    
    resTTNT_uni_row <- data.frame(
      Drug = drug,
      cox = res
    )
    resTTNT_uni_all <- rbind(resTTNT_uni_all, resTTNT_uni_row)
    
  }
resTTNT_uni_all <- resTTNT_uni_all |> 
  arrange(Drug, cox.p) |> 
  arrange(cox.p) %>% 
  mutate(p.adj = p.adjust(cox.p, method = "BH"))

#univariate analysis for ibr response
parameters <- c("del17p", "TP53", "IGHV.status")
drugs <- c("Ibrutinib")

resTTNT_uni <- data.frame()

for (par in parameters) {
  for (drug in drugs) {
    res <- filter(ibr_cox_surv_viab, !is.na(del17p)) %>%
      do(com_uni(.[[par]], .$ttnt, .$next_tx, FALSE))
    
    resTTNT_uni_row <- data.frame(
      parameter = par,
      cox = res,
      Drug = drug
    )
    resTTNT_uni <- rbind(resTTNT_uni, resTTNT_uni_row)
    
  }
}
resTTNT_uni <- resTTNT_uni |> 
  arrange(Drug, cox.p)

#plot
ranking_conc <- resTTNT_uni_all |> 
  arrange(desc(cox.concordance)) |> 
  pull(Drug)
resTTNT_uni_all$Drug <- factor(resTTNT_uni_all$Drug, levels = ranking_conc)

SFig13f <- resTTNT_uni_all |> 
  filter(cox.concordance > 0.75 | Drug %in% c("Ibrutinib", "ONO.4059")) |> 
  ggplot(aes(x=Drug, y=cox.concordance))+
  geom_hline(yintercept=0.5, color="grey", linetype="dashed")+
  geom_hline(yintercept=0.75, color="grey", linetype="dashed")+
  geom_point()+
  geom_errorbar(aes(ymin=cox.concordance-cox.concordance_se, ymax=cox.concordance+cox.concordance_se), width = 0.2)+
  theme_classic()+
  labs(title="Univariate Cox TTNT:\nDrug response during ibrutinib treatment", y="C-index", x="Drug")+
  theme(legend.position = "right", legend.title=element_text(size=8, face="bold"), 
        legend.text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1, color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), 
        axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
        plot.title = element_text(hjust=0.5, size=8, face="bold"), 
        legend.key.size = unit(0.5, 'cm'))+
  scale_y_continuous(limits = c(0,1.1), breaks=c(0.25,0.5,0.75,1))

SFig13f <- SFig13f + theme(plot.margin = margin(t=1, l=1, r=1, b=3, unit = "cm"))

#maxstat: ibr and R2
tw37_surv_viab <- tw37_surv_viab |> 
  group_by(patientID) |> 
  mutate(viab.mean = mean(viability)) |> 
  filter(conc_level == 1)

res.cut <- surv_cutpoint(tw37_surv_viab, time = "ttnt", event = "next_tx",
                         variables = c("viab.mean"))

summary(res.cut)
ttnt_value <- summary(res.cut)$cutpoint
plot(res.cut, "viab.mean")
res.cat <- surv_categorize(res.cut)
head(res.cat)

fit <- survfit(Surv(ttnt, next_tx) ~ viab.mean, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval = TRUE, title="TW37 response")+
  labs(y="Freedom from\nsubsequent treatment", x="Time since starting ibrutinib (years)")



#same for flu: nr of patients with chemo
#treatments
#patients with flu tx, select last flu monotherapy
flu <- df_tx |> 
  filter(regimen %in% c("F", "FC", "FCR", "BR")) |> 
  dplyr::select(-end) |> 
  group_by(patientID) |> 
  arrange(patientID, start) |> 
  slice_tail(n=1) |> 
  rename(flu_start=start)

#nr of individual flu patients
flu_tx_patID <- unique(flu$patientID)

#calculate TTNT
flu_start <- df_tx |> 
  filter(regimen %in% c("F", "FC", "FCR", "BR")) |> 
  dplyr::select(patientID, flu_start = start) |> 
  group_by(patientID) |> 
  arrange(patientID, flu_start) |> 
  slice_tail(n=1)

#filter treatments after flu and calculate ttnt
flu_ttnt <- df_tx %>%
  filter(patientID %in% flu_tx_patID) |> 
  merge(flu_start, by = "patientID") |> 
  filter(start > flu_start) |> 
  mutate(ttnt = as.numeric(difftime(as.POSIXct(start), as.POSIXct(flu_start), units = "days"))) |> 
  group_by(patientID) |> 
  arrange(start) |> 
  slice_tail(n=1) |> 
  dplyr::select(patientID, ttnt)

#combine dataframes
flu <- merge(flu, flu_ttnt, by="patientID", all.x = TRUE)

flu_surv <- merge(flu, df_surv, by.x="patientID", by.y="patID", all.x = TRUE) |> 
  mutate(next_tx = ifelse(
    is.na(ttnt), FALSE, TRUE),
    ttnt = ttnt/365,
    ttnt = ifelse(
      next_tx==TRUE, ttnt, as.numeric(difftime(as.POSIXct(LKA), as.POSIXct(flu_start), units = "days"))/365)
  ) |> 
  relocate(next_tx, .after = "ttnt")

#samples
#select untreated samples before flu start, select sample closest to flu start (???)
samples_flu_tx <- viability_table_all |> 
  filter(patientID %in% flu_tx_patID) |> 
  arrange(patientID) |> 
  dplyr::select (patientID, sampleID) |> 
  left_join(survT[, colnames(survT) %in% c("sampleID", "sampleDate")],
            by = "sampleID") |> 
  merge(flu_surv[,c("patientID", "flu_start")], by="patientID") |> 
  mutate(diff_d_flu = as.numeric(difftime(as.POSIXct(flu_start), as.POSIXct(sampleDate), units = "days"))) |> 
  filter(!is.na(sampleDate))

flu_pat_treatments <- df_tx |> 
  filter(patientID %in% flu_tx_patID) |> 
  merge(flu_surv[,c("patientID", "flu_start")], by="patientID")

samples_flu_tx2 <- samples_flu_tx %>%
  dplyr::select(-flu_start) |> 
  left_join(flu_pat_treatments, by = c("patientID")) %>%
  mutate(during_treatment = sampleDate >= start & sampleDate <= end) |> 
  arrange(patientID, diff_d_flu)

#filter for samples without treatment, include samples taken after treatment
samples_flu_tx3 <- samples_flu_tx2 |> 
  group_by(patientID, sampleID, sampleDate, diff_d_flu) |>
  mutate(during_treatment_any = ifelse(sum(as.numeric(during_treatment) >0), TRUE, FALSE)) |> 
  filter(
    diff_d_flu >= 0, 
    during_treatment_any == FALSE) |> 
  dplyr::select(patientID, sampleID, sampleDate, diff_d_flu, flu_start, during_treatment) |> 
  distinct() |> 
  arrange(patientID, diff_d_flu) |> 
  group_by(patientID) |> 
  slice_head(n=1)

#combine with viability data
diseases <- samples_flu_tx3[,c("patientID", "sampleID", "sampleDate")]

# merge the tables
merged_table <- merge(diseases, viability_table_all, by=c("patientID", "sampleID"))

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:length(merged_table)), names_to = "drug", values_to = "viability")
merged_table_long

#remove last two characters from each string in 'drug' column
merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
merged_table_long <- merged_table_long |> 
  mutate(conc_level=rep(seq(1,5), times=nrow(merged_table_long)/5)) 

#cap supravital viability results (cut-off <1.2)
merged_table_long$viability[merged_table_long$viability >1.2] <-1.2

#create means and AUC, add IGHV and TP53
viabTab.flu <- merged_table_long |> 
  dplyr::filter(drug == "Fludarabine") |> 
  mutate(conc = rep(c(2, 0.4, 0.08, 0.016, 0.0032), times=nrow(merged_table_long)/315)) |> 
  group_by(patientID, sampleDate, drug) |> 
  arrange(patientID, drug, conc) |>
  summarise(viab.mean = mean(viability), 
            viab.auc = AUC(conc, viability)*0.5) |> 
  merge(patMeta[,c("Patient.ID", "IGHV.status", "del17p", "TP53")], by.x="patientID", by.y="Patient.ID")

#combine treatments with viability
flu_surv_viab <- merge(flu_surv, viabTab.flu, by="patientID") |> 
  filter(!is.na(ttnt)) |> 
  mutate(date_diff_d = as.numeric(difftime(as.POSIXct(flu_start), as.POSIXct(sampleDate), units = "days")))

#plot 
flu_surv_viab |> 
  ggplot(aes(x=factor(next_tx), y=viab.mean))+
  geom_boxplot(outliers=FALSE)+
  geom_point()+
  stat_compare_means(method="t.test")+
  facet_wrap(~IGHV.status)
  
#plot dose-response curves
viabTab.flu.allconc <- merged_table_long |> 
  dplyr::filter(drug == "Fludarabine") |> 
  mutate(conc = rep(c(2, 0.4, 0.08, 0.016, 0.0032), times=nrow(merged_table_long)/315)) |> 
  group_by(patientID, sampleDate, drug) |> 
  arrange(patientID, drug, conc)
  
  # Fit the dose-response model
#summary


  # flu
  conc_seq <- data.frame(conc = seq(min(viabTab.flu.allconc$conc), max(viabTab.flu.allconc$conc), length.out = 500))
  predictions_list <- list()
  
  for (pat in unique(viabTab.flu.allconc$patientID)) {
  
  #fit model    
  model_cll <- drm(viability ~ conc, data = viabTab.flu.allconc[viabTab.flu.allconc$patientID==pat,], fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50")))
  
  # Create a new data frame for prediction
  new_data_cll <- cbind(data.frame(pat = predict(model_cll, newdata = new_data_cll)))
  
  # Create predictions
  new_data_cll <- data.frame(conc = conc_seq)
  new_data_cll$viability <- predict(model_cll, newdata = new_data_cll)
  new_data_cll$patientID <- pat
  
  # Store in list
  predictions_list[[pat]] <- new_data_cll
  }
  
  # Combine all predictions
  all_predictions_flu <- bind_rows(predictions_list)

  
  #plot drug response curves
  #flu
  p01 <- all_predictions_flu |> 
    merge(flu_surv_viab, by="patientID", all.x = TRUE) |> 
    filter(!is.na(next_tx)) |> 
    ggplot(aes(x=conc, y=viability, color=next_tx, group=patientID)) +
    geom_line()+
    theme_classic()+
    #scale_color_manual(values=mycolors)+
    labs(title = "Fludarabine", y="Viability", x=bquote("Concentration ("*mu*"M)"), color = "Next treatment\nduring follow-up")+
    theme(legend.position = "right", legend.title=element_text(size=8, face="bold"), 
          legend.text = element_text(size=8), 
          axis.text = element_text(color="black", size=8), axis.line = element_line(size = 0.5), 
          axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
          plot.title = element_text(hjust=0.5, size=8, face="bold"), 
          legend.key.size = unit(0.5, 'cm'))+
    scale_y_continuous(breaks=seq(0.0,1.0, 0.2), limits=c(0.0,1.1))+
    scale_x_log10(
      breaks = 10^(-2:0),
      labels = trans_format("log10", math_format(10^.x)))  
  SFig13d <- p01 + theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))
 
  #follow-up 
  flu_fu <- all_predictions_flu |> 
    merge(flu_surv_viab, by="patientID", all.x = TRUE) |> 
    filter(!is.na(next_tx)) |> 
    mutate(fu = as.numeric(difftime(as.POSIXct(LKA), as.POSIXct(flu_start), units = "days"))/365) |> 
    filter(fu > 0) |> 
    group_by(patientID) |> 
    slice_head(n=1) |> 
    dplyr::select(patientID, fu)
  summary(flu_fu)
  
  #boxplots
  p11 <- viabTab.flu.allconc |> 
    merge(flu_surv_viab, by="patientID", all.x = TRUE) |> 
    filter(!is.na(next_tx)) |> 
    ggplot(aes(x=factor(conc), y=viability, color=next_tx))+
    geom_boxplot(aes(group=interaction(conc, next_tx)), outliers = FALSE)+
    geom_jitter(size=1, alpha = 0.25, position = position_jitterdodge(jitter.width=0.1))+
    scale_y_continuous(breaks=seq(0.4,1.0, 0.2), limits=c(0,1.05))+
    theme_classic()+
    labs(title="Fludarabine response", y="Viability", x=bquote("Concentration ("*mu*"M)"), color = "Next treatment\nduring follow-up")+
    theme(legend.position = "right", legend.title=element_text(size=8, face="bold"), 
          legend.text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1, color="black", size=8), 
          axis.text.y = element_text(color="black", size=8), axis.line = element_line(size = 0.5), 
          axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
          plot.title = element_text(hjust=0.5, size=8, face="bold"), 
          legend.key.size = unit(0.5, 'cm'))+
    facet_wrap(del17p~IGHV.status)
  p11
  
  #drug response parameters: run junyans drug curve script on new samples
  #Fit sigmoid curve for each drug and each sample
  
  #Fit sigmoid curve for each drug and each sample
  sumData <- viabTab.flu.allconc |> 
    dplyr::rename(Drug = drug, viab = viability) |> 
    dplyr::select(patientID, sampleID, Drug, conc, viab)
  
  testF <- function(m0, m1) {
    Fstat <- ((sum(residuals(m0)^2) - sum(residuals(m1)^2))/(m0$df.residual-m1$df.residual))/(sum(residuals(m1)^2)/m1$df.residual)
    1 - pf(Fstat, m0$df.residual-m1$df.residual, m1$df.residual)
  }
  
  calcAUC <- function(m1, minConc =0, maxConc=15, n=100) {
    concSeq <- data.frame(conc = seq(minConc, maxConc, length.out = n))
    valueConc <- data.frame(viab = predict(m1, concSeq),
                            conc = concSeq[,1])
    valueConc <- valueConc[order(valueConc$conc),]
    areaTotal <- 0
    for (i in seq(nrow(valueConc)-1)) {
      areaEach <- (valueConc$viab[i] + valueConc$viab[i+1])*(valueConc$conc[i+1]-valueConc$conc[i])*0.5
      areaTotal <- areaTotal + areaEach
    }
    areaTotal/(valueConc[nrow(valueConc),]$conc-valueConc[1,]$conc)
  }
  
  calcParm <- function(df_input) {
    df <- as.data.frame(df_input)
    df$conc <- as.numeric(df$conc)
    df$viab <- as.numeric(df$viab)
    
    formula_obj <- viab ~ conc
    
    out <- tryCatch({
      m1 <- eval(substitute(
        drm(formula_obj,
            data = df,
            fct = LL2.3u(),
            robust = "tukey",
            lowerl = c(0, 0, NA),
            upperl = c(4, 1, NA)),
        list(formula_obj = formula_obj, df = df)
      ))
      
      fitRes <- m1$coefficients
      m0 <- lm(viab ~ 1, data = df)
      pred_vals <- predict(m1)
      r2 <- cor(pred_vals, df$viab, use = "pairwise.complete.obs")^2
      auc <- calcAUC(m1)
      
      tibble(
        HS = fitRes[[1]],
        Einf = fitRes[[2]],
        IC50 = fitRes[[3]],
        pF = testF(m0, m1),
        R2 = r2,
        AUC = auc
      )
    }, error = function(e) {
      message("Fit failed for a group: ", e$message)
      tibble(
        HS = NA_real_,
        Einf = NA_real_,
        IC50 = NA_real_,
        pF = NA_real_,
        R2 = NA_real_,
        AUC = NA_real_
      )
    })
    
    return(out)
  }
  
  ic50Tab_flu <- sumData %>%
    group_by(patientID, Drug) %>% nest() %>%
    mutate(mm = map(data, calcParm)) %>%
    unnest(mm) %>%
    dplyr::select(-data) %>% ungroup()
  
  #fitting ok??
  
  #combine
  flu_surv_viab <- merge(flu_surv_viab, ic50Tab_flu, by="patientID")
  
  #cox regression
  #univariate
  
  #univariate analysis for flu response
  parameters <- colnames(ic50Tab_flu[,c(3:8)])
  parameters <- append(parameters, c("viab.mean", "viab.auc"))
  drugs <- c("Fludarabine")
  
  resTTNT_uni_flu <- data.frame()
  
  for (par in parameters) {
    for (drug in drugs) {
      res <- filter(flu_surv_viab, !is.na(HS)) %>%
        do(com_uni(.[[par]], .$ttnt, .$next_tx, TRUE))
      
      resTTNT_uni_row <- data.frame(
        parameter = par,
        cox = res,
        Drug = drug
      )
      resTTNT_uni_flu <- rbind(resTTNT_uni_flu, resTTNT_uni_row)
      
    }
  }
  resTTNT_uni_flu <- resTTNT_uni_flu |> 
    arrange(Drug, cox.p)

#plot
p2 <- resTTNT_uni_flu |> 
  filter(parameter %in% c("Einf", "HS", "IC50", "R2", "AUC", "viab.mean")) |> 
  ggplot(aes(x=parameter, y=cox.concordance))+
  geom_hline(yintercept=0.5, color="grey", linetype="dashed")+
  geom_hline(yintercept=0.75, color="grey", linetype="dashed")+
  geom_point()+
  geom_errorbar(aes(ymin=cox.concordance-cox.concordance_se, ymax=cox.concordance+cox.concordance_se), width = 0.2)+
  theme_classic()+
  labs(title="Univariate Cox TTNT:\nFludarabine drug-response metrics", y="C-index", x="Parameter")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), 
        axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
        plot.title = element_text(hjust=0.5, size=8, face="bold"), 
        legend.key.size = unit(0.5, 'cm'))+
  scale_y_continuous(limits=c(0,1.1), breaks=c(0,0.25,0.5,0.75,1))+
  scale_x_discrete(labels = c(
    "AUC",
    expression(E[inf]),   # subscript
    "Hill coefficient", 
    expression(IC[50]),  # subscript
    expression(R^2),    # superscript
    "Mean viability"))

SFig13e <- p2 + theme(plot.margin = margin(t=1, l=1, r=1, b=3, unit = "cm"))

#univariate analysis for genetic markers on flu response
parameters <- c("del17p", "TP53", "IGHV.status")
drugs <- c("Fludarabine")

for (par in parameters) {
  for (drug in drugs) {
    res <- filter(flu_surv_viab, !is.na(del17p)) %>%
      do(com_uni(.[[par]], .$ttnt, .$next_tx, FALSE))
    
    resTTNT_uni_row <- data.frame(
      parameter = par,
      cox = res,
      Drug = drug
    )
    resTTNT_uni_flu <- rbind(resTTNT_uni_flu, resTTNT_uni_row)
    
  }
}
resTTNT_uni_flu <- resTTNT_uni_flu |> 
  arrange(Drug, cox.p)  

#univariate analysis for all drugs
all_drugs_surv_flu <- merge(
  viabTab.all.mean, flu_surv_viab[,-which(names(flu_surv_viab) %in% c("drug", "sampleDate"))], by="patientID") |> 
  pivot_wider(names_from = "drug", values_from = "mean.viab")

resTTNT_uni_all_flu <- data.frame()

drugs <- unique(viabTab.all.mean$drug)

for (drug in drugs) {
  res <- filter(all_drugs_surv_flu, !is.na(ttnt)) %>%
    do(com_uni(.[[drug]], .$ttnt, .$next_tx, TRUE))
  
  resTTNT_uni_row <- data.frame(
    Drug = drug,
    cox = res
  )
  resTTNT_uni_all_flu <- rbind(resTTNT_uni_all_flu, resTTNT_uni_row)
  
}
resTTNT_uni_all_flu <- resTTNT_uni_all_flu |> 
  arrange(Drug, cox.p) |> 
  arrange(cox.p) %>% 
  mutate(p.adj = p.adjust(cox.p, method = "BH"))


#assess cytarabine and doxorubicine in AML
aml_pat <- usedSamples_all |> 
  filter(diagnosis == "AML") |> 
  pull(patientID)

#for AraC and Doxo: nr of patients with chemo, select samples before chemo
patAML <- survT |> 
  filter(patID %in% aml_pat) |> 
  group_by(patID) |> 
  arrange(sampleDate) |> 
  slice_head(n=1) |> 
  mutate(sampleBeforeTx = ifelse(sampleDate <= firstTreatDate, TRUE, FALSE), .before ="firstTreatDate") |> 
  filter(sampleBeforeTx == TRUE) |> 
  rename(patientID = patID)

#combine with viability data
aml_merged <- merge(patAML, viability_table, by=c("patientID", "sampleID"))

#transfrom in long data frame
merged_table_long <- aml_merged |>
  pivot_longer(cols = c(17:length(aml_merged)), names_to = "drug", values_to = "viability")
merged_table_long

#remove last two characters from each string in 'drug' column
merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
merged_table_long <- merged_table_long |> 
  mutate(conc_level=rep(seq(1,5), times=nrow(merged_table_long)/5)) 

#cap supravital viability results (cut-off <1.2)
merged_table_long$viability[merged_table_long$viability >1.2] <-1.2

#create means and AUC
viabTab.aml.dox <- merged_table_long |> 
  dplyr::filter(drug %in% c("Doxorubicine")) |> 
  mutate(conc = rep(c(0.4, 0.08, 0.016, 0.0032, 0.00064), times=nrow(merged_table_long)/315))

viabTab.aml.cyt <- merged_table_long |> 
  dplyr::filter(drug %in% c("Cytarabine")) |> 
  mutate(conc = rep(c(10, 2, 0.4, 0.08, 0.016), times=nrow(merged_table_long)/315))

viabTab.aml <- rbind(viabTab.aml.cyt, viabTab.aml.dox) |> 
  group_by(patientID, drug) |> 
  mutate(viab.mean = mean(viability), 
         viab.auc = AUC(conc, viability)*0.5) |> 
  arrange(patientID, drug, conc)

#plot 
viabTab.aml |> 
  filter(!is.na(died)) |> 
  ggplot(aes(x=factor(died), y=viab.mean))+
  geom_boxplot(outliers=FALSE)+
  geom_point()+
  stat_compare_means(method="t.test")+
  facet_wrap(~drug)


# Fit the dose-response model
#summary

#AraC
conc_seq <- data.frame(conc = seq(min(viabTab.aml[viabTab.aml$drug=="Cytarabine",]$conc), max(viabTab.aml[viabTab.aml$drug=="Cytarabine",]$conc), length.out = 500))
predictions_list <- list()

for (pat in unique(viabTab.aml$patientID)) {
  
  #fit model    
  model_aml_cyt <- drm(viability ~ conc, data = viabTab.aml[viabTab.aml$patientID==pat & viabTab.aml$drug=="Cytarabine",], 
                       fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50")))
  
  # Create a new data frame for prediction
  new_data_aml_cyt <- data.frame(conc = conc_seq$conc)
  
  # Create predictions
  new_data_aml_cyt$viability <- predict(model_aml_cyt, newdata = new_data_aml_cyt)
  new_data_aml_cyt$patientID <- pat
  
  # Store in list
  predictions_list[[pat]] <- new_data_aml_cyt
}

# Combine all predictions
all_predictions_aml_cyt <- bind_rows(predictions_list)

#plot drug response curves
#AraC
p02 <- all_predictions_aml_cyt |> 
  merge(viabTab.aml[,c("patientID", "died", "OS")], by="patientID", all.x = TRUE) |> 
  filter(!is.na(died)) |> 
  ggplot(aes(x=conc, y=viability, color=patientID, linetype=died, group=patientID)) +
  geom_line()+
  scale_x_log10()+
  theme_classic()+
  #scale_color_manual(values=mycolors)+
  labs(title = "Cytarabine in AML", y="Viability", x=bquote("Log Concentration ("*mu*"M)"), color = "PatientID")+
  theme(legend.position = "right", legend.title=element_text(size=8, face="bold"), 
        legend.text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1, color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), axis.line = element_line(size = 0.5), 
        axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
        plot.title = element_text(hjust=0.5, size=8, face="bold"), 
        legend.key.size = unit(0.5, 'cm'))
p02

#same for Doxo
conc_seq <- data.frame(conc = seq(min(viabTab.aml[viabTab.aml$drug=="Doxorubicine",]$conc), max(viabTab.aml[viabTab.aml$drug=="Doxorubicine",]$conc), length.out = 500))
predictions_list <- list()

for (pat in unique(viabTab.aml$patientID)) {
  
  #fit model    
  model_aml_dox <- drm(viability ~ conc, data = viabTab.aml[viabTab.aml$patientID==pat & viabTab.aml$drug=="Doxorubicine",], 
                       fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50")))
  
  # Create a new data frame for prediction
  new_data_aml_dox <- data.frame(conc = conc_seq$conc)
  
  # Create predictions
  new_data_aml_dox$viability <- predict(model_aml_dox, newdata = new_data_aml_dox)
  new_data_aml_dox$patientID <- pat
  
  # Store in list
  predictions_list[[pat]] <- new_data_aml_dox
}

# Combine all predictions
all_predictions_aml_dox <- bind_rows(predictions_list)

#plot drug response curves
#AraC
p03 <- all_predictions_aml_dox |> 
  merge(viabTab.aml[,c("patientID", "died", "OS")], by="patientID", all.x = TRUE) |> 
  filter(!is.na(died)) |> 
  ggplot(aes(x=conc, y=viability, color=patientID, linetype=died, group=patientID)) +
  geom_line()+
  scale_x_log10()+
  theme_classic()+
  #scale_color_manual(values=mycolors)+
  labs(title = "Doxorubicine in AML", y="Viability", x=bquote("Log Concentration ("*mu*"M)"), color = "PatientID")+
  theme(legend.position = "right", legend.title=element_text(size=8, face="bold"), 
        legend.text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1, color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), axis.line = element_line(size = 0.5), 
        axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
        plot.title = element_text(hjust=0.5, size=8, face="bold"), 
        legend.key.size = unit(0.5, 'cm'))
p03


#Fit sigmoid curve for each drug and each sample

sumData <- viabTab.aml |> 
  dplyr::rename(Drug = drug, viab = viability) |> 
  dplyr::select(patientID, sampleID, Drug, conc, viab)

testF <- function(m0, m1) {
  Fstat <- ((sum(residuals(m0)^2) - sum(residuals(m1)^2))/(m0$df.residual-m1$df.residual))/(sum(residuals(m1)^2)/m1$df.residual)
  1 - pf(Fstat, m0$df.residual-m1$df.residual, m1$df.residual)
}

calcAUC <- function(m1, minConc =0, maxConc=15, n=100) {
  concSeq <- data.frame(conc = seq(minConc, maxConc, length.out = n))
  valueConc <- data.frame(viab = predict(m1, concSeq),
                          conc = concSeq[,1])
  valueConc <- valueConc[order(valueConc$conc),]
  areaTotal <- 0
  for (i in seq(nrow(valueConc)-1)) {
    areaEach <- (valueConc$viab[i] + valueConc$viab[i+1])*(valueConc$conc[i+1]-valueConc$conc[i])*0.5
    areaTotal <- areaTotal + areaEach
  }
  areaTotal/(valueConc[nrow(valueConc),]$conc-valueConc[1,]$conc)
}

calcParm <- function(df_input) {
  df <- as.data.frame(df_input)
  df$conc <- as.numeric(df$conc)
  df$viab <- as.numeric(df$viab)
  
  formula_obj <- viab ~ conc
  
  out <- tryCatch({
    m1 <- eval(substitute(
      drm(formula_obj,
          data = df,
          fct = LL2.3u(),
          robust = "tukey",
          lowerl = c(0, 0, NA),
          upperl = c(4, 1, NA)),
      list(formula_obj = formula_obj, df = df)
    ))
    
    fitRes <- m1$coefficients
    m0 <- lm(viab ~ 1, data = df)
    pred_vals <- predict(m1)
    r2 <- cor(pred_vals, df$viab, use = "pairwise.complete.obs")^2
    auc <- calcAUC(m1)
    
    tibble(
      HS = fitRes[[1]],
      Einf = fitRes[[2]],
      IC50 = fitRes[[3]],
      pF = testF(m0, m1),
      R2 = r2,
      AUC = auc
    )
  }, error = function(e) {
    message("Fit failed for a group: ", e$message)
    tibble(
      HS = NA_real_,
      Einf = NA_real_,
      IC50 = NA_real_,
      pF = NA_real_,
      R2 = NA_real_,
      AUC = NA_real_
    )
  })
  
  return(out)
}

ic50Tab_aml <- sumData %>%
  group_by(patientID, Drug) %>% nest() %>%
  mutate(mm = map(data, calcParm)) %>%
  unnest(mm) %>%
  dplyr::select(-data) %>% ungroup()

#poor fit for doxo

#combine
aml_surv_viab <- merge(viabTab.aml, ic50Tab_aml, by.x=c("patientID", "drug"), by.y=c("patientID", "Drug"))

aml_surv_viab_cyt <- aml_surv_viab |> filter(drug=="Cytarabine")

#cox regression
#univariate

#univariate analysis for cyt response
parameters <- colnames(aml_surv_viab_cyt[,c(23:28)])
parameters <- append(parameters, c("viab.mean", "viab.auc"))
drugs <- c("Cytarabine")

resOS_uni_aml <- data.frame()

for (par in parameters) {
  for (drug in drugs) {
    res <- filter(aml_surv_viab_cyt, !is.na(HS), !is.na(Einf), !is.na(IC50), !is.na(pF), !is.na(R2), !is.na(AUC)) %>%
      do(com_uni(.[[par]], .$OS, .$died, TRUE))
    
    resOS_uni_row <- data.frame(
      parameter = par,
      cox = res,
      Drug = drug
    )
    resOS_uni_aml <- rbind(resOS_uni_aml, resOS_uni_row)
    
  }
}
resOS_uni_aml <- resOS_uni_aml |> 
  arrange(Drug, cox.p)

#plot
mycolors_cox_cyt <- setNames(c("#2171b5"), c("Cytarabine"))

SFig13h <- resOS_uni_aml |> 
  filter(!(parameter %in%c("viab.auc", "pF"))) |> 
  ggplot(aes(y=cox.HR, x=parameter, color=Drug))+
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = cox.lower, ymax = cox.higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size = 3,
            aes(y = cox.higher + 0.75, 
                label = ifelse(cox.p < 0.001, 
                               "italic(P)<0.001", 
                               paste0("italic(P)==", sprintf("%1.3f", cox.p)))),
            parse = TRUE)+
  theme_classic()+
  labs(y="Hazard ratio", x="")+
  coord_flip() + theme_classic() +
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none",
        plot.title = element_text(hjust=0.6, size=8, face="bold"))+
  scale_color_manual(values=mycolors_cox_cyt)+
  guides(color = guide_legend(override.aes = aes(label = "")))+
  ylim(0,3)+
  labs(title="Univariate Cox OS:\nCytarabine drug-response metrics", y="Hazard ratio", x="")+
  scale_x_discrete(labels = c(
    "AUC",
    expression(E[inf]),
    "Hill coefficient", 
    expression(IC[50]),  # subscript
    expression(R^2),      # superscript
    "Mean viability"))

SFig13h <- SFig13h + theme(plot.margin = margin(t=1, r = 1, b=2, unit = "cm"))
SFig13h


#kaplan-meier
aml_km <- aml_surv_viab_cyt |> 
  filter(conc_level == 1)

#maxstat
#OS
res.cut <- surv_cutpoint(aml_km, time = "OS", event = "died",
                         variables = "viab.mean")

summary(res.cut)
os_value <- summary(res.cut)$cutpoint
plot(res.cut, "viab.mean")
res.cat <- surv_categorize(res.cut)
head(res.cat)

fit <- survfit(Surv(OS, died) ~ viab.mean, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)

aml_km <- aml_km |> 
  mutate(
    cyt_resistance_os = ifelse(
      viab.mean < os_value, "Low mean viability", "High mean viability")
  )

os <- survfit2(Surv(OS, died) ~ cyt_resistance_os, data = aml_km)

SFig13i<- os |>
  ggsurvfit() +
  labs(x = "Years", y = "Overall survival", title = 
         "Cytarabine") + 
  #add_confidence_interval() +
  add_censor_mark(size = 1, alpha = 0.5) +
  #add_risktable(size = 3)+
  scale_ggsurvfit() +
  labs(x="Years")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), 
        axis.title = element_text(size=8), plot.title = element_text(hjust=0.5, size=8, face="bold"),
        legend.position = "bottom")+
  annotate("text", 
           x = 1,  # or specify a numeric value like x = 5
           y = 0.05, # position near bottom (or try 0.1, 0.15, etc.)
           label = glue::glue("Log-rank {survfit2_p(os)}"),
           hjust = 1,  # right-align (use 0 for left-align)
           vjust = 0,  # bottom-align (use 1 for top-align)
           size = 3)+
  scale_x_continuous(limits=c(0,3))
SFig13i <- SFig13i + theme(plot.margin = margin(t=1, l = 1, unit = "cm"))
SFig13i

#validate in SMART, select for cytarabine treatments (inverse correlation for )
smart_cyt <- df_chemo |> 
  filter(diagnosis == "AML", name == "Cytarabine") |> 
  merge(df_invivo, by="patientID") |> 
  filter(grepl("DA-Induction", treatment_spec), !(Relationship_treatment_death %in% c("Related", "Possibly related")))

smart_cyt |> 
  ggplot(aes(x=factor(Death), y=normVal))+
  geom_boxplot(outliers=FALSE)+
  geom_point()+
  stat_compare_means(method="t.test")

#Fit sigmoid curve for each sample

sumData <- smart_cyt |> 
  dplyr::rename(Drug = name, viab = normVal, conc = concentration) |> 
  dplyr::select(patientID, Drug, conc, viab)

testF <- function(m0, m1) {
  Fstat <- ((sum(residuals(m0)^2) - sum(residuals(m1)^2))/(m0$df.residual-m1$df.residual))/(sum(residuals(m1)^2)/m1$df.residual)
  1 - pf(Fstat, m0$df.residual-m1$df.residual, m1$df.residual)
}

calcAUC <- function(m1, minConc =0, maxConc=15, n=100) {
  concSeq <- data.frame(conc = seq(minConc, maxConc, length.out = n))
  valueConc <- data.frame(viab = predict(m1, concSeq),
                          conc = concSeq[,1])
  valueConc <- valueConc[order(valueConc$conc),]
  areaTotal <- 0
  for (i in seq(nrow(valueConc)-1)) {
    areaEach <- (valueConc$viab[i] + valueConc$viab[i+1])*(valueConc$conc[i+1]-valueConc$conc[i])*0.5
    areaTotal <- areaTotal + areaEach
  }
  areaTotal/(valueConc[nrow(valueConc),]$conc-valueConc[1,]$conc)
}

calcParm <- function(df_input) {
  df <- as.data.frame(df_input)
  df$conc <- as.numeric(df$conc)
  df$viab <- as.numeric(df$viab)
  
  formula_obj <- viab ~ conc
  
  out <- tryCatch({
    m1 <- eval(substitute(
      drm(formula_obj,
          data = df,
          fct = LL2.3u(),
          robust = "tukey",
          lowerl = c(0, 0, NA),
          upperl = c(4, 1, NA)),
      list(formula_obj = formula_obj, df = df)
    ))
    
    fitRes <- m1$coefficients
    m0 <- lm(viab ~ 1, data = df)
    pred_vals <- predict(m1)
    r2 <- cor(pred_vals, df$viab, use = "pairwise.complete.obs")^2
    auc <- calcAUC(m1)
    
    tibble(
      HS = fitRes[[1]],
      Einf = fitRes[[2]],
      IC50 = fitRes[[3]],
      pF = testF(m0, m1),
      R2 = r2,
      AUC = auc
    )
  }, error = function(e) {
    message("Fit failed for a group: ", e$message)
    tibble(
      HS = NA_real_,
      Einf = NA_real_,
      IC50 = NA_real_,
      pF = NA_real_,
      R2 = NA_real_,
      AUC = NA_real_
    )
  })
  
  return(out)
}

ic50Tab_smart <- sumData %>%
  group_by(patientID, Drug) %>% nest() %>%
  mutate(mm = map(data, calcParm)) %>%
  unnest(mm) %>%
  dplyr::select(-data) %>% ungroup()


#combine
smart_surv_viab <- merge(smart_cyt, ic50Tab_smart, by.x=c("patientID", "name"), by.y=c("patientID", "Drug"))

#cox regression
#univariate

#univariate analysis for cyt response
parameters <- colnames(smart_surv_viab[,c(63:68)])
parameters <- append(parameters, c("AUC.x", "IC50.x", ""))

drugs <- c("Cytarabine")

resOS_uni_smart <- data.frame()

for (par in parameters) {
  for (drug in drugs) {
    res <- filter(smart_surv_viab, !is.na(HS), !is.na(Einf), !is.na(IC50.x), !is.na(pF), !is.na(R2), !is.na(AUC.x)) %>%
      do(com_uni(.[[par]], .$time_OS, .$Death, TRUE))
    
    resOS_uni_row <- data.frame(
      parameter = par,
      cox = res,
      Drug = drug
    )
    resOS_uni_smart <- rbind(resOS_uni_smart, resOS_uni_row)
    
  }
}
resOS_uni_smart <- resOS_uni_smart |> 
  arrange(Drug, cox.p)

#plot
FigS11y <- resOS_uni_smart |> 
  ggplot(aes(y=cox.HR, x=parameter))+
  geom_hline(yintercept=1, color="grey", linetype="dashed")+
  geom_point()+
  geom_errorbar(aes(ymin=cox.lower, ymax=cox.higher))+
  coord_flip()+
  facet_wrap(~Drug)+
  theme_classic()+
  labs(title="Univariate Cox OS\n(Cytarabine drug-response metrics)", y="HR", x="Parameter")+
  theme(axis.text.x = element_text(color="black", size=8), 
        axis.text.y = element_text(color="black", size=8), 
        axis.title.y = element_text(size=8), axis.title.x = element_text(size=8), 
        plot.title = element_text(hjust=0.5, size=8, face="bold"), 
        legend.key.size = unit(0.5, 'cm'))
FigS11y

