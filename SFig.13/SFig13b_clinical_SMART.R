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
library(readxl)
library(ggpmisc)

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

#load clinical responses from trials
ORR_AML <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/clinical_responses_AML_CLL_20250707.xlsx", sheet = "AML")

ORR_AML <- ORR_AML[,1:7]

#filter for drugs with clinical data
ORR_AML <- ORR_AML |> 
  dplyr::rename(ORR = "ORR%", CR = "CR%") |> 
  dplyr::filter(ORR != "NA")
ORR_AML$ORR <- as.numeric(ORR_AML$ORR)
ORR_AML$CR <- as.numeric(ORR_AML$CR)  

drug_list <- unique(ORR_AML$drug)

#load drug class
category <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/analysis/drugCategory_modified.csv", sep=";") |> dplyr::select(drug, Category)

#filter for drugs in SMART
smart <- data_all |> 
  filter(Drug %in% drug_list, diagnosis == "AML")

#make plot: mean viab
smart_means <- smart |> 
  group_by(Drug) |> 
  summarize(viab_mean = mean(Drug_response_mean, na.rm=TRUE),
            viab_auc = mean(AUC, na.rm=TRUE)) |> 
  merge(ORR_AML, by.x="Drug", by.y="drug") |> 
  merge(category, by.x="Drug", by.y="drug")

mycolors <- c("Chemo" = "gray50", 
              "Kinase" = "#6A3D9A", 
              "Other" = "#CAB2D6")

smart_means$nr_pat <- as.numeric(smart_means$nr_pat)

aml_orr <- smart_means |> 
  ggplot(aes(x=viab_mean, y=ORR)) +
  geom_point(aes(size=nr_pat, color=Category)) +
  geom_text_repel(aes(label=Drug), size=3, color="black",
                  box.padding = 0.5)+
  stat_correlation(mapping = use_label("R", "R2", "P", "n"), label.x = "left", size=3)+
  theme_bw()+
  labs(title = "ORR ~ Mean viability", x = "Mean viability", y = "ORR", color="Drug type", size="Nr. of patients\nin clinical trial")+
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), axis.title = element_text(size=8), plot.title = element_text(hjust=0.5, size=8, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', 
                                    fill = NA, 
                                    size = 1))+
  scale_color_manual(values=mycolors, labels = c("Kinase inhibitor", "Chemotherapy", "Other"))+
  scale_size_continuous(
    breaks=c(50,100,150,200),
    labels =c(50,100,150,200))+
  ylim(0,100) 

aml_crr <- smart_means |> 
  ggplot(aes(x=viab_mean, y=CR)) +
  geom_point(aes(size=nr_pat, color=Category)) +
  geom_text_repel(aes(label=Drug), size=3, color="black",
                  box.padding = 0.5)+
  stat_correlation(mapping = use_label("R", "R2", "P", "n"), label.x = "left", size=3)+
  theme_bw()+
  labs(title = "CRR ~ Mean viability", x = "Mean viability", y = "CRR", color="Drug type", size="Nr. of patients\nin clinical trial")+
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), axis.title = element_text(size=8), plot.title = element_text(hjust=0.5, size=8, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', 
                                    fill = NA, 
                                    size = 1))+
  scale_color_manual(values=mycolors, labels = c("Kinase inhibitor", "Chemotherapy", "Other"))+
  scale_size_continuous(
    breaks=c(50,100,150,200),
    labels =c(50,100,150,200))+
  ylim(0,100)

SFig13b <- aml_orr+aml_crr + plot_annotation("AML (Liebers et al. 2023)", theme = theme(
  plot.margin = margin(t=1, b=1, l=1, unit = "cm"),
  plot.title = element_text(hjust=0.05, size=8, face="bold"))) +
  plot_layout(guides = "collect")
SFig13b

#...

#correlation of cytarabine parameters with CR rate
smart_cyt <- smart |> 
  filter(diagnosis == "AML", Drug == "Cytarabine")

test <- screenData |> 
  filter(diagnosis == "AML", Drug == "Cytarabine", treatment_spec == "DA-Induction")


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