library(tidyverse)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(maxstat)
library(survminer)
library(readxl)

#load data
usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/patmeta_210324.RData")

dob <- read_excel("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/clinical annotation_HD/clinical_data_cps1000_from_clinical.xlsx") |> 
  dplyr::select("Patient ID", "Date of birth", diagnosis) |> 
  filter(diagnosis == "CLL")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/survival_190516.RData")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/20250625_timeline.RData")

#calculate age at sampling
dob <- dob %>%
  mutate(dob = dmy(`Date of birth`))

survT <- survT |> 
  dplyr::select(patID, sampleID, sampleDate)

age <- merge(dob, survT, by.x="Patient ID", by.y="patID") |> 
  dplyr::select("Patient ID", dob, sampleID, sampleDate) |> 
  mutate(age = as.numeric(difftime(as.POSIXct(sampleDate), as.POSIXct(dob), units = "days"))/365) |> 
  rename(patientID = "Patient ID")

#load table with disease entities  
diseases <- usedSamples_all |> 
  arrange(patientID) |> 
  dplyr::select (patientID, sampleID, diagnosis) 
diseases

# merge the tables
merged_table <- cbind(diseases, viability_table) |> 
  relocate(diagnosis, .before = X10058.F4_1)
merged_table

#transfrom in long data frame
merged_table_long <- merged_table |>
  pivot_longer(cols = c(4:318), names_to = "drug", values_to = "viability")
merged_table_long

merged_table_long$drug <- str_sub(merged_table_long$drug, end = -3)
merged_table_long

#create means, select CLL
viabTab.mean <- merged_table_long |> 
  group_by(patientID, sampleID ,diagnosis, drug) |> 
  summarize(viab = mean(viability)) |> 
  filter(diagnosis == "CLL")

# 1. untreated
#overview of treatment conditions
treatmentTab |> 
  filter(sampleID %in% viabTab.mean$sampleID) |> 
  group_by(treatStatus) |> 
  summarize(n=n())

#select only untreated samples
untreated <- treatmentTab |> 
  dplyr::select(Patient.ID, sampleID, treatStatus) |> 
  filter(treatStatus %in% c("untreated"))

viabTab.mean.untreated <- viabTab.mean |> 
  filter(sampleID %in% untreated$sampleID)

length(unique(viabTab.mean.untreated$patientID))

#merge with genomic variantions
#genetic background
geneBack <- patMeta[match(unique(viabTab.mean$patientID), patMeta$Patient.ID),] %>%
  dplyr::select(-HIPO.ID, -diagnosis, -project, -date.of.diagnosis, -treatment, -date.of.first.treatment,
                -gender, -Methylation_Cluster) %>%
  mutate(IGHV.status = ifelse(is.na(IGHV.status), NA, 
                              ifelse(IGHV.status == "M",1,0))) %>%
  data.frame() %>% column_to_rownames("Patient.ID") %>%
  mutate_all(as.character) %>% mutate_all(as.integer) 

#only select variants that have at least 5 mutated cases
geneBack <- geneBack[,colSums((!is.na(geneBack) & geneBack == 1)) >=5]
geneBack$Patient.ID <- unique(viabTab.mean$patientID)

geneTab <- gather(geneBack, key = "gene", value = "status", -Patient.ID)

#combine
combined <- merge(viabTab.mean.untreated, geneBack, by.x="patientID", by.y="Patient.ID") |> 
  dplyr::select(-sampleID, -diagnosis, -viab)

####

#Function for cox regression
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
                     higher = summary(surv)[[8]][,4])
  }, error = function(err) {
    resTab <- tibble(p = NA, 
                     HR = NA, 
                     lower = NA, 
                     higher = NA)
  })
  
}

#include age?

#multivariate
com <- function(response, ighv, TP53_del17p, age, time, endpoint, scale =FALSE) {  
  
  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
    age <- (age - mean(age, na.rm = TRUE))/sd(age, na.rm=TRUE)
  }
  
  tryCatch({
    surv <- coxph(Surv(time, endpoint) ~ response+ighv+TP53_del17p+age) 
    resTab <- tibble(variable = c("response", "ighv", "TP53_del17p", "age"),
                     p = summary(surv)[[7]][,5], 
                     HR = summary(surv)[[7]][,2], 
                     lower = summary(surv)[[8]][,3], 
                     higher = summary(surv)[[8]][,4])
  }, error = function(err) {
    resTab <- tibble(variable = c("response", "ighv", "TP53_del17p", "age"),
                     p = NA, 
                     HR = NA, 
                     lower = NA, 
                     higher = NA)
  })
  
}

#Perform cox regression for each drug and clinical feature
#load survival table
load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/survival_190516.RData")
survT <- dplyr::rename(survT, patientID = patID)

#combine drug response and clinical table
testTab <- left_join(viabTab.mean.untreated, survT, by = c("sampleID","patientID"))

#combine with genetic data
testTab2 <- merge(testTab, combined, by = c("patientID","drug"))



#univariate analysis for drugs only, separated by IGHV
#for OS
resOS_uni <- filter(testTab2, !is.na(OS), !is.na(IGHV.status)) %>%
  group_by(drug, IGHV.status) %>%
  do(com_uni(.$viab, .$OS, .$died, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "OS")

#for TTT
resTTT_uni <- filter(testTab2, !is.na(TTT), !is.na(IGHV.status)) %>%
  group_by(drug, IGHV.status) %>%
  do(com_uni(.$viab, .$TTT, .$treatedAfter, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "TTT")

#Plot Hazard ratios and p values
plotTab_uni <- bind_rows(resOS_uni, resTTT_uni) %>%
  filter(drug %in% unique(filter(.,p.adj < 0.05)$drug))

plotTab_uni |> 
  mutate(IGHV.status = ifelse(IGHV.status == 1, "M-CLL", "U-CLL")) |> 
           ggplot(aes(x=drug, y = HR, col = Endpoint, dodge = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size=3,
            aes(y = higher + 0.5, label = sprintf("p=%1.3f",p))) +
  facet_wrap(~IGHV.status, scales = "free") + 
  ylim(0,7) +
  coord_flip() +
  labs(title= "Untreated CLL", y="Hazard ratio", x="")+ 
  theme_bw()+
  theme(axis.text = element_text(color="black", size=8),
                                      legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8),
        plot.title = element_text(size=8, color="black", face="bold", hjust=0.5))

# 2. select all samples
#select samples
viabTab.mean.all <- viabTab.mean

length(unique(viabTab.mean.all$patientID))

#merge with genomic variantions
#combine
combined <- merge(viabTab.mean.all, geneBack, by.x="patientID", by.y="Patient.ID") |> 
  dplyr::select(-sampleID, -diagnosis, -viab)

####

#Perform cox regression for each drug and clinical feature
#combine drug response and clinical table
testTab <- left_join(viabTab.mean.all, survT, by = c("sampleID","patientID"))

#combine with genetic data
testTab2 <- merge(testTab, combined, by = c("patientID","drug"))

#summary
length(unique(testTab2$patientID))

testTab2 |> 
  group_by(treatedAfter) |> 
  summarize(n=n()/63)

#univariate analysis for drugs only, separated by IGHV
#for OS
resOS_uni <- filter(testTab2, !is.na(OS), !is.na(IGHV.status)) %>%
  group_by(drug, IGHV.status) %>%
  do(com_uni(.$viab, .$OS, .$died, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "OS")

#for TTT
resTTT_uni <- filter(testTab2, !is.na(TTT), !is.na(IGHV.status)) %>%
  group_by(drug, IGHV.status) %>%
  do(com_uni(.$viab, .$TTT, .$treatedAfter, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "TTT")

#Plot Hazard ratios and p values
plotTab_uni <- bind_rows(resOS_uni, resTTT_uni) %>%
  filter(drug %in% unique(filter(.,p.adj < 0.05)$drug))

#no sign. p.adj

FigS11d <- plotTab_uni |> 
  mutate(IGHV.status = ifelse(IGHV.status == 1, "M-CLL", "U-CLL")) |> 
  ggplot(aes(x=drug, y = HR, col = Endpoint, dodge = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size=3,
            aes(y = higher + 0.5, label = sprintf("p=%1.3f",p))) +
  facet_wrap(~IGHV.status, scales = "free") + 
  ylim(0,5) +
  coord_flip() +
  labs(y="Hazard ratio", x="")+ 
  theme_bw()+
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8))
FigS11d

#multivariate
#for OS
resOS <- testTab2 |> 
  mutate(TP53_del17p = ifelse(TP53==1|del17p==1, 1,0)) |> 
  merge(age, by=c("patientID", "sampleID"), .before = 3) |> 
  filter(!is.na(OS), !is.na(IGHV.status), !is.na(TP53_del17p), !is.na(age)) %>%
  group_by(drug) %>%
  do(com(.$viab, .$IGHV.status, .$TP53_del17p, .$age, .$OS, .$died, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "OS")

#for TTT
resTTT <- testTab2 |> 
  mutate(TP53_del17p = ifelse(TP53==1|del17p==1, 1,0)) |> 
  merge(age, by=c("patientID", "sampleID"), .before = 3) |> 
  filter(!is.na(TTT), !is.na(IGHV.status), !is.na(TP53_del17p), !is.na(age)) %>%
  group_by(drug) %>%
  do(com(.$viab, .$IGHV.status, .$TP53_del17p,.$age,.$TTT, .$treatedAfter, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "TTT")

#Plot Hazard ratios and p values
plotTab <- bind_rows(resOS, resTTT) %>%
  filter(drug %in% unique(filter(.,p.adj < 0.05)$drug))

#filter drugs with significant impact of response
drug_list <- plotTab |> filter(variable == "response", p.adj < 0.05) |> 
  pull(drug)

plotTab |> 
  filter(drug %in% drug_list, variable == "response") |> 
ggplot(aes(x=drug, y = HR, col = Endpoint, dodge = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width =0.8), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.8),
            aes(y = higher + 0.25, label = sprintf("p=%1.2f",p))) +
  ylim(-0.5,4) + 
  coord_flip() + theme_bw()+
  labs(title = "Drug responses with sign. effect on OS/TT\nin Cox proportional hazards model (IGHV, TP53, del17p)")

#plot HR
mycolors_cox <- setNames(c("#2171b5", "#cb181d"), c("TTT", "OS"))

Fig7e <- plotTab |> 
  filter(drug == "Ibrutinib") |> 
  ggplot(aes(x=variable, y = HR, col = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size = 3,
            aes(y = higher + 0.6, 
                label = paste0("italic(P)==", sprintf("%1.3f", p))),
            parse = TRUE) +
  #facet_wrap(~Endpoint)+
  labs(y="Hazard ratio", x="")+
  coord_flip() + theme_classic() +
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8))+
  scale_x_discrete(labels = c("TP53_del17p" = "TP53/del(17p)",
                              "response" = "Ex vivo\nresistance to\nibrutinib",
                              "ighv" = "IGHV", 
                              "age" = "Age at\nsampling"))+
  scale_color_manual(values=mycolors_cox)+
  guides(color = guide_legend(override.aes = aes(label = "")))+
  ylim(0,6)

Fig7e + theme(plot.margin = margin(t=3, b=1, l=1, r=1, unit = "cm"))
Fig7e

#maxstat
testTab3 <- testTab2 |> 
  filter(!is.na(IGHV.status), drug == "Ibrutinib") 

#OS
res.cut <- surv_cutpoint(testTab3, time = "OS", event = "died",
                         variables = c("viab"))

summary(res.cut)
os_value <- summary(res.cut)$cutpoint
plot(res.cut, "viab")
res.cat <- surv_categorize(res.cut)
head(res.cat)

fit <- survfit(Surv(OS, died) ~ viab, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)

#TTT
res.cut <- surv_cutpoint(testTab3, time = "TTT", event = "treatment",
                         variables = c("viab"))

summary(res.cut)
ttt_value <- summary(res.cut)$cutpoint
plot(res.cut, "viab")
res.cat <- surv_categorize(res.cut)
head(res.cat)

fit <- survfit(Surv(TTT, treatment) ~ viab, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)

#Kaplan-meier curve
testTab3 <- testTab3 |> 
  mutate(
    ono_resistance_os = ifelse(
    viab < os_value, 0, 1),
    ono_resistance_ttt = ifelse(
      viab < ttt_value, 0, 1)
  )
  
os <- survfit2(Surv(OS, died) ~ ono_resistance_os+IGHV.status, data = testTab3)

Fig7d1 <- os |>
  ggsurvfit() +
  labs(x = "Years", y = "Overall survival", title = 
         "Overall survival") + 
  #add_confidence_interval() +
  add_censor_mark(size = 1, alpha = 0.5) +
  #add_risktable(size = 3)+
  scale_ggsurvfit() +
  labs(x="Years from sampling")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8),
              legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), 
        axis.title = element_text(size=8), plot.title = element_text(hjust=0.5, size=8, face="bold"),
        legend.position = "bottom")+
  scale_color_manual(labels=c("U-CLL, low ibrutinib viability", "U-CLL, high ibrutinib viability", 
                     "M-CLL, low ibrutinib viability", "M-CLL, high ibrutinib viability"), 
                     values = c("#CA054D", "#B96D40", "#A4D4B4", "#3B1C32"))+
  annotate("text", 
           x = 2,  # or specify a numeric value like x = 5
           y = 0.05, # position near bottom (or try 0.1, 0.15, etc.)
           label = glue::glue("Log-rank {survfit2_p(os)}"),
           hjust = 1,  # right-align (use 0 for left-align)
           vjust = 0,  # bottom-align (use 1 for top-align)
           size = 3)+
  guides(color = guide_legend(ncol = 2))
Fig7d1


ttt <- survfit2(Surv(TTT) ~ ono_resistance_ttt+IGHV.status, data = testTab3)

Fig7d2 <- ttt |> 
  ggsurvfit() +
  labs(x = "Years", y = "Treatment-free survival", title = 
         "Time to treatment") + 
  #add_confidence_interval() +
  add_censor_mark(size = 1, alpha = 0.5) +
  #add_risktable() +
  scale_ggsurvfit() +
  labs(x="Years from sampling")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=8), 
        axis.title = element_text(size=8), plot.title = element_text(hjust=0.5, size=8, face="bold"),
        legend.position = "bottom")+
  scale_color_manual(labels=c("U-CLL, low ibrutinib viability", "U-CLL, high ibrutinib viability", 
                              "M-CLL, low ibrutinib viability", "M-CLL, high ibrutinib viability"), 
                     values = c("#CA054D", "#B96D40", "#A4D4B4", "#3B1C32"))+
  annotate("text", 
           x = 2,  # or specify a numeric value like x = 5
           y = 0.05, # position near bottom (or try 0.1, 0.15, etc.)
           label = glue::glue("Log-rank {survfit2_p(os)}"),
           hjust = 1,  # right-align (use 0 for left-align)
           vjust = 0,  # bottom-align (use 1 for top-align)
           size = 3)+
  guides(color = guide_legend(ncol = 2))
Fig7d2

Fig7f <- (Fig7d1/Fig7d2) + plot_layout(ncol=1, guides = 'collect') & theme(legend.position='bottom', 
                                                                           plot.margin = margin(t=1, b=1, l=1, r=1, unit = "cm"))
Fig7f

#validation in EMBL2016
load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/EMBL2016/newEMBL_processed_20220506.RData")
load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/EMBL2016/patmeta_210324.RData")

patMeta[patMeta$Patient.ID %in% "P1499",]$diagnosis <- "HCL-v"

##define some useful function
absmax <- function(x) { x[which.max( abs(x) )]}

# Sample filtering
#1. Remove counter screen and CV plates 
screenData <- filter(emblNew, patID != "ATP", !is.na(name))

#2. Remove low-quality samples and problematic drugs based on previous QA analysis
screenData <- filter(screenData, ! lowQuality)

#3. Use incubation effect adjusted viabilities
screenData <- mutate(screenData, viab = normVal.sigm, viab.auc = normVal_auc.sigm)

#4. Remove the two lowest concentrations
#screenData <- filter(screenData, concIndex %in% seq(1,9))

#5. Use AUC (mean) to summarise drug effect across 9 concentrations
screenData <- dplyr::rename(screenData, viab.mp = normVal_auc.sigm, patientID = patID) #AUC of edge effect corrected values

#6. Show drugs that were not tested for all samples
nTotal <- length(unique(screenData$patID))

#drug that are not screened in all samples
drugTab <- distinct(screenData, patientID, name) %>% mutate_all(as.character) %>%
  group_by(name) %>% summarise(n = length(patientID)) %>% filter(n<nTotal) %>%
  arrange(n)

screenData <- screenData %>% filter(!str_detect(name, "PRO|CHEMBL"))

#Annoate target class
tarAnno <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/EMBL2016/targetAnnotation_all.csv")
screenData <- mutate(screenData, class = tarAnno[match(emblObjID, tarAnno$drugID),]$Class,
                     target = tarAnno[match(emblObjID, tarAnno$drugID),]$target)

# Format table
screenData <- dplyr::select(screenData, well, val, patientID, batch, 
                            sampleID, drugID, conc, concIndex, viab, viab.auc, name, diagnosis) %>%
  mutate(patientID = as.character(patientID), drugID = as.character(drugID)) %>%
  dplyr::rename(Drug = name) %>% 
  mutate(concIndex = as.factor(as.character(concIndex)),
         diagnosis = patMeta[match(patientID, patMeta$Patient.ID),]$diagnosis)

#exclude common patients with CPS1000
cps_pat <- intersect(screenData$patientID, viabTab.all.mean$patientID)

screenData <- screenData |> 
  filter(!(patientID %in% cps_pat), diagnosis == "CLL")

#patient samples
length(unique(screenData$patientID))

#drugs
length(unique(screenData$Drug))

#merge with genomic variantions
#genetic background
geneBack_embl <- patMeta[match(unique(screenData$patientID), patMeta$Patient.ID),] %>%
  dplyr::select(-HIPO.ID, -diagnosis, -project, -date.of.diagnosis, -treatment, -date.of.first.treatment,
                -gender, -Methylation_Cluster) %>%
  mutate(IGHV.status = ifelse(is.na(IGHV.status), NA, 
                              ifelse(IGHV.status == "M",1,0))) %>%
  data.frame() %>% column_to_rownames("Patient.ID") %>%
  mutate_all(as.character) %>% mutate_all(as.integer) 

#only select variants that have at least 5 mutated cases
geneBack_embl <- geneBack_embl[,colSums((!is.na(geneBack_embl) & geneBack_embl == 1)) >=5]
geneBack_embl$Patient.ID <- unique(screenData$patientID)

geneTab_embl <- gather(geneBack_embl, key = "gene", value = "status", -Patient.ID)

#combine
combined_embl <- screenData |> 
  dplyr::select(patientID, sampleID, conc, concIndex, viab, viab.auc, Drug) |> 
  filter(Drug %in% c("Ibrutinib")) |> 
  group_by(patientID, sampleID, Drug) |> 
  summarize(viab.mean = mean(viab)) |> 
  merge(geneBack_embl[,c("Patient.ID", "IGHV.status", "del17p", "TP53")], by.x="patientID", by.y="Patient.ID")

#merge with treatment and survival data
surv.all.embl <- combined_embl |> 
  merge(treatmentTab[,c("Patient.ID", "sampleID", "sampleDate", "treatStatus")], by.x=c("patientID", "sampleID"), by.y=c("Patient.ID", "sampleID")) |> 
  merge(survT, by=c("patientID", "sampleID")) |> 
  filter(treatStatus %in% c("untreated", "treatedBefore"))

#summary
length(unique(surv.all.embl$patientID))

surv.all.embl |> 
  group_by(treatStatus) |> 
  summarize(n=n())

#univariate analysis
#for OS
resOS_uni_embl <- filter(surv.all.embl, !is.na(OS)
                         #, !is.na(IGHV.status)
                         ) %>%
  group_by(Drug
           #, IGHV.status
           ) %>%
  do(com_uni(.$viab.mean, .$OS, .$died, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "OS", Predictor = "viab.mean")

resOS_uni_embl_gen <- data.frame()
genetics <- c("IGHV.status", "del17p", "TP53")

for (gen in genetics) {
  
res <- filter(surv.all.embl, !is.na(OS)) %>%
  group_by(Drug) %>%
  do(com_uni(.[[gen]], .$OS, .$died, FALSE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "OS", Predictor = gen)

resOS_uni_embl_gen <- rbind(resOS_uni_embl_gen, res)
}

resOS_uni_embl <- rbind(resOS_uni_embl, resOS_uni_embl_gen)
resOS_uni_embl

#for TTT
resTTT_uni_embl <- filter(surv.all.embl, !is.na(TTT), !is.na(IGHV.status)) %>%
  group_by(Drug
           #, IGHV.status
           ) %>%
  do(com_uni(.$viab.mean, .$TTT, .$treatedAfter, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "TTT", Predictor = "viab.mean")

resTTT_uni_embl_gen <- data.frame()
genetics <- c("IGHV.status", "del17p", "TP53")

for (gen in genetics) {
  
  res <- filter(surv.all.embl, !is.na(TTT)) %>%
    group_by(Drug) %>%
    do(com_uni(.[[gen]], .$TTT, .$treatedAfter, FALSE)) %>% ungroup() %>%
    arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
    mutate(Endpoint = "TTT", Predictor = gen)
  
  resTTT_uni_embl_gen <- rbind(resTTT_uni_embl_gen, res)
}

resTTT_uni_embl <- rbind(resTTT_uni_embl, resTTT_uni_embl_gen)
resTTT_uni_embl

#Plot Hazard ratios and p values
plotTab_uni_embl <- bind_rows(resOS_uni_embl, resTTT_uni_embl)

FigS11h <- plotTab_uni_embl |> 
  ggplot(aes(x=Predictor, y = HR, col = Endpoint, dodge = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size=3,
            aes(y = higher + 1.5, label = sprintf("p=%1.3f",p))) +
  coord_flip() +
  labs(title= "Univariate Cox (validation: EMBL2016)", y="Hazard ratio", x="")+ 
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8),
        plot.title = element_text(size=8, color="black", face="bold", hjust=0.5))+
  scale_x_discrete(labels=c("del(17p)", "IGHV", "TP53", "Response to\nIbrutinib"))+
  ylim(0,10)
FigS11h

#multivariate
#for OS
resOS_embl <- filter(surv.all.embl, !is.na(OS), !is.na(IGHV.status), !is.na(TP53), !is.na(del17p)) %>%
  group_by(Drug) %>%
  do(com(.$viab.mean, .$IGHV.status, .$TP53, .$del17p, .$OS, .$died, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "OS")

#for TTT
resTTT_embl <- filter(surv.all.embl, !is.na(TTT), !is.na(IGHV.status), !is.na(TP53), !is.na(del17p)) %>%
  group_by(Drug) %>%
  do(com(.$viab.mean, .$IGHV.status, .$TP53, .$del17p,.$TTT, .$treatedAfter, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "TTT")

plotTab_multi_embl <- rbind(resOS_embl, resTTT_embl)

#plot factors
FigS11i <- plotTab_multi_embl |> 
  ggplot(aes(x=variable, y = HR, col = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size = 3,
            aes(y = higher + 1.5, label = sprintf("p=%1.3f",p))) +
  labs(title="Multivariate Cox (validation: EMBL2016)", y="Hazard ratio", x="")+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8),
        plot.title = element_text(size=8, color="black", face="bold", hjust=0.5))+
  scale_x_discrete(labels=c("del(17p)", "IGHV", "Response to\nIbrutinib", "TP53"))+
  ylim(0,10)
FigS11i


