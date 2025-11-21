library(tidyverse)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(maxstat)
library(survminer)

#load data
usedSamples_all <- read.csv("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/usedSamples_all.csv", sep=";")

viability_table <- read.csv2("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/viability_table.csv")

load("~/Dropbox/Medizin/Forschung/Aktuelle Forschungsprojekte/CLL CPS1000/raw data/patmeta_210324.RData")

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
combined <- merge(viabTab.mean, geneBack, by.x="patientID", by.y="Patient.ID") |> 
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


#multivariate
com <- function(response, ighv, tp53, del17p, time, endpoint, scale =FALSE) {  
  
  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
  }
  
  tryCatch({
    surv <- coxph(Surv(time, endpoint) ~ response+ighv+tp53+del17p) 
    resTab <- tibble(variable = c("response", "ighv", "tp53", "del17p"),
                     p = summary(surv)[[7]][,5], 
                     HR = summary(surv)[[7]][,2], 
                     lower = summary(surv)[[8]][,3], 
                     higher = summary(surv)[[8]][,4])
  }, error = function(err) {
    resTab <- tibble(variable = c("response", "ighv", "tp53", "del17p"),
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
testTab <- left_join(viabTab.mean, survT, by = c("sampleID","patientID"))

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

FigS11d <- plotTab_uni |> 
  mutate(IGHV.status = ifelse(IGHV.status == 1, "M-CLL", "U-CLL")) |> 
           ggplot(aes(x=drug, y = HR, col = Endpoint, dodge = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size=3,
            aes(y = higher + 0.5, label = sprintf("p=%1.3f",p))) +
  facet_wrap(~IGHV.status, scales = "free") + ylim(0,4) +
  coord_flip() +
  labs(y="Hazard ratio", x="")+ 
  theme_bw()+
  theme(axis.text = element_text(color="black", size=8),
                                      legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8))
FigS11d


#multivariate
#for OS
resOS <- filter(testTab2, !is.na(OS), !is.na(IGHV.status), !is.na(TP53), !is.na(del17p)) %>%
  group_by(drug) %>%
  do(com(.$viab, .$IGHV.status, .$TP53, .$del17p, .$OS, .$died, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "OS")

#for TTT
resTTT <- filter(testTab2, !is.na(TTT), !is.na(IGHV.status), !is.na(TP53), !is.na(del17p)) %>%
  group_by(drug) %>%
  do(com(.$viab, .$IGHV.status, .$TP53, .$del17p,.$TTT, .$treatedAfter, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(Endpoint = "TTT")

#Plot Hazard ratios and p values
plotTab <- bind_rows(resOS, resTTT) %>%
  filter(drug %in% unique(filter(.,p.adj < 0.05)$drug))

#filter drugs with significant impact of response
drug_list <- plotTab |> filter(variable == "response", p.adj < 0.1) |> 
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

#plot factors for ONO
Fig7c <- plotTab |> 
  filter(drug == "ONO.4059") |> 
  ggplot(aes(x=variable, y = HR, col = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(position = position_dodge(width =0.75), 
                aes(ymin = lower, ymax = higher), width = 0.3) + 
  geom_text(position = position_dodge(width =0.75), size = 3,
            aes(y = higher + 0.5, label = sprintf("p=%1.3f",p))) +
  #facet_wrap(~Endpoint)+
  labs(y="Hazard ratio", x="")+
  coord_flip() + theme_classic() +
  theme(axis.text = element_text(color="black", size=8),
        legend.title = element_text(size=8, face="bold"), 
        legend.text=element_text(size=8), axis.title = element_text(size=8))+
  scale_x_discrete(labels = c("tp53" = "TP53",
                              "response" = "ONO-4059\nresponse",
                              "ighv" = "IGHV",
                              "del17p" = "del(17p)"))+
  ylim(0,6)

Fig7c + theme(plot.margin = margin(t=1, b=1, l=1, r=1, unit = "cm"))
Fig7c

#maxstat
testTab3 <- testTab2 |> 
  filter(!is.na(IGHV.status), drug == "ONO.4059") 

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
  labs(x = "Years", y = "Overall survival probability", title = 
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
  scale_color_manual(labels=c("U-CLL, low ONO-4059 viability", "U-CLL, high ONO-4059 viability", 
                     "M-CLL, low ONO-4059 viability", "M-CLL, high ONO-4059 viability"), 
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
  labs(x = "Years", y = "Treatment-free survival probability", title = 
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
  scale_color_manual(labels=c("U-CLL, low ONO-4059 viability", "U-CLL, high ONO-4059 viability", 
                              "M-CLL, low ONO-4059 viability", "M-CLL, high ONO-4059 viability"), 
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

Fig7d <- Fig7d1 + Fig7d2 + plot_layout(ncol=2, guides = 'collect') & theme(legend.position='bottom', 
                                                                           plot.margin = margin(t=1, b=1, l=1, r=1, unit = "cm"))
Fig7d
