##################################################################################################

### MEDICATION AND MULTIMORBIDITY ADDITIONAL ANALYSES TESTS

##################################################################################################

# srun -p interactive --pty bash

covs <- read.csv("covariates_additional.csv")
med <- covs[complete.cases(covs$Med),]
med$bp_linked <- as.character(med$bp_linked)
med$db_linked <- as.character(med$db_linked)
med$bp_linked[is.na(med$bp_linked)] <- 0
med$db_linked[is.na(med$db_linked)] <- 0

bp <- med[which(colnames(med) %in% c('SampleID', 'bp_linked', 'Drug_name'))]
sub <- bp[which(bp$bp_linked %in% 1),]
unique(sub$Drug_name)
table(sub$Drug_name)

db <- med[which(colnames(med) %in% c('SampleID', 'db_linked', 'Drug_name'))]
sub <- db[which(db$db_linked %in% 1),]
unique(sub$Drug_name)
table(sub$Drug_name)

########################################################################################

#### Run sensitivities for full models

# IHD
library(survival)
library(survminer)
library(tidyverse)

# Read in phenotypes
# Read in phenotypes
d1 <- readRDS("d1_202110_diseases_and_cancer.rds")

# Add PA as covariate 
PA <- read.csv("Censor_test/PA.csv")
# PA <- PA[which(colnames(PA) %in% c('SampleID', 'PA'))]
d1 <- left_join(d1, PA, by = 'SampleID')

# Add bp med as covariate 
med <- read.csv("Censor_test/binary_bp_variable.csv")
# d1 <- d1[which(d1$SampleID %in% med$SampleID),]
d1 <- left_join(d1, med, by = 'SampleID')
d1$BP_med <- ifelse(d1$BP_med == 1, 1, 0)
d1$BP_med[is.na(d1$BP_med)] <- 0

# Subset to those with med linkage available
covs <- read.csv("/Censor_test/binary_med_link.csv")
d1 <- d1[which(d1$SampleID %in% covs$SampleID),]
table(d1$BP_med)
clock <- names(d1)[56:1523] 
name <- 'IHD_FO'

  # Read linkage codes
  location_codes <- 'prepped_traits_used_with_dates/'
  codes <- read.csv(paste0(location_codes, name, '.csv'))

  # Save location
  location_prev <- '01_individual/Run_050723/Prevalent/'
  
  # Basic models per protein
  d1_data <- d1
  
  mat_hazard <- matrix(nrow=length(clock),ncol=10)
  output_hazard <- as.data.frame(mat_hazard)
  # for(j in 1){ 
  #   tryCatch({ 
  dat1= d1_data
  tmp1 = codes[which(codes$SampleID %in% dat1$SampleID),] 
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$SampleID %in% tmp1$SampleID),] 
  age_onset = codes[,c("first", "SampleID")]
  affected = merge(age_onset, affected, by= "SampleID")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1,4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$MOB))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$YOB)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$SampleID %in% codes$SampleID),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$SampleID.y <- NULL
  names(affected)[names(affected)=="SampleID"] <- "SampleID"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age_assessment
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < 0, "NA", cox$tte)
  # cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  
  # Run for IHD
  mat_hazard <- matrix(nrow=length(clock),ncol=10)
  output_hazard <- as.data.frame(mat_hazard)
  for(j in 1:length(clock)){ 
    tryCatch({ 
            print('Trait in neither list - running as standard.')  
            mod = coxph(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Sex) + cox$Age_assessment + cox$BMI + cox$Dep + as.factor(cox$Alc) + as.factor(cox$Smo) + as.factor(cox$Edu) + as.factor(cox$PA), data = cox)
            print(j)
            output_hazard[j,1] <- as.character(clock[[j]])
            output_hazard[j,2] <- as.character(name)
            output_hazard[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
            output_hazard[j,6] <- summary(mod)$coefficients[1,5]
            output_hazard[j,8] <- mod$n[1] - mod$nevent[1]
            output_hazard[j,7] <- mod$nevent[1]
            table <- cox.zph(mod)
            p1 <- table$table[,"p"]
            output_hazard[j,9] <-p1[1]
            output_hazard[j,10] <-p1[10]
            # cox$Event <- ifelse(cox$tte < -1, "NA", cox$Event)
            # write.csv(cox, paste0(location_tables, name ,'.csv'), row.names = F)
    }, error = function(e) cat("skipped"))
  } 
  
  # Save results
  comb <- output_hazard
  names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph_protein", "cox.zph_global")
  comb <- na.omit(comb)
  comb <- comb[order(comb$P.Value),]
  write.csv(comb, paste0('/Cox/med_comparison/', name, '_FULL_BP_NON_MED_SUBSET_V2.csv'), row.names = F)
 
  
  ########################################################################################
  
  #### Run sensitivities for full models
  
  # Run for IHD
  mat_hazard <- matrix(nrow=length(clock),ncol=10)
  output_hazard <- as.data.frame(mat_hazard)
  for(j in 1:length(clock)){ 
    tryCatch({ 
      print('Trait in neither list - running as standard.')  
      mod = coxph(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Sex) + cox$Age_assessment + cox$BMI + cox$Dep + as.factor(cox$Alc) + as.factor(cox$Smo) + as.factor(cox$Edu) + as.factor(cox$PA) + as.factor(cox$BP_med), data = cox)
      print(j)
      output_hazard[j,1] <- as.character(clock[[j]])
      output_hazard[j,2] <- as.character(name)
      output_hazard[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
      output_hazard[j,6] <- summary(mod)$coefficients[1,5]
      output_hazard[j,8] <- mod$n[1] - mod$nevent[1]
      output_hazard[j,7] <- mod$nevent[1]
      table <- cox.zph(mod)
      p1 <- table$table[,"p"]
      output_hazard[j,9] <-p1[1]
      output_hazard[j,10] <-p1[11]
      # cox$Event <- ifelse(cox$tte < -1, "NA", cox$Event)
      # write.csv(cox, paste0(location_tables, name ,'.csv'), row.names = F)
    }, error = function(e) cat("skipped"))
  } 
  
  # Save results
  comb <- output_hazard
  names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph_protein", "cox.zph_global")
  comb <- na.omit(comb)
  comb <- comb[order(comb$P.Value),]
  write.csv(comb, paste0('Cox/med_comparison/', name, 'FULL_BP_V2.csv'), row.names = F)
  
##########################################################################################
  
# Create a table of incident morbidity status for individuals

diseases <- c('AL_FO', 'ALS_FO', 'COPD_FO', 'CYS_FO', 'DEP_FO', 'Diab_FO', 'ENDO_FO', 
                'IBD_FO', 'IHD_FO', 'LIV_FO', 'LUP_FO', 'MS_FO', 'PD_FO', 'RA_FO', 'SCZ_FO', 'ST_FO', 'VD_FO',
              'DEATH', 'Prostate_FO', 'Breast_FO', 'BRAIN_FO', 'GYN_FO', 'LUNG_FO', 'Colorectal_FO')
  
files <- list.files('01_individual/Run_050723/Tables_cut/')

# IDs
d1 <- readRDS("Cox_preps/d1_cancer_reg_202012.rds")
ID <- d1[1]

library(tidyverse)

for(i in 1:length(diseases)){
  tryCatch({
    print(i)
    name <- as.character(diseases[i])
    print(name)
    # name <- sub('.csv', '', name)
    file <- read.csv(paste0('Censor_test/tables_cut/', diseases[i], '.csv'))
    file <- file[which(file$Event %in% 1),]
    file <- file[which(file$tte >= 0),]
    file <- file[c(1,16)]
    file$name <- 1
    file <- file[,-2]
    names(file)[2] <- name
    file[,2] <- as.character(file[,2])
    
    ID <- left_join(ID, file, by = 'SampleID')
  }, error = function(e) cat("skipped"))
}

ID[is.na(ID)] <- 0
write.csv(ID, 'Censor_test/Morbidity/covs_incident.csv', row.names = F)

# Tabulate with inclusion of sex specific contexts 
res <- data.frame(Individual = 1:length(rownames(ID)), Count = 1:length(rownames(ID)))

for(i in 1:length(rownames(ID))){
  individual <- ID[i,]
  individual <- individual[,-1]
  library(dplyr)
  individual <- individual %>% mutate_if(is.character, as.numeric)
  individual <- as.data.frame(t(individual))
  names(individual) <- 'X'
  count <- length(which(individual$X == 1))
  
  res[i,1] <- ID[i,1]
  res[i,2]<- count
  print(i)
}

res2 <- res[which(res$Individual %in% d1$SampleID),]
res <- res[order(-res$Count),]

length(which(res$Count >= 2))
length(which(res$Count >= 3))
length(which(res$Count >= 4))
length(which(res$Count >= 5))
length(which(res$Count >= 6))
length(which(res$Count >= 7))
write.csv(res, 'covs_incident_counts_table.csv', row.names = F)

####################################################################################

# Run models with m.status vs protein

res <- read.csv("Censor_test/Morbidity/covs_incident_counts_table.csv")
m_3 <- res[which(res$Count >= 3),]
olink <- read.csv('knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')
dat <- read.csv("parquet_47600_040723.csv")
d1 <- dat[,which(colnames(dat) %in% c('f.eid','f.21022.0.0', 'f.52.0.0', 'f.34.0.0', 'f.31.0.0', 'f.40000.0.0', 'f.40007.0.0'))]
names(d1) <- c('SampleID', 'Sex', 'YOB', 'MOB', 'Age_recruitment', 'DOD', 'Age_death')
dat = NULL
d1 <- d1[which(d1$SampleID %in% olink$SampleID),] # 47600
names(m_3)[1] <- 'SampleID'
d1 <- left_join(d1, m_3, by = 'SampleID')
d1$morb <- ifelse(d1$Count > 2, 1, 0)
d1$morb[is.na(d1$morb)] <- 0
d1 <- left_join(d1, olink, by = 'SampleID')
clock <- colnames(d1[11:1478])
library(glm2)
library(lme4)
d1$morb <- as.factor(d1$morb)
results <- data.frame(Prot = 1:length(clock), Effect = 1:length(clock), OR = 1:length(clock), SE = 1:length(clock), p = 1:length(clock))

for(i in 1:length(clock)){
  name <- as.character(clock[i])
  
  m <- glm(d1$morb ~ scale(d1[,name]) + Age_recruitment + factor(Sex), data = d1, family = binomial(link='logit'))
  
  Eff <- coef(summary(m))[2,1]
  SE <- coef(summary(m))[2,2]
  p <- coef(summary(m))[2,4]
  OR <- exp(Eff)
 
  # Save results
  results[i,1] <- name
  results[i,2] <- Eff
  results[i,3] <- OR
  results[i,4] <- SE
  results[i,5] <- p
  
  print(name)
}

res <- results[order(results$p),]
write.csv(res, 'inc_morbidity_results_raw.csv', row.names = F)

# 0.05/1468 = 0.000034

res <- read.csv('Morbidity/inc_morbidity_results_raw.csv')

# top <- res[which(res$p < 0.000034),]
top <- res[which(res$p < 0.0000031),]
# 765 associations

# Were the 54 proteins associated with >8 morbidities in incident models present 

prots <- c("GDF15", "IL6_1", "IL6_2", "PLAUR", 'NEFL', 'ASGR1', 'CHI3L1', 'IL6_3', 'IL6_4',
           "TNFRSF1A", 'CSF1', 'TNFRSF1B', "LGALS9", "CD74", "HAVCR2","CD300E" , "TNFRSF4",
           "CD274", "CD27", "TNF_2", 'TNF_1', 'TNF_3', 'TNF_4', 'CCL7', 'ST6GAL1', 
           'WFDC2', 'PRSS8', "IGFBP4", "BST2", 'HGF', 
           "TNFRSF10A", "TIMP1",     "VSIG4",     "MMP12" ,    "MSR1" ,     "CDCP1",
           "PGF" ,      "IL2RA" ,    "IL18BP",    "LAIR1"  ,   "LAMP3" ,    "CST3",
           "CCL3" ,     "CXCL9" ,    "TNFRSF9" ,  "LILRB4"  ,  "CXCL13" ,   "MDK",
           "TNFSF13",   "MZB1" ,     "TNFRSF14",  "ZBTB17" ,   "ITGA11" ,   "ITGAV"
)

top$prot <- sub("\\..*", "", top$Prot)
length(which(top$prot %in% prots))
prot_not <- prots[-which(prots %in% top$prot)]
res$Thr <- ifelse(res$p < 0.0000031, 1, 0)
res$prot <- sub("\\..*", "", res$Prot)
res$Morbid_prot <- ifelse(res$prot %in% prots, 1, 0)

# Examine largest effect size
eff <- res[order(-res$OR),]
write.csv(res, 'inc_morbidity_results_formatted.csv', row.names = F)

################################################################################v

# Process medications results

non <- read.csv("IHD_FO_FULL_BP_NON_MED_SUBSET_V2.csv")
non <- non[which(non$P.Value < 0.0000031),]

# Check that these were present previously 
assoc <- read.csv("Cox/Associations_retained_Bon.csv")
assoc <- assoc[which(assoc$Outcome %in% 'IHD_FO'),]
assoc$assoc <- paste0(assoc$Outcome, assoc$Predictor)
non$assoc <- paste0(non$Outcome, non$Predictor)
non2 <- non[which(non$assoc %in% assoc$assoc),]

# Read in with bp added
med <- read.csv("med_comparison/IHD_FOFULL_BP_V2.csv")
med <- med[which(med$P.Value < 0.0000031),]
med$assoc <- paste0(med$Outcome, med$Predictor)
med2 <- med[which(med$assoc %in% assoc$assoc),]

# Are all the same assocs present
length(which(med2$assoc %in% non2$assoc))
length(which(non2$assoc %in% med2$assoc))
length(which(!non2$assoc %in% med2$assoc))
diff <- non2[which(!non2$assoc %in% med2$assoc),]
write.csv(diff, 'binary_med_36_attenuated.csv', row.names = F)

# Plot the concordance between effect sizes in top 100
# Top 100 from main set
med <- read.csv("Cox/med_comparison/IHD_FOFULL_BP_V2.csv")
med$assoc <- paste0(med$Outcome, med$Predictor)
non3 <- left_join(non2, med, by = 'assoc')
non3$attenuated <- ifelse(non3$assoc %in% diff$assoc, 'yes', 'no')
library(ggplot2)

pdf("IHD_rep_V2.pdf", width = 5, height = 5)
ggplot(non3, aes(x = jitter(Hazard.Ratio.x), Hazard.Ratio.y, color = attenuated)) +
  geom_vline(xintercept = 1, linetype = 1) +
  geom_hline(yintercept = 1, linetype = 1) +
  geom_point(size = 0.8, alpha = 0.4) +
  theme_minimal() +
  scale_color_manual(values = c('steelblue4', 'red3')) + 
  xlab("HR") +
  ylab("HR medication") +
  geom_abline(slope=1, intercept=0, linetype = 1)  + 
  theme(legend.position = 'none')
dev.off()

write.csv(non3, 'med_comparison/non3_all.csv', row.names = F)
cor.test(abs(log(non3$Hazard.Ratio.x)), abs(log(non3$Hazard.Ratio.y)))
