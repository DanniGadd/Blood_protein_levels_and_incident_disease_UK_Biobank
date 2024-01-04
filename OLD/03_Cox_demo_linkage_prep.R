###############################################################################################

## Health linkage preps for cox

###############################################################################################

# srun -p interactive --pty bash

# module load R

R

### PREP BASIS FOR COX MODELS

library(survival)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

# Process phenotype (d1) and death / alive date data 
dat <- read.csv("/path_to_file.../parquet_47600_040723.csv")
d1 <- dat[,which(colnames(dat) %in% c('f.eid','f.21022.0.0', 'f.52.0.0', 'f.34.0.0', 'f.31.0.0', 'f.40000.0.0', 'f.40007.0.0'))]
names(d1) <- c('SampleID', 'Sex', 'YOB', 'MOB', 'Age_recruitment', 'DOD', 'Age_death')
dat = NULL

# Add date of assessment cenre from Eric's extraction
t <- read.csv("/home/dgadd/PPP_core_input_files/eric_individual_fieldIDs_extracted.csv")
date <- t[which(colnames(t) %in% c('f.eid', 'f.53.0.0'))]
names(date) <- c('SampleID', 'DOA')
t = NULL

d1 <- left_join(d1, date, by = 'SampleID')

# Subset to complete protein data individuals - prepped data
olink_internal <- read.csv('/path_to_file.../knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')
d1 <- d1[which(d1$SampleID %in% olink_internal$SampleID),] # 47600

# Format DOA date info 
d1$DOA <- gsub('-', '', d1$DOA)
d1$DOA <- substr(d1$DOA,1,6)

# Format DOD date info 
d1$DOD <- gsub('-', '', d1$DOD)
d1$DOD <- substr(d1$DOD,1,6)

# Format so that mob has '01' rather than '1' for months less than 10, calculating DOB
subset1 <- d1[which(d1$MOB < 10),]
subset2 <- d1[which(d1$MOB > 9),]
subset1$MOB <- sub("^", "0", subset1$MOB) # adds 0 in front of any single digit month 
d1 <- rbind(subset1, subset2)
d1$DOB <- paste0(d1$YOB, d1$MOB)

# Work out those who were reported as having died during UKB follow up 
age_dead <- d1[complete.cases(d1$Age_death),] 
age_dead$dead <- 1
age_alive <- d1[-which(d1$SampleID %in% age_dead$SampleID),]
age_alive$dead <- 0

## extract year/month of death as separate variables ##
age_dead$y_dead <- as.numeric(substr(age_dead$DOD, 1, 4))
age_dead$m_dead <- as.numeric(substr(age_dead$DOD, 5, 6))

# Filter death data to earliest sampling date (i.e. check for any death exclusions)
age_dead_include <- age_dead[which(age_dead$DOD <= 202206),] # died before sampling date
age_dead_exclude <- age_dead[which(age_dead$DOD > 202206),] # died after sampling date and have been excluded 

# Calculate a more exact estimate (by year and month) for age of death in the included individuals that died 
age_dead_include$y_diff <- age_dead_include$y_dead - as.numeric(age_dead_include$YOB)
age_dead_include$m_diff <- (age_dead_include$m_dead - as.numeric(age_dead_include$MOB))/12
age_dead_include$diff <- age_dead_include$y_diff + age_dead_include$m_diff

# Work out those who are alive (i.e. everytone not in the list of included dead people from above)
age_alive <- d1[-which(d1$SampleID %in% age_dead_include$SampleID),] # people alive at censor date

# Ensure that all individuals in the samples are coded as alive and all individuals in the dead file are coded as such 
age_alive$dead <- 0
age_dead_include$dead <- 1

# Find age the 'alive' people were at censor
age_alive$y_diff <- 2022 - as.numeric(age_alive$YOB)
age_alive$m_diff <- (06 - as.numeric(age_alive$MOB))/12
age_alive$diff <- age_alive$y_diff + age_alive$m_diff

# So for those that died, we are taking forward their age at death 
# And for those that were alive, we are taking forward their age at censor date

# Subset to just the cols needed for joining dead and alive data
age_alive <- age_alive[c('SampleID', 'Sex', 'Age_recruitment', 'DOD', 'Age_death', 'DOA', 'YOB', 'MOB', 'DOB', 'dead', 'diff')] 
age_dead_include <- age_dead_include[c('SampleID', 'Sex', 'Age_recruitment', 'DOD', 'Age_death', 'DOA', 'YOB', 'MOB', 'DOB', 'dead', 'diff')] 
names(age_alive)[11] <- 'aged'
names(age_dead_include)[11] <- 'aged'

# Bind dead and alive individuals together
d1 <- rbind(age_alive, age_dead_include)
dim(d1) # 47600

# Get a more precise (year and month) estimate of age at recruitment (using the date of assessment) and DOB
d1$DOA_year <- as.numeric(substr(d1$DOA, 1, 4))
d1$DOA_month <- as.numeric(substr(d1$DOA, 5, 6))

d1$Age_y <- d1$DOA_year - as.numeric(d1$YOB)
d1$Age_m <- (d1$DOA_month - as.numeric(d1$MOB))/12
d1$Age_assessment <- d1$Age_y + d1$Age_m

# Remove old age recrutiment variable 
d1 <- d1[-c(3)]

# Merge protein data into phenotypes file  
d1 <- left_join(d1, olink_internal, by = 'SampleID')

### Add covariates
covs <- read.csv("/path_to_file.../covariates_joint.csv")
d2 <- left_join(covs, d1, by = "SampleID")
saveRDS(d2, '/path_to_file.../d1_202206_test.rds')

### Now feed into the cox models as the phenotype file (non-cancer diseases)               

                                                          
###############################################################################################

## PREP FOR MORTALITY MODELS - CENSOR DATE 

# srun -p interactive --pty bash

### PREP BASIS FOR COX MODELS

library(survival)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

# Process phenotype (d1) and death / alive date data 
dat <- read.csv("/path_to_file.../parquet_47600_040723.csv")
d1 <- dat[,which(colnames(dat) %in% c('f.eid','f.21022.0.0', 'f.52.0.0', 'f.34.0.0', 'f.31.0.0', 'f.40000.0.0', 'f.40007.0.0'))]
names(d1) <- c('SampleID', 'Sex', 'YOB', 'MOB', 'Age_recruitment', 'DOD', 'Age_death')
dat = NULL

# Add date of assessment cenre from Eric's extraction
t <- read.csv("/path_to_file.../eric_individual_fieldIDs_extracted.csv")
date <- t[which(colnames(t) %in% c('f.eid', 'f.53.0.0'))]
names(date) <- c('SampleID', 'DOA')
t = NULL

d1 <- left_join(d1, date, by = 'SampleID')

# Subset to complete protein data individuals - prepped data
olink_internal <- read.csv('/path_to_file.../knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')
d1 <- d1[which(d1$SampleID %in% olink_internal$SampleID),] 

# Format DOA date info 
d1$DOA <- gsub('-', '', d1$DOA)
d1$DOA <- substr(d1$DOA,1,6)

# Format DOD date info 
d1$DOD <- gsub('-', '', d1$DOD)
d1$DOD <- substr(d1$DOD,1,6)

# Format so that mob has '01' rather than '1' for months less than 10, calculating DOB
subset1 <- d1[which(d1$MOB < 10),]
subset2 <- d1[which(d1$MOB > 9),]
subset1$MOB <- sub("^", "0", subset1$MOB) # adds 0 in front of any single digit month 
d1 <- rbind(subset1, subset2)
d1$DOB <- paste0(d1$YOB, d1$MOB)

# Work out those who were reported as having died during UKB follow up 
age_dead <- d1[complete.cases(d1$Age_death),] 
age_dead$dead <- 1
age_alive <- d1[-which(d1$SampleID %in% age_dead$SampleID),] 
age_alive$dead <- 0

## extract year/month of death as separate variables ##
age_dead$y_dead <- as.numeric(substr(age_dead$DOD, 1, 4))
age_dead$m_dead <- as.numeric(substr(age_dead$DOD, 5, 6))

# Filter death data to earliest sampling date (i.e. check for any death exclusions)
age_dead_include <- age_dead[which(age_dead$DOD <= 202211),]
age_dead_exclude <- age_dead[which(age_dead$DOD > 202211),] 

# Calculate a more exact estimate (by year and month) for age of death in the included individuals that died 
age_dead_include$y_diff <- age_dead_include$y_dead - as.numeric(age_dead_include$YOB)
age_dead_include$m_diff <- (age_dead_include$m_dead - as.numeric(age_dead_include$MOB))/12
age_dead_include$diff <- age_dead_include$y_diff + age_dead_include$m_diff

# Work out those who are alive (i.e. everytone not in the list of included dead people from above)
age_alive <- d1[-which(d1$SampleID %in% age_dead_include$SampleID),] 

# Ensure that all individuals in the samples are coded as alive and all individuals in the dead file are coded as such 
age_alive$dead <- 0
age_dead_include$dead <- 1

# Find age the 'alive' people were at censor
age_alive$y_diff <- 2022 - as.numeric(age_alive$YOB)
age_alive$m_diff <- (11 - as.numeric(age_alive$MOB))/12
age_alive$diff <- age_alive$y_diff + age_alive$m_diff

# So for those that died, we are taking forward their age at death 
# And for those that were alive, we are taking forward their age at censor date

# Subset to just the cols needed for joining dead and alive data
age_alive <- age_alive[c('SampleID', 'Sex', 'Age_recruitment', 'DOD', 'Age_death', 'DOA', 'YOB', 'MOB', 'DOB', 'dead', 'diff')] 
age_dead_include <- age_dead_include[c('SampleID', 'Sex', 'Age_recruitment', 'DOD', 'Age_death', 'DOA', 'YOB', 'MOB', 'DOB', 'dead', 'diff')] 
names(age_alive)[11] <- 'aged'
names(age_dead_include)[11] <- 'aged'

# Bind dead and alive individuals together
d1 <- rbind(age_alive, age_dead_include)

# Get a more precise (year and month) estimate of age at recruitment (using the date of assessment) and DOB
d1$DOA_year <- as.numeric(substr(d1$DOA, 1, 4))
d1$DOA_month <- as.numeric(substr(d1$DOA, 5, 6))

d1$Age_y <- d1$DOA_year - as.numeric(d1$YOB)
d1$Age_m <- (d1$DOA_month - as.numeric(d1$MOB))/12
d1$Age_assessment <- d1$Age_y + d1$Age_m

d1 <- d1[-c(3)]

# Merge protein data into phenotypes file  
d1 <- left_join(d1, olink_internal, by = 'SampleID')

### Add covariates
covs <- read.csv("/path_to_file.../covariates_joint.csv")
d1 <- left_join(covs, d1, by = "SampleID")
saveRDS(d1, '/path_to_file.../d1_mortality_202112.rds')

### Now feed into the cox models as the phenotype file for mortality

###############################################################################################

## PREP FOR CANCER REGISTRY MODELS - CENSOR DATE 

# srun -p interactive --pty bash

### PREP BASIS FOR COX MODELS

library(survival)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

# Process phenotype (d1) and death / alive date data 
dat <- read.csv("/path_to_file.../parquet_47600_040723.csv")
d1 <- dat[,which(colnames(dat) %in% c('f.eid','f.21022.0.0', 'f.52.0.0', 'f.34.0.0', 'f.31.0.0', 'f.40000.0.0', 'f.40007.0.0'))]
names(d1) <- c('SampleID', 'Sex', 'YOB', 'MOB', 'Age_recruitment', 'DOD', 'Age_death')
dat = NULL

# Add date of assessment cenre from Eric's extraction
t <- read.csv("/path_to_file.../eric_individual_fieldIDs_extracted.csv")
date <- t[which(colnames(t) %in% c('f.eid', 'f.53.0.0'))]
names(date) <- c('SampleID', 'DOA')
t = NULL

d1 <- left_join(d1, date, by = 'SampleID')

# Subset to complete protein data individuals - prepped data
olink_internal <- read.csv('/path_to_file.../knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')
d1 <- d1[which(d1$SampleID %in% olink_internal$SampleID),] 

# Format DOA date info 
d1$DOA <- gsub('-', '', d1$DOA)
d1$DOA <- substr(d1$DOA,1,6)

# Format DOD date info 
d1$DOD <- gsub('-', '', d1$DOD)
d1$DOD <- substr(d1$DOD,1,6)

# Format so that mob has '01' rather than '1' for months less than 10, calculating DOB
subset1 <- d1[which(d1$MOB < 10),]
subset2 <- d1[which(d1$MOB > 9),]
subset1$MOB <- sub("^", "0", subset1$MOB) # adds 0 in front of any single digit month 
d1 <- rbind(subset1, subset2)
d1$DOB <- paste0(d1$YOB, d1$MOB)

# Work out those who were reported as having died during UKB follow up 
age_dead <- d1[complete.cases(d1$Age_death),] 
age_dead$dead <- 1
age_alive <- d1[-which(d1$SampleID %in% age_dead$SampleID),] 
age_alive$dead <- 0

## extract year/month of death as separate variables ##
age_dead$y_dead <- as.numeric(substr(age_dead$DOD, 1, 4))
age_dead$m_dead <- as.numeric(substr(age_dead$DOD, 5, 6))

# Filter death data to earliest sampling date (i.e. check for any death exclusions)
age_dead_include <- age_dead[which(age_dead$DOD <= 201612),] 
age_dead_exclude <- age_dead[which(age_dead$DOD > 201612),] 

# Calculate a more exact estimate (by year and month) for age of death in the included individuals that died 
age_dead_include$y_diff <- age_dead_include$y_dead - as.numeric(age_dead_include$YOB)
age_dead_include$m_diff <- (age_dead_include$m_dead - as.numeric(age_dead_include$MOB))/12
age_dead_include$diff <- age_dead_include$y_diff + age_dead_include$m_diff

# Work out those who are alive (i.e. everytone not in the list of included dead people from above)
age_alive <- d1[-which(d1$SampleID %in% age_dead_include$SampleID),] 

# Ensure that all individuals in the samples are coded as alive and all individuals in the dead file are coded as such 
age_alive$dead <- 0
age_dead_include$dead <- 1

# Find age the 'alive' people were at censor
age_alive$y_diff <- 2016 - as.numeric(age_alive$YOB)
age_alive$m_diff <- (12 - as.numeric(age_alive$MOB))/12
age_alive$diff <- age_alive$y_diff + age_alive$m_diff

# So for those that died, we are taking forward their age at death 
# And for those that were alive, we are taking forward their age at censor date

# Subset to just the cols needed for joining dead and alive data
age_alive <- age_alive[c('SampleID', 'Sex', 'Age_recruitment', 'DOD', 'Age_death', 'DOA', 'YOB', 'MOB', 'DOB', 'dead', 'diff')] 
age_dead_include <- age_dead_include[c('SampleID', 'Sex', 'Age_recruitment', 'DOD', 'Age_death', 'DOA', 'YOB', 'MOB', 'DOB', 'dead', 'diff')] 
names(age_alive)[11] <- 'aged'
names(age_dead_include)[11] <- 'aged'

# Bind dead and alive individuals together
d1 <- rbind(age_alive, age_dead_include)
dim(d1) # 47600  11

# Get a more precise (year and month) estimate of age at recruitment (using the date of assessment) and DOB
d1$DOA_year <- as.numeric(substr(d1$DOA, 1, 4))
d1$DOA_month <- as.numeric(substr(d1$DOA, 5, 6))

d1$Age_y <- d1$DOA_year - as.numeric(d1$YOB)
d1$Age_m <- (d1$DOA_month - as.numeric(d1$MOB))/12
d1$Age_assessment <- d1$Age_y + d1$Age_m

d1 <- d1[-c(3)]

# Merge protein data into phenotypes file  
d1 <- left_join(d1, olink_internal, by = 'SampleID')

### Add covariates
covs <- read.csv("/path_to_file.../covariates_joint.csv")
d1 <- left_join(covs, d1, by = "SampleID")
saveRDS(d1, '/path_to_file.../d1_cancer_reg_202012.rds')

### Now feed into the cox models as the phenotype file for cancer traits                    

                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
