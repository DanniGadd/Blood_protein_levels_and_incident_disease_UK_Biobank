###############################################################################################

## Health linkage preps for cox

###############################################################################################

# srun -p interactive --pty bash

# module load R

R

# Location prepped lists used
out_1 <- '/path_to_file.../prepped_lists_used/'
# Location prepped traits diseases
out_2 <- '/path_to_file.../prepped_traits_used_with_dates/'
# Location prepped traits cancer
out_3 <- '/path_to_file.../prepped_traits_cancer_reg/'

###############################################################################################

### First occurrence traits - preps

###############################################################################################

### LUP first occurrence

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'M32'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'LUP_FO.csv'), row.names = F)

# How many cases occurred after 2007
sub <- dis[which(dis$first > 200701),]


### SCZ first occurrence

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'F20'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'SCZ_FO.csv'), row.names = F)



### RA first occurrence

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis1 <- dis[which(dis$icd10_2digit %in% 'M06'),]
dis2 <- dis[which(dis$icd10_2digit %in% 'M05'),]
dis <- rbind(dis1, dis2)
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_2,'RA_FO.csv'), row.names = F)


### CYS first occurrence

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'N30'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'CYS_FO.csv'), row.names = F)

### MS first occurrence

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'G35'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'MS_FO.csv'), row.names = F)


### ENDO first occurrence

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'N80'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'ENDO_FO.csv'), row.names = F)


### VD first occurrence

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'F01'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'VD_FO.csv'), row.names = F)


### MND

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'G12'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'ALS_FO.csv'), row.names = F)


### DEP

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'F33'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'DEP_FO.csv'), row.names = F)


### AL

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis1 <- dis[which(dis$icd10_2digit %in% 'F00'),]
dis2 <- dis[which(dis$icd10_2digit %in% 'G30'),]
dis <- rbind(dis1, dis2)
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_2,'AL_FO.csv'), row.names = F)

### PD

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'G20'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'PD_FO.csv'), row.names = F)


### LIV

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis1 <- dis[which(dis$icd10_2digit %in% 'K70'),]
dis2 <- dis[which(dis$icd10_2digit %in% 'K72'),]
dis3 <- dis[which(dis$icd10_2digit %in% 'K74'),]
dis4 <- dis[which(dis$icd10_2digit %in% 'K71'),]
dis5 <- dis[which(dis$icd10_2digit %in% 'K73'),]
dis <- rbind(dis1, dis2)
dis <- rbind(dis, dis3)
dis <- rbind(dis, dis4)
dis <- rbind(dis, dis5)
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_2,'LIV_FO.csv'), row.names = F)


### IBD

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis1 <- dis[which(dis$icd10_2digit %in% 'K51'),]
dis2 <- dis[which(dis$icd10_2digit %in% 'K50'),]
dis <- rbind(dis1, dis2)
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_2,'IBD_FO.csv'), row.names = F)


### ST

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'I63'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'ST_FO.csv'), row.names = F)


### CKD

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis <- dis[which(dis$icd10_2digit %in% 'N18'),]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_2,'CKD_FO.csv'), row.names = F)

### COPD 

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

# dis1 <- dis[which(dis$icd10_2digit %in% 'J40'),]
dis2 <- dis[which(dis$icd10_2digit %in% 'J41'),]
dis3 <- dis[which(dis$icd10_2digit %in% 'J42'),]
dis4 <- dis[which(dis$icd10_2digit %in% 'J43'),]
dis5 <- dis[which(dis$icd10_2digit %in% 'J44'),]
# dis <- rbind(dis1, dis2)
dis <- rbind(dis2, dis3)
dis <- rbind(dis, dis4)
dis <- rbind(dis, dis5)
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_2,'COPD_FO.csv'), row.names = F)


### CHD 

dis <- readRDS("/path_to_file.../ukbb_icd10first_diag_date.rds")

dis1 <- dis[which(dis$icd10_2digit %in% 'I20'),]
dis2 <- dis[which(dis$icd10_2digit %in% 'I21'),]
dis3 <- dis[which(dis$icd10_2digit %in% 'I22'),]
dis4 <- dis[which(dis$icd10_2digit %in% 'I23'),]
dis5 <- dis[which(dis$icd10_2digit %in% 'I24'),]
dis6 <- dis[which(dis$icd10_2digit %in% 'I25'),]
dis <- rbind(dis1, dis2)
dis <- rbind(dis, dis3)
dis <- rbind(dis, dis4)
dis <- rbind(dis, dis5)
dis <- rbind(dis, dis6)
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_2,'IHD_FO.csv'), row.names = F)


### Cancer registry data  

### BC  

dis <- readRDS("/path_to_file.../ukbb_cancer_registry_2digit_earliestdate.rds")

dis <- dis[which(dis$icd10_2digit %in% 'C50'),]
dis <- dis[c(1,3,2)]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_3,'Breast_FO.csv'), row.names = F)


### Prostate 

dis <- readRDS("/path_to_file.../ukbb_cancer_registry_2digit_earliestdate.rds")

dis <- dis[which(dis$icd10_2digit %in% 'C61'),]
dis <- dis[c(1,3,2)]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_3,'Prostate_FO.csv'), row.names = F)

# Lung cancer FO

dis <- readRDS("/path_to_file.../ukbb_cancer_registry_2digit_earliestdate.rds")

dis <- dis[which(dis$icd10_2digit %in% 'C34'),]
dis <- dis[c(1,3,2)]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
write.csv(dis, paste0(out_3,'LUNG_FO.csv'), row.names = F)


# Colorectal

dis <- readRDS("/path_to_file.../ukbb_cancer_registry_2digit_earliestdate.rds")

dis1 <- dis[which(dis$icd10_2digit %in% 'C18'),]
dis2 <- dis[which(dis$icd10_2digit %in% 'C20'),]
dis3 <- dis[which(dis$icd10_2digit %in% 'C21'),]

dis <- rbind(dis1, dis2)
dis <- rbind(dis, dis3)

dis <- dis[c(1,3,2)]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_3,'Colorectal_FO.csv'), row.names = F)


# Gyn

dis <- readRDS("/path_to_file.../ukbb_cancer_registry_2digit_earliestdate.rds")

list <- c('C51', 'C52', 'C53', 'C54', 'C55', 'C56', 'C57', 'C58')

dis <- dis[which(dis$icd10_2digit %in% list),]

# C51  C52  C53  C54  C55  C56  C57  C58
# 205   39  570 2389   59 1615  156   13

dis <- dis[c(1,3,2)]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_3,'GYN_FO.csv'), row.names = F)


# CNS / brain

dis <- readRDS("/path_to_file.../ukbb_cancer_registry_2digit_earliestdate.rds")

list <- c('C71', 'C72')

# C71 C72
# 894  29

dis <- dis[which(dis$icd10_2digit %in% list),]

dis <- dis[c(1,3,2)]
names(dis) <- c('SampleID', 'first', 'code')
dis <- as.data.frame(dis)
dis$first <- sub('-', '', dis$first)
dis$first <- sub('-', '', dis$first)
dis$first <- substr(dis$first,1,6)
dis <- dis[which(nchar(dis$first) == 6),]
dis <- dis[order(dis$first),]
# dis <- dis[-which(duplicated(dis$SampleID)),]
write.csv(dis, paste0(out_3,'BRAIN_FO.csv'), row.names = F)


###############################################################################################

##### MORTALITY DATA 

###############################################################################################

dat <- read.csv("/path_to_file.../parquet_54189_260722.csv")
t <- read.csv("/path_to_file.../eric_individual_fieldIDs_extracted.csv")

# Read in olink file prepped for cox models 
proteins <- read.csv('/path_to_file.../knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')


# Mortality
dat <- read.csv("/path_to_file.../parquet_54189_260722.csv")
d1 <- dat[,which(colnames(dat) %in% c('f.eid','f.21022.0.0', 'f.52.0.0', 'f.34.0.0', 'f.31.0.0', 'f.40000.0.0', 'f.40007.0.0'))]
names(d1) <- c('SampleID', 'Sex', 'YOB', 'MOB', 'Age_recruitment', 'DOD', 'Age_death')

d1 <- d1[which(d1$SampleID %in% proteins$SampleID),]

death <- d1[c(1,6)]
death <- na.omit(death)
death$DOD <- gsub('-', '', death$DOD)
names(death) <- c('SampleID', 'first')
death <- death[order(death$first),]
death$first <- substr(death$first,1,6) 
DEATH <- death # 4445
write.csv(DEATH, paste0('/path_to_file.../prepped_death/', 'DEATH.csv'), row.names = F)

# ########################################################################################

### ALS - check against manual extraction of codes from GP/ICD sources

# library(readxl)
clinical <- read.csv("/path_to_file.../clinical_r2_conv.csv")

library(readxl)
phen <- read_excel(paste0(phen_location, "ALS_phenotype_PH62_ver_124_concepts_20220719T110535_Rob.xlsx"))
phen <- as.data.frame(phen)
phen <- phen[which(phen$Rob_include %in% 'Y'),]
phen$code <- gsub('.{2}$', '', phen$code)
patterns <- phen$code
sub <- clinical[which(clinical$read_2 %in% 'F1520'),]
data.frame(table(sub$read_2))
sub <- sub[c(1,4,3)]
names(sub) <- c('SampleID', 'code', 'first')
sub$first <- sub('/', '', sub$first)
sub$first <- sub('/', '', sub$first)
sub$y <- substr(sub$first,5,8)
sub$m <- substr(sub$first,3,4)
sub$first <- paste0(sub$y, sub$m)
sub$y <- NULL
sub$m <- NULL
save <- phen[which(phen$code %in% clinical$read_2),]
write.csv(phen, paste0(out_1,'ALS.csv'), row.names = F)

ICD9 <- read.csv(paste0(ICD_location, "ICD9_formatted.csv"))
ICD10 <- read.csv(paste0(ICD_location,"ICD10_formatted.csv"))
ICD10 <- ICD10[grep('G122', ICD10$code),]
ICD9 <- ICD9[grep('3352', ICD9$code),]
ICD9$first <- gsub('-', '', ICD9$first)
ICD10$first <- gsub('-', '', ICD10$first)
ICD9$first <- substr(ICD9$first,1,6)
ICD10$first <- substr(ICD10$first,1,6)

# check to ensure dates have correct number of characters in strings
which(nchar(ICD10$first) > 6)
which(nchar(ICD9$first) > 6)

# Check ordering of dates to spot any odd formats
sub <- sub[order(sub$first),]

sub$type <- 'GP'
ICD9$type <- 'ICD9'
ICD10$type <- 'ICD10'

dat <- rbind(sub, ICD10)
dat <- rbind(dat, ICD9)

subset <- dat[order(dat$first),]
subset <- subset[-which(nchar(subset$first) < 6),]
subset$year <- substr(subset$first, 1,4)
subset <-subset[-which(duplicated(subset$SampleID)),]
ALS <- subset

# write.csv(ALS, paste0(out_2, 'ALS.csv'), row.names = F)
write.csv(ALS, '/path_to_file.../ALS_int_MND.csv', row.names = F)
write.csv(ALS, '/path_to_file.../ALS_int_ALS_only.csv', row.names = F)

# Now look at which of the individuals from FO are ALS
ALS <- read.csv("/path_to_file.../ALS.csv")
ALS_FO <- read.csv("/path_to_file.../ALS_FO.csv")

case1 <- ALS[which(ALS$Event %in% 1),]
case2 <- ALS_FO[which(ALS_FO$Event %in% 1),]
which(case2$SampleID %in% case1$SampleID)

# ALS only codes
ALS_codes <- read.csv('/path_to_file.../ALS_int_ALS_only.csv')
case2 <- left_join(case2, ALS_codes, by = "SampleID")
# G12.2 and F1520 read3 codes - checked against (exclude non specific)

# ########################################################################################
