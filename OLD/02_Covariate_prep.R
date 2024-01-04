###############################################################################################

## Covariate prep

###############################################################################################

# srun -p interactive --pty bash

# module load R 

R

library(ggplot2)
library(data.table)
library(tidyverse)

# Prepped data
olink_internal <- read.csv('/path_to_file.../knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')

# Extraction (v1)
dat <- read.csv("/path_to_file.../parquet_47600_040723.csv")
dat <- dat[which(dat$f.eid %in% olink_internal$SampleID),]

not <- olink_internal[-which(olink_internal$SampleID %in% dat$f.eid),]

covs <- dat[,which(colnames(dat) %in% c('f.eid','f.189.0.0', 'f.20116.0.0', 'f.1558.0.0', 'f.21001.0.0'))]
names(covs) <- c('SampleID', 'Dep', 'Alc', 'Smo', 'BMI')

# Load in additional covariates from ben GWAS
cov <- fread("/path_to_file.../combined_covars_v1.tsv")
cov <- cov[which(cov$IID %in% olink_internal$SampleID),]
table(is.na(cov$UKBPC_16)) # FALSE- 48821 with genetic data

# Remove prefer no to answer for alc and smo
length(which(covs$Smo < 0)) # 186
covs$Smo[covs$Smo < 0] <- NA

length(which(covs$Alc < 0)) # 57
covs$Alc[covs$Alc < 0] <- NA

# Join up covariates into one set
library(tidyverse)
names(cov)[2] <- 'SampleID'
joint <- left_join(covs, cov, by = 'SampleID')

# Add updated data linkage to education qualifications
t <- read.csv("/path_to_file.../eric_individual_fieldIDs_extracted.csv")
t <- t[c(1,10:33)]
ed <- t[which(t[,1] %in% olink_internal[,1]),]

# Recode
#1	College or University degree
#2	A levels/AS levels or equivalent
#3	O levels/GCSEs or equivalent
#4	CSEs or equivalent
#5	NVQ or HND or HNC or equivalent
#6	Other professional qualifications eg: nursing, teaching
#-7	None of the above
#-3	Prefer not to answer

# For each individual, search for whether they had a 1 (college or uni degree) and index their IDs

for(i in 2:25){
 print(colnames(ed)[i])
 name <- as.character(colnames(ed)[i])
 print(table(ed[,name]))
}

sub1 <- ed[,c("f.eid","f.6138.0.0")]
sub1 <- sub1[which(sub1[,2] %in% 1),]
names(sub1)[2] <- 'X'
sub2 <- ed[,c("f.eid","f.6138.1.0")]
sub2 <- sub2[which(sub2[,2] %in% 1),]
names(sub2)[2] <- 'X'
sub3 <- ed[,c("f.eid","f.6138.2.0")]
sub3 <- sub3[which(sub3[,2] %in% 1),]
names(sub3)[2] <- 'X'
sub4 <- ed[,c("f.eid","f.6138.3.0")]
sub4 <- sub4[which(sub4[,2] %in% 1),]
names(sub4)[2] <- 'X'
bind <- rbind(sub1, sub2)
bind <- rbind(bind, sub3)
bind <- rbind(bind, sub4)
names(bind)[1] <- 'SampleID'
bind <- bind[which(bind$SampleID %in% olink_internal$SampleID),]
length(unique(bind[,1])) 
names(ed) <- c('SampleID', 'Educ') # 
names(ed)[1] <- 'SampleID'
ed$Edu <- ifelse(ed$SampleID %in% bind$SampleID, 1, 0)
ed <- ed[c(1,26)]
joint <- left_join(joint, ed, by = 'SampleID')

write.csv(joint, '/path_to_file.../all_quant_20200904.pheno.csv', row.names = F)

###############################################################################################

## Additional covariate preps

###############################################################################################

### Load in 47600 individuals present in analyses
olink <- readRDS('/path_to_file.../knn_imputed_processed_olink_internal_proteins.rds')

# Get core info to join in covariates 
ID <- olink[c(1:2)]

###############################################################################################

## PGS

PGS <- read.table("/path_to_file.../PGS.tab", header = T)
names(PGS)[1] <- 'SampleID'

# Extract relevant PGS for this study
PGS <- PGS[which(colnames(PGS) %in% c('SampleID', 'f.26206.0.0',
                                      'f.26223.0.0','f.26248.0.0','f.26260.0.0','f.26273.0.0','f.26285.0.0'))]

names(PGS) <- c('SampleID', 'PGS_AL', 'PGS_CVD', 'PGS_IST', 'PGS_PD', 'PGS_RA','PGS_T2D')

library(tidyverse)
ID <- left_join(ID, PGS, by = 'SampleID')

###############################################################################################

## Extraction variables
ex <- read.delim("/path_to_file.../added_phenos_june_23.txt")
names(ex)[1] <- 'SampleID'

ex <- ex[which(colnames(ex) %in% c('SampleID', 'f.884.0.0',
                                  'f.20003.0.0'))]
names(ex) <- c('SampleID', 'PA', 'Med')

# Subset to overlap with PPP data
ID <- left_join(ID, ex, by = 'SampleID')
table(is.na(ID$PA))
table(is.na(ID$Med))
# 47544 PA
# 35073 Med

# Process PA
ID$PA <- ifelse(ID$PA == '-1', NA, ID$PA)
ID$PA <- ifelse(ID$PA == '-3', NA, ID$PA)
ID$PA <- ifelse(ID$PA == 0 | ID$PA == 1 | ID$PA == 2, 0, ifelse(ID$PA == 3 | ID$PA == 4, 1, 2))

###############################################################################################

## Phenotypes file (majority of covariates extracted)

# Load phenotypes file
phen <- read.table(gzfile("/path_to_file.../all_quant_20200904.pheno.v2.tsv.gz"), header = T)
names(phen)[1] <- 'SampleID'

phen <- phen[which(phen$SampleID %in% ID$SampleID),]

covs <- phen[which(colnames(phen) %in% c('SampleID', 'f_30750_0_0', 'f_30710_0_0',
                                         'f_30690_0_0', 'f_30760_0_0', 'f_30780_0_0',
                                         'f_48_0_0', 'f_49_0_0', 'f_4080_0_0', 'f_30870_0_0', 'f_30740_0_0',
                                         'f_30700_0_0', 'f_30720_0_0', 'f_30670_0_0', 'f_30880_0_0',
                                         'f_30650_0_0', 'f_30650_0_0', 'f_30620_0_0','f_30600_0_0',
                                         'f_30020_0_0', 'f_30050_0_0', 'f_30060_0_0',
                                         'f_30040_0_0', 'f_30080_0_0', 'f_30010_0_0', 'f_30000_0_0'))]


names(covs) <- c('SampleID', 'waist_circum', 'hip_circum', 'sys_bp_mmHg', 'leukocyte', 'erythrocyte', 
                 'Hb_g_dec', 'mean_corp_vol', 'mean_corp_hb', 'mean_crop_hb_conc', 'platelet_count', 
                 'albumin', 'alanine_aminotransferase', 'aspartate aminotransferase', 'urea',
                 'cholesterol_mmol_L', 'creatinine_umol_L', 'CRP', 'cystatin_C', 'glucose',
                 'hba1c', 'HDL_mmol_L', 'LDL_mmol_L', 'triglycerides_mmol_L', 'urate_umol_L')

ID <- left_join(ID, covs, by = 'SampleID')

################################################################################### 

### Format medications info needed

library(readxl)
cat <- read_excel("/path_to_file.../Med_categories_UKB_3code.xlsx")
index <- read_excel("/path_to_file.../Med_table_index.xlsx")

cat <- as.data.frame(cat)
index <- as.data.frame(index)

nap <- index[which(index$Category %in% 'naproxen'),]

para <- index[which(index$Category %in% 'paracetamol'),]

# Join index into medications data by UKB coding
names(index)[2] <- 'Med'
names(index) <- c('Category', 'Med', 'ATC_code', 'Drug_name')
test <- left_join(ID, index, by = 'Med')

test$ATC <- substr(test$ATC, 1, 3) 

# ATC codes - C02 for antihyperintensives
hyp <- test[grep('C02', test$ATC),]
unique(hyp$Drug_name)

# Cumulative BP medications
bp <- test[which(test$ATC %in% c('C02', 'C03', 'C07', 'C08', 'C09')),]
# Website ATC: https://www.whocc.no/atc_ddd_index/

# Add lowering drugs into table of covariates, while also retaining NAs for individuals without med reporting

test$Med_available <- test$Med
test$Med_available <- ifelse(test$Med_available > 1, 0, NA)

# FALSE  TRUE
# 35073 12527

bp <- bp[c(1,47)]
names(bp)[2] <- 'bp_ATC'
bp$bp_linked <- 1

test <- left_join(test, bp, by = 'SampleID')

# Save out covariates for inclusion

write.csv(test, '/path_to_file.../covariates_additional.csv', row.names = F)

################################################################################### 

### Process covariates for inclusion in analyses (imputation for complete ProteinScore testing subset of the 47,600)

covs <- read.csv("/path_to_file.../covariates_joint.csv")
added <- read.csv("/path_to_file.../covariates_additional.csv")
join <- left_join(covs, added, by = 'SampleID')
join <- join[-c(9:39)]
join <- join[-c(6,7,8,10)]
PRS_Med <- join[c(1,7:21,23,48:56)]
Other <- join[c(1,2,5,24:47,3,4,6,22)]

## Check people's missingness for the cov measurements for imputation
olink_internal <- Other
length <- length(rownames(olink_internal))
res <- data.frame(SampleID = 1:length, Complete = 1:length, Missing = 1:length)

for (i in 1:length){
  variable <- as.character(olink_internal$SampleID[i])
  individual <- olink_internal[which(olink_internal$SampleID %in% variable),]
  individual <- individual[2:31]
  missing <- sum(is.na(individual))
  complete <- 30 - missing
  res[i,1] <- variable
  res[i,2] <-  complete
  res[i,3] <- missing
  print(i)
}

# Order by those that have the most missing data
res <- res[order(res$Complete),]
res$prop <- (res$Missing / 30) * 100
res$ex <- ifelse(res$prop >= 10, 1, 0)

write.csv(res, '/path_to_file.../missingness_people_covs.csv', row.names = F)

# Subset to remove NMAR individuals
exclude <- res[which(res$ex %in% 1),]
covs <- olink_internal[-which(olink_internal$SampleID %in% exclude$SampleID),]
  
# Check to ensure all variables have <10% missingness in measurements
join <- covs
res <- data.frame(Trait = 1:length(colnames(covs)),
                  Com = 1:length(colnames(covs)),
                  Miss = 1:length(colnames(covs)))

for(i in 1:length(colnames(join))){
  col <- join[i]
  missing <- sum(is.na(col))
  com <- 40257 - missing
  print(names(col))
  print(missing)
  res[i,1] <- names(col)
  res[i,2] <- com
  res[i,3] <- missing
}

res <- res[order(res$Miss),]
res$prop <- (res$Miss / 40257 ) * 100
res$ex <- ifelse(res$prop >= 10, 1, 0 )

# All fine - proceed
write.csv(res, '/path_to_file.../missing_added_covariates.csv', row.names = F)

# Impute categorical variables via median imputation
join$Alc[is.na(join$Alc)] <- median(join$Alc, na.rm = TRUE)
join$Smo[is.na(join$Smo)] <- median(join$Smo, na.rm = TRUE)
join$Ed[is.na(join$Ed)] <- median(join$Ed, na.rm = TRUE)
join$PA[is.na(join$PA)] <- median(join$PA, na.rm = TRUE)

# Impute continuous variables via knn
vars <- colnames(join)[2:27]
impute_knn <- join[which(colnames(join) %in% c('SampleID', vars))]

library(tidyverse)
library(impute)
IDs <- impute_knn[1]
data <- impute_knn[-1]
data <- as.matrix(data)
rownames(data) <- IDs$SampleID
data <- t(data)
print(dim(data))
imputed <- impute.knn(data)
imputed_data <- as.data.frame(t(imputed$data))
identical(as.character(rownames(imputed_data)), as.character(IDs$SampleID))
imputed_data <- cbind(IDs, imputed_data)

# Add categorical covs back in
imputed_data$Alc <- join$Alc
imputed_data$Smo <- join$Smo
imputed_data$Ed <- join$Ed
imputed_data$PA <- join$PA

# Save out untransformed
write.csv(imputed_data, '/path_to_file.../imputed_covs.csv', row.names = F)

# Save unimputed
write.csv(PRS_Med, '/path_to_file.../PRS_MED_unimputed_covs.csv', row.names = F)

# Transform continuous imputed covariates (log)
for(i in colnames(imputed_data)[3:27]){ 
  imputed_data[,i]<- log(imputed_data[,i])
}

# Save transformed 
write.csv(imputed_data, '/path_to_file.../imputed_covs_transformed.csv', row.names = F)



































































