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
olink_internal <- read.csv('knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')

# Extraction (v1)
dat <- read.csv("parquet_47600_040723.csv")
dat <- dat[which(dat$f.eid %in% olink_internal$SampleID),]

not <- olink_internal[-which(olink_internal$SampleID %in% dat$f.eid),]

covs <- dat[,which(colnames(dat) %in% c('f.eid','f.189.0.0', 'f.20116.0.0', 'f.1558.0.0', 'f.21001.0.0'))]
names(covs) <- c('SampleID', 'Dep', 'Alc', 'Smo', 'BMI')

# Plot distributions
pdf('BMI_unadjusted.pdf')
ggplot(covs, aes(x=BMI)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 
dev.off() # 229 missing

pdf('Dep_unadjusted.pdf')
ggplot(covs, aes(x=Dep)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
dev.off() # 59 missing

pdf('Smo_unadjusted.pdf')
ggplot(covs, aes(x=Smo)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
dev.off() # 57 missing

pdf('Alc_unadjusted.pdf')
ggplot(covs, aes(x=Alc)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
dev.off() # 57 missing

# Load in additional covariates from ben GWAS
cov <- fread("combined_covars_v1.tsv")
cov <- cov[which(cov$IID %in% olink_internal$SampleID),]
table(is.na(cov$UKBPC_16)) # FALSE- 48821 with genetic data
table(cov$ukb_centre_fct)

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
t <- read.csv("eric_individual_fieldIDs_extracted.csv")
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
table(ed$Edu)
ed <- ed[c(1,26)]
joint <- left_join(joint, ed, by = 'SampleID')
write.csv(joint, 'covariates_joint.csv', row.names = F)

#################################################

# ### Add additional covariates extracted during revision
# Load phenotypes file
phen <- read.table(gzfile("all_quant_20200904.pheno.v2.tsv.gz"), header = T)
names(phen)[1] <- 'SampleID'
phen <- phen[which(phen$SampleID %in% olink_internal$SampleID),]

# Save both sets of covariate files out 
write.csv(phen, 'all_quant_20200904.pheno.csv', row.names = F)
write.csv(joint, 'f', row.names = F)




























































