###############################################################################################

## Metabolomics UKB

###############################################################################################

# srun -p interactive --pty bash

# module load R

R

# Load packages
library(arrow)
library(data.table)
library(tidyverse)
library(ukbnmr)

ukbb_metab <-data.table::rbindlist(lapply(Sys.glob("/path_to_files.../Metabolomics/metabolomics/part-*.parquet"), arrow::read_parquet))
ukbb_metab$f.eid <- as.character(ukbb_metab$f.eid)

# Subset to olink proteins
olink_internal <- read.csv('/path_to_files.../knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')
ukbb_metab <- ukbb_metab[which(ukbb_metab$f.eid %in% olink_internal$SampleID),]

# Use ukbnmr package to extract data with annotations
nmr <- extract_biomarkers(ukbb_metab) 
ukbb_metab = NULL

# Read in eric extraction QC metab variables
fields <- read.delim("/path_to_files.../categories_221_222/fields.ukb", header = FALSE)
test <- read.delim("/path_to_files.../categories_221_222/ukb673344.txt", header = FALSE, sep = "\t")
names <- test[1,]
test <- test[-1]
names <- names[-243]
names(test) <- names[1,]
test <- test[-1,]
test <- test[,-242]
names(test) <- sub('-', '.', names(test))
names(test) <- paste0('f.', names(test))

# Subset to olink proteins
test <- test[which(test$f.eid %in% olink_internal$SampleID),]
olink_internal = NULL

# Extract fields: info needed: 23700-23948 and 23649-23660
biomarker_qc_flags <- extract_biomarker_qc_flags(test)
sample_qc_flags <- extract_sample_qc_flags(test)

# Save out the extracted files needed
fwrite(nmr, file="/path_to_files.../nmr_biomarker_data.csv")
fwrite(biomarker_qc_flags, file="/path_to_files.../eric_metab_extraction/nmr_biomarker_qc_flags.csv")
fwrite(sample_qc_flags, file="/path_to_files.../eric_metab_extraction/nmr_sample_qc_flags.csv")


### Load metabolite data, untransformed (12059 overlap with proteins)
metab <- read.csv('/path_to_files.../nmr_biomarker_data.csv')
metab <- metab[which(metab$visit_index %in% 0),]
length(which(metab[,1] %in% olink$SampleID))
names(metab)[1] <- 'SampleID'
metab <- metab[-2]

## Check people's missingness for the metab measurements
olink_internal <- metab
length <- length(rownames(olink_internal))
res <- data.frame(SampleID = 1:length, Complete = 1:length, Missing = 1:length)

for (i in 1:length){
  variable <- as.character(olink_internal$SampleID[i])
  individual <- olink_internal[which(olink_internal$SampleID %in% variable),]
  individual <- individual[1:250]
  missing <- sum(is.na(individual))
  complete <- 250 - missing
  res[i,1] <- variable
  res[i,2] <-  complete
  res[i,3] <- missing
  print(i)
}

# Order by those that have the most missing data
res <- res[order(res$Complete),]

# Index how many people have >10% missingness as a proportion
res$prop <- res$Missing / 250
res$exclude <- ifelse(res$prop > 0.1, '1', '0')

table(res$exclude)

keep <- res[which(res$exclude %in% 0),]
exclude <- res[which(res$exclude %in% 1),]
olink_internal <- olink_internal[which(olink_internal$SampleID %in% keep$SampleID),]


# Run missingness assessment on metabolite overlap
ukb <- as.data.frame(olink_internal)
length <- length(colnames(ukb))
res2 <- data.frame(Variable = 1:length, Complete = 1:length, Missing = 1:length)
names <- colnames(ukb)
for(i in 1:length(names)){
  variable <- as.character(names[i])
  complete <- ukb[which(complete.cases(ukb[,variable])),]
  incomplete <- ukb[-which(complete.cases(ukb[,variable])),]
  index_present <- dim(complete)[1]
  index_missing <- dim(incomplete)[1]
  res2[i,1] <- variable
  res2[i,2] <-  index_present
  res2[i,3] <- index_missing
  print(i)
}

res2$prop <- res2$Missing / res2$Complete
res2$exclude <- ifelse(res2$prop > 0.1, '1', '0')

# Run knn imputation for missing values

library(tidyverse)
library(impute)

olink_internal_sub2 <- ukb
IDs <- olink_internal_sub2[c(1)]
data <- olink_internal_sub2[c(2:250)]
data <- as.matrix(data)
rownames(data) <- IDs$SampleID
data <- t(data)
print(dim(data))

imputed <- impute.knn(data)
imputed_data <- as.data.frame(t(imputed$data))

identical(as.character(rownames(imputed_data)), as.character(IDs$SampleID))

imputed_data <- cbind(IDs, imputed_data)

write.csv(imputed_data, '/path_to_files.../knn_imputed_processed_metab_data.csv', row.names = F)
dim(imputed_data)







