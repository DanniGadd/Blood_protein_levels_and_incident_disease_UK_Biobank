###############################################################################################

## Metabolomics UKB - scoring analyses

###############################################################################################

# # ### Load sbatch information when run as array per disease
# args <- commandArgs(trailingOnly=TRUE)
# taskid <- as.numeric(args[1])
# print(paste0('taskid for this iteration is ', taskid))

### Interactive session and R required for script troubleshooting
# srun -p interactive --pty bash
# module load R/4.2.0-foss-2021b
# R

### Load packages
library(readxl)
library(tidyverse)
library(MethylPipeR)
library(glmnet)
library(survival)
library(gbm)
library(data.table)
library(pacman)

###############################################################################################

### Use the trait cox tables from the incident disease summary table script:

# Read subset files back in 
results <- list()
traits_list <- list()
location <- 'Results/Cox/tables_cut/'
files <- list.files(location, '.csv')

for(i in 1:length(files)){
  name <- files[i]
  name <- gsub("\\..*", "", name)  
  file <- read.csv(paste0(location, name, '.csv'))
  traits_list[[i]] <- name
  results[[i]] <- file
  # print(name)
  # print(dim(file))
}

# Set the array iteration for disease selection
i <- 8
format <- sub(".csv", "", files)

### Define the traits that will be assessed over 5yr and 10yr follow up 
list_5yr <- c('ALS', 'Brain_FO', 'Dep', 'LUP')
list_10yr <- format[-which(format %in% list_5yr)]

### Set the disease for this array 
name <- as.character(traits_list[[i]])
cox <- results[[i]]
print(paste0('The disease trait for this iteration is ', name))

### Isolate cases and controls available in the full sample for this disease
cox <- cox[which(cox$tte > 0),] # ensure no negative tte between 0 and -1 for predictor models
names(cox)[15] <- 'time_to_event'
cases <- cox[which(cox$Event == '1'),]
controls <- cox[which(cox$Event == '0'),]

###############################################################################################

### Load imputed proteins from imputation prep script - (i.e. that have not been transformed or scaled with related excluded - transformations and scaling are done per train/test set below)
olink <- readRDS('knn_imputed_processed_olink_internal_proteins.rds')

### Load metabolite data, untransformed (12059 overlap with proteins)
metab <- read.csv('nmr_biomarker_data.csv')
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
write.csv(imputed_data, 'knn_imputed_processed_metab_data.csv', row.names = F)
dim(imputed_data)

# Read in metabolite imputed data
metab <- read.csv('/knn_imputed_processed_metab_data.csv')

# Subset olink to 12059 sample
olink <- olink[which(olink$SampleID %in% metab$SampleID),]
   
# Assess variables in relation to metab data (as per protein data assessment)
covs <- read.csv('covariates_joint.csv')
covs <- covs[complete.cases(covs$Batch),]

# metab <- olink_internal_sub2[which(olink_internal_sub2$SampleID %in% covs$SampleID),]
## Regress Proteins onto Covariates
phenotypes_residualised <- left_join(metab, covs, by = 'SampleID')
for(i in colnames(phenotypes_residualised)[2:250]){
  phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ factor(ukb_centre_fct),
                                   na.action = na.exclude, data = phenotypes_residualised)$residuals
}

# Correlate residuals against protein levels
met <-metab[which(metab$SampleID %in% phenotypes_residualised$SampleID),]
length <- c(2:250)
names <- colnames(met)[2:250]
res <- data.frame(MET = 1:3, Coef = 1:3, p = 1:3)
identical(colnames(met[2:250]), colnames(phenotypes_residualised[,2:250]))

for (i in 1:249){
  metabolite <- as.character(names[i])
  cor <- cor.test(met[,metabolite], phenotypes_residualised[,metabolite], na.rm = T)
  r <- cor$estimate
  p <- cor$p.value
  res[i,1] <- metabolite
  res[i,2] <- r
  res[i,3] <- p
}

res <- res[order(res$Coef),]
write.csv(res, 'correlations_pre_post_resid.csv', row.names = F)

# Create a joint file with both sets of markers
joint <- left_join(olink, metab, by = 'SampleID')

# Subset to diabetes test set sample and join metabolomics in
library(tidyverse)
library(survival)
library(gbm)
library(precrec)
library(ggplot2)

# Load in additional covariates required in rerun
merge <- read.csv('imputed_covs_transformed.csv')

### Define diseases from chosen models
location_out <- "Results/Cox/Proteinscores/"
core_models_output <- 'Results/Cox/Score_processing/covariate_assessment/'
chosen <- read.csv("Results/Cox/Score_processing/models_chosen_accounted_features.csv")
location_out <- "Results/Cox/Proteinscores/"

seed <- 2103
set.seed(seed)
iteration <- 13
threshold <- 10

# Load test data for specific disease iteration
y_train <- readRDS("Proteinscores/DEATH/files/y_train_13.rds")
y_test <- readRDS("Proteinscores/DEATH/files/y_test_13.rds")

x_train <- readRDS("Proteinscores/DEATH/files/x_train_13.rds")
x_test <- readRDS("Proteinscores/DEATH/files/x_test_13.rds")

# Merge in added covariates to y test 
y_test <- y_test[-c(2:6)]
y_test <- left_join(y_test, merge, by = 'SampleID')

# x training 
olink2 <- olink[which(olink$SampleID %in% y_train$SampleID),]
metab2 <- metab[which(metab$SampleID %in% y_train$SampleID),]
joint2 <- joint[which(joint$SampleID %in% y_train$SampleID),]

# x testing
olink3 <- olink[which(olink$SampleID %in% y_test$SampleID),]
metab3 <- metab[which(metab$SampleID %in% y_test$SampleID),]
joint3 <- joint[which(joint$SampleID %in% y_test$SampleID),]

# y files subset from cox, then subset to metab overlap
train <- cox[which(cox$SampleID %in% y_train$SampleID),]
test <- cox[which(cox$SampleID %in% y_test$SampleID),]
train <- train[which(train$SampleID %in% olink2$SampleID),]
test <- test[which(test$SampleID %in% olink3$SampleID),]

which(train$SampleID %in% test$SampleID)# 0

### Subset cases and controls to the population overlap 
cases_train <- cases[which(cases$SampleID %in% y_train$SampleID),]
cases_test <- cases[which(cases$SampleID %in% y_test$SampleID),]
cases_tr <- cases_train[which(cases_train$SampleID %in% olink2$SampleID),]
cases_te <- cases_test[which(cases_test$SampleID %in% olink3$SampleID),]

#################

### Proteins

y_train <- train
y_test <- test
x_train <- olink2
x_test <- olink3

# match files x and y 
x_train <- x_train[match(y_train$SampleID, x_train$SampleID),]
x_test <- x_test[match(y_test$SampleID,  x_test$SampleID),]

# Sanity check on matching - printed to terminal
identical(as.character(x_train$SampleID), as.character(y_train$SampleID))
identical(as.character(x_test$SampleID), as.character(y_test$SampleID))

x_train <- as.matrix(x_train[3:1468])
x_test <- as.matrix(x_test[3:1468])

# Transform and scale protein inputs in train and test populations
for(k in colnames(x_train)){
  x_train[,i] <- qnorm((rank(x_train[,k], na.last='keep')-0.5)/sum(!is.na(x_train[,k])))
}

x_train <- apply(x_train, 2, scale)

for(k in colnames(x_test)){
  x_test[,k] <- qnorm((rank(x_test[,k], na.last='keep')-0.5)/sum(!is.na(x_test[,k])))
}

x_test <- apply(x_test, 2, scale)

# Save train and test files with iteration identifiers
location_out1 <- 'Results/Cox/metab_analysis/Death/prot/'
saveRDS(x_train, paste0(location_out1, '/files/x_train_', '13', '.rds'))
saveRDS(x_test, paste0(location_out1, '/files/x_test_', '13', '.rds'))
saveRDS(y_train, paste0(location_out1, '/files/y_train_', '13', '.rds'))
saveRDS(y_test, paste0(location_out1, '/files/y_test_', '13', '.rds'))

###############################################################################################

# Metab

y_train <- train
y_test <- test
x_train <- metab2
x_test <- metab3

# match files x and y 
x_train <- x_train[match(y_train$SampleID, x_train$SampleID),]
x_test <- x_test[match(y_test$SampleID,  x_test$SampleID),]

# Sanity check on matching - printed to terminal
identical(as.character(x_train$SampleID), as.character(y_train$SampleID))
identical(as.character(x_test$SampleID), as.character(y_test$SampleID))

x_train <- as.matrix(x_train[2:250])
x_test <- as.matrix(x_test[2:250])

# Transform and scale protein inputs in train and test populations
for(k in colnames(x_train)){
  x_train[,i] <- qnorm((rank(x_train[,k], na.last='keep')-0.5)/sum(!is.na(x_train[,k])))
}

x_train <- apply(x_train, 2, scale)

for(k in colnames(x_test)){
  x_test[,k] <- qnorm((rank(x_test[,k], na.last='keep')-0.5)/sum(!is.na(x_test[,k])))
}

x_test <- apply(x_test, 2, scale)

# Save train and test files with iteration identifiers
location_out1 <- 'Results/Cox/metab_analysis/Death/metab/'
saveRDS(x_train, paste0(location_out1, '/files/x_train_', '13', '.rds'))
saveRDS(x_test, paste0(location_out1,  '/files/x_test_', '13', '.rds'))
saveRDS(y_train, paste0(location_out1,  '/files/y_train_', '13', '.rds'))
saveRDS(y_test, paste0(location_out1,  '/files/y_test_', '13', '.rds'))

###############################################################################################

# Joint

y_train <- train
y_test <- test

x_train <- joint2
x_test <- joint3

# match files x and y 
x_train <- x_train[match(y_train$SampleID, x_train$SampleID),]
x_test <- x_test[match(y_test$SampleID,  x_test$SampleID),]

# Sanity check on matching - printed to terminal
identical(as.character(x_train$SampleID), as.character(y_train$SampleID))
identical(as.character(x_test$SampleID), as.character(y_test$SampleID))
x_train <- as.matrix(x_train[3:1719])
x_test <- as.matrix(x_test[3:1719])

# Transform and scale protein inputs in train and test populations
for(k in colnames(x_train)){
  x_train[,i] <- qnorm((rank(x_train[,k], na.last='keep')-0.5)/sum(!is.na(x_train[,k])))
}

x_train <- apply(x_train, 2, scale)

for(k in colnames(x_test)){
  x_test[,k] <- qnorm((rank(x_test[,k], na.last='keep')-0.5)/sum(!is.na(x_test[,k])))
}

x_test <- apply(x_test, 2, scale)

# Save train and test files with iteration identifiers
location_out1 <- 'Results/Cox/metab_analysis/Death/joint/'
saveRDS(x_train, paste0(location_out1, '/files/x_train_', '13', '.rds'))
saveRDS(x_test, paste0(location_out1,  '/files/x_test_', '13', '.rds'))
saveRDS(y_train, paste0(location_out1,  '/files/y_train_', '13', '.rds'))
saveRDS(y_test, paste0(location_out1,  '/files/y_test_', '13', '.rds'))

#####################################################################################

### Run training tests for each omics type 

### LOAD TRAIN MODEL FUNCTION
train <- function(trainingData, trainingTarget, seed_iteration, cv, iter, name_dis, output){
  tryCatch({
    # Create iteration specific output folder with results subfolder
    model_logs <- paste0(output, '/model_logs/', iter, '/')
    dir.create(paste0(model_logs))
    
    initLogs(model_logs, note = paste0('Model for ', name, ' and iteration ', iter, '.'))
    
    set.seed(seed_iteration)
    
    trainingTarget <- trainingTarget[-c(2:6)]
    
    # Fit the training cox glmnet
    mprModel <- fitMPRModelCV(type = 'survival',
                              method = 'glmnet',
                              trainXs = trainingData,
                              trainY = trainingTarget,
                              seed = seed_iteration,
                              nFolds = cv,
                              save = TRUE)
    
    # Extract and save weights for selected proteins
    model_summary <-  mprModel$model
    res <- coef(mprModel$model, s= 'lambda.min')
    res <- res[res[,1]!=0,]
    res <- as.data.frame(res)
    names(res)[1] <- 'Coefficient'
    res$Protein <- rownames(res)
    res <- res[c(2,1)]
    
    # save weights
    write.csv(res, paste0(output, '/files/weights_', iter, '_rerun_3cv.csv'), row.names = F)
    print(dim(res))
    
    # save model
    saveRDS(mprModel, paste0(output, '/files/model_', iter, '_rerun_3cv.rds'))
    
    sessionStartTimestamp <- getOption('mprSessionStartTimestamp')
    
    gc()
    sink()
  }, error = function(e) cat("skipped"))
}

###############################################################################################

# Protein training
location_out1 <- 'Results/Cox/metab_analysis/Death/prot/'
x_train <- readRDS(paste0(location_out1, '/files/x_train_', '18', '.rds'))
y_train <- readRDS(paste0(location_out1, '/files/y_train_', '18', '.rds'))

cv <- 5
cross <- 5
seed <- 2103
train(x_train, y_train, seed, cross, iteration, name, location_out1)

###############################################################################################

# Metab training

location_out1 <- 'Results/Cox/metab_analysis/Death/metab/'
x_train <- readRDS(paste0(location_out1, '/files/x_train_', '13', '.rds'))
y_train <- readRDS(paste0(location_out1, '/files/y_train_', '13', '.rds'))

cv <- 3
cross <- 3
seed <- 2103
train(x_train, y_train, seed, cross, iteration, name, location_out1)

###############################################################################################

# Joint training

location_out1 <- 'Results/Cox/metab_analysis/Death/joint/'
x_train <- readRDS(paste0(location_out1, '/files/x_train_', '18', '.rds'))
y_train <- readRDS(paste0(location_out1, '/files/y_train_', '18', '.rds'))

cv <- 5
cross <- 5
seed <- 2103
train(x_train, y_train, seed, cross, iteration, name, location_out1)

###############################################################################################



















