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
location <- '01_paper/Results/Cox/tables_cut/'

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
i <- 10

format <- sub(".csv", "", files)

### Define the traits that will be assessed over 5yr and 10yr follow up 
list_5yr <- c('ALS', 'Brain_FO', 'Dep', 'LUP')
list_10yr <- format[-which(format %in% list_5yr)]

##########################################################################

### Set the disease for this array 
name <- as.character(traits_list[[i]])
cox <- results[[i]]
print(paste0('The disease trait for this iteration is ', name))

### Isolate cases and controls available in the full sample for this disease
cox <- cox[which(cox$tte > 0),] # ensure no negative tte between 0 and -1 for predictor models
names(cox)[15] <- 'time_to_event'
cases <- cox[which(cox$Event == '1'),]
controls <- cox[which(cox$Event == '0'),]


### Load imputed proteins from imputation prep script - (i.e. that have not been transformed or scaled with related excluded - transformations and scaling are done per train/test set below)
olink <- readRDS('PPP_core_input_files/knn_imputed_processed_olink_internal_proteins.rds')

# Read in metabolite imputed data (from diab script)
metab <- read.csv('Revision_analyses/PPP_core_input_files/knn_imputed_processed_metab_data.csv')

# Subset olink to 12059 sample
olink <- olink[which(olink$SampleID %in% metab$SampleID),]

# Create a joint file with both sets of markers
joint <- left_join(olink, metab, by = 'SampleID')

# Subset to diabetes test set sample and join metabolomics in
library(tidyverse)
library(survival)
library(gbm)
library(precrec)
library(ggplot2)

# Load in additional covariates required in rerun
merge <- read.csv('Covariate_preps/imputed_covs_transformed.csv')

### Define diseases from chosen models

location_out <- "00_Run_260823/01_paper/Results/Cox/Proteinscores/"
core_models_output <- '00_Run_260823/01_paper/Results/Cox/Score_processing/covariate_assessment/'

chosen <- read.csv("00_Run_260823/01_paper/Results/Cox/Score_processing/models_chosen_accounted_features.csv")

location_out <- "00_Run_260823/01_paper/Results/Cox/Proteinscores/"

seed <- 3104
set.seed(seed)
iteration <- 18
threshold <- 10

# Load test data for specific disease iteration
y_train <- readRDS("01_paper/Results/Cox/Proteinscores/Diab_FO/files/y_train_18.rds")
y_test <- readRDS("01_paper/Results/Cox/Proteinscores/Diab_FO/files/y_test_18.rds")

x_train <- readRDS("01_paper/Results/Cox/Proteinscores/Diab_FO/files/x_train_18.rds")
x_test <- readRDS("01_paper/Results/Cox/Proteinscores/Diab_FO/files/x_test_18.rds")

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
location_out1 <- '01_paper/Results/Cox/metab_analysis/Diab/prot/'
saveRDS(x_train, paste0(location_out1, '/files/x_train_', '18', '.rds'))
saveRDS(x_test, paste0(location_out1, '/files/x_test_', '18', '.rds'))
saveRDS(y_train, paste0(location_out1, '/files/y_train_', '18', '.rds'))
saveRDS(y_test, paste0(location_out1, '/files/y_test_', '18', '.rds'))

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
location_out1 <- '01_paper/Results/Cox/metab_analysis/Diab/metab/'
saveRDS(x_train, paste0(location_out1, '/files/x_train_', '18', '.rds'))
saveRDS(x_test, paste0(location_out1,  '/files/x_test_', '18', '.rds'))
saveRDS(y_train, paste0(location_out1,  '/files/y_train_', '18', '.rds'))
saveRDS(y_test, paste0(location_out1,  '/files/y_test_', '18', '.rds'))


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
location_out1 <- '01_paper/Results/Cox/metab_analysis/Diab/joint/'
saveRDS(x_train, paste0(location_out1, '/files/x_train_', '18', '.rds'))
saveRDS(x_test, paste0(location_out1,  '/files/x_test_', '18', '.rds'))
saveRDS(y_train, paste0(location_out1,  '/files/y_train_', '18', '.rds'))
saveRDS(y_test, paste0(location_out1,  '/files/y_test_', '18', '.rds'))

#####################################################################################

### Run training tests for each omics type 

### LOAD TRAIN MODEL FUNCTION
train <- function(trainingData, trainingTarget, seed_iteration, cv, iter, name_dis, output){
  tryCatch({
    
    # trainingData <- x_train
    # trainingTarget <- y_train
    # seed_iteration <- seed
    # cv <- cross
    # iter <- iteration
    # name_dis <- name
    # output <- location_out1
    
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
    write.csv(res, paste0(output, '/files/weights_', iter, '.csv'), row.names = F)
    print(dim(res))
    
    # save model
    saveRDS(mprModel, paste0(output, '/files/model_', iter, '.rds'))
    
    sessionStartTimestamp <- getOption('mprSessionStartTimestamp')
    
    gc()
    sink()
  }, error = function(e) cat("skipped"))
}

###############################################################################################

# Protein training

location_out1 <- '01_paper/Results/Cox/metab_analysis/Diab/prot/'
x_train <- readRDS(paste0(location_out1, '/files/x_train_', '18', '.rds'))
y_train <- readRDS(paste0(location_out1, '/files/y_train_', '18', '.rds'))

cv <- 3
cross <- 3
seed <- 3104

train(x_train, y_train, seed, cross, iteration, name, location_out1)

###############################################################################################

# Metab training

location_out1 <- '01_paper/Results/Cox/metab_analysis/Diab/metab/'
x_train <- readRDS(paste0(location_out1, '/files/x_train_', '18', '.rds'))
y_train <- readRDS(paste0(location_out1, '/files/y_train_', '18', '.rds'))

cv <- 3
cross <- 3
seed <- 3104

train(x_train, y_train, seed, cross, iteration, name, location_out1)

###############################################################################################

# Joint training

location_out1 <- '01_paper/Results/Cox/metab_analysis/Diab/joint/'
x_train <- readRDS(paste0(location_out1, '/files/x_train_', '18', '.rds'))
y_train <- readRDS(paste0(location_out1, '/files/y_train_', '18', '.rds'))

cv <- 3
cross <- 3
seed <- 3104

train(x_train, y_train, seed, cross, iteration, name, location_out1)


















