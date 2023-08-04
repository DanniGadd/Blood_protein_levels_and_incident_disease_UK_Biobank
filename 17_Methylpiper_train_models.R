###############################################################################################

## Train scores

###############################################################################################

### Load sbatch information when run as array per disease
args <- commandArgs(trailingOnly=TRUE)
taskid <- as.numeric(args[1])
print(paste0('taskid for this iteration is ', taskid))

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

results <- list()
traits_list <- list()
location <- '/path_to_file.../Censor_test/tables_cut/'

files <- list.files(location, '.csv')

for(i in 1:length(files)){
  name <- files[i]
  name <- gsub("\\..*", "", name)  
  file <- read.csv(paste0(location, name, '.csv'))
  traits_list[[i]] <- name
  results[[i]] <- file
}

# Set the array iteration for disease selection
i <- taskid
format <- sub(".csv", "", files)

### Define the traits that will be assessed over 5yr and 10yr follow up 
list_5yr <- c( 'ALS_FO', 'CYS_FO', 'ENDO_FO')
list_10yr <- format[-which(format %in% list_5yr)]
i <- taskid

###############################################################################################

### Load 50 seeds randomly sampled between 1 and 5000
# seed <- 4237
# set.seed(seed)
# seed_list <- floor(runif(50, min=1, max=5000))
# write.csv(seed_list, '/path_to_file.../seed_list_50_sampled_with_seed_4237.csv')
seed_list <- read.csv('/path_to_file.../seed_list_50_sampled_with_seed_4237.csv')
seed_list <- seed_list$x

###############################################################################################

### Set the disease for this array 
name <- as.character(traits_list[[i]])
cox <- results[[i]]
print(paste0('The disease trait for this iteration is ', name))

# # Create output folder for disease with file subdirectory
location_out <- '/path_to_file.../ProteinScores/Run_210723/'
dir.create(paste0(location_out, name, '/model_logs'))

### Isolate cases and controls available in the full sample for this disease
cox <- cox[which(cox$tte > 0),] 
names(cox)[15] <- 'time_to_event'
cases <- cox[which(cox$Event == '1'),]
controls <- cox[which(cox$Event == '0'),]

length <- dim(cases)[1]
length <- length/2

# Assign fold number based on availability of cases in training sample
if(length >= 1000){
  cv = 10
} else {
  if(length < 1000 & length >= 500){
    cv = 5
  } else {
    cv = 3
  }
}
cross <- cv
print(paste0('Using ', cross, ' cross-fold validation.'))

###############################################################################################

### Load imputed proteins from imputation prep script - (i.e. that have not been transformed or scaled with related excluded - transformations and scaling are done per train/test set below)
olink <- readRDS('/path_to_file.../knn_imputed_processed_olink_internal_proteins.rds')
print('Protein data loaded. Looping analyses through seeds now.')

###############################################################################################

### LOAD TRAIN MODEL FUNCTION
train <- function(trainingData, trainingTarget, seed_iteration, cv, iter, name_dis, output){
  tryCatch({

    model_logs <- paste0(output, name_dis, '/model_logs/', iter, '/')
    dir.create(paste0(model_logs))
    initLogs(model_logs, note = paste0('Model for ', name, ' and iteration ', iter, '.'))
    set.seed(seed_iteration)

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
    write.csv(res, paste0(output, name_dis, '/files/weights_', iter, '.csv'), row.names = F)
    print(dim(res))

    # save model
    saveRDS(mprModel, paste0(output, name_dis, '/files/model_', iter, '.rds'))

    sessionStartTimestamp <- getOption('mprSessionStartTimestamp')
    
    gc()
    sink()
  }, error = function(e) cat("skipped"))
}


###############################################################################################

### Run training and save weights for each iteration

print(paste0('Allocations prepped - trying model runs.'))

list_weights <- list()

  for(i in 34:50){
    seed <- seed_list[i]
    iteration <- i
    print(i)

    x_train <- readRDS(paste0(location_out, name, '/files/x_train_', i, '.rds'))
    y_train <- readRDS(paste0(location_out, name, '/files/y_train_', i, '.rds'))

    train(x_train, y_train, seed, cross, iteration, name, location_out)
  }

print(paste0('Model runs complete.'))


























