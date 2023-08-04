###############################################################################################

## Sample allocations for Methylpiper models 

###############################################################################################

args <- commandArgs(trailingOnly=TRUE)
taskid <- as.numeric(args[1])
print(paste0('taskid for this iteration is ', taskid))

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

results <- list()
traits_list <- list()
location <- '/path_to_file.../tables_cut/'

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

###############################################################################################

### Load 50 pre-generated seeds randomly sampled between 1 and 5000
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

# Create output folder for disease with file subdirectory
location_out <- '/path_to_file.../ProteinScores/Run_210723/'
dir.create(paste0(location_out, name))
dir.create(paste0(location_out, name, '/files'))

### Isolate cases and controls available in the full sample for this disease
cox <- cox[which(cox$tte > 0),] 
names(cox)[15] <- 'time_to_event'
cases <- cox[which(cox$Event == '1'),]
controls <- cox[which(cox$Event == '0'),]

###############################################################################################

### Load imputed proteins from imputation prep script - (i.e. that have not been transformed or scaled with related excluded - transformations and scaling are done per train/test set below)
olink <- readRDS('/path_to_file.../knn_imputed_processed_olink_internal_proteins.rds')
print('Protein data loaded. Sampling train and test populations now.')

###############################################################################################

### Sample train and test samples and save down for recall

# Set ratio for case:controls
ratio <- 3

# Set number of proteins
N_prot <- dim(olink)[2]

# Assign threshold for follow-up testing of scores based on lists above 
if(name %in% list_10yr == TRUE){
  threshold = 10
} else {
  if(name %in% list_5yr == TRUE){
    threshold = 5
  } else {
    print(paste0('Threshold cannot be ascribed, due to trait name not being recognised.'))
  }
}

thr <- threshold 
print(paste0('Threshold set to ', thr, ' for test sampling.'))

# Create list for cvs
list <- list()

# Create table as summary of what has been sampled
record <- data.frame(iter = 1:50, total_cases = 1:50, total_controls = 1:50,
                     cases_train = 1:50, cases_test = 1:50, controls_train = 1:50, controls_test = 1:50,
                     cases_test_under = 1:50, cases_test_over = 1:50, controls_train_sampled = 1:50, controls_test_sampled = 1:50,
                     ytrain = 1:50, ytest = 1:50, recoded = 1:50,
                     seed = seed_list)

# Assign starting number of cases and controls
record[,2] <- dim(cases)[1]
record[,3] <- dim(controls)[1]

# Sample train and test populations x 50 iterations by random seed
for(i in 1:length(seed_list)){
  tryCatch({
    # Extract seed for this iteration
    seed <- seed_list[i]
    iteration <- i
    print(seed)
    set.seed(seed)

    # Randomly sample 50% cases
    case_IDs <- cases$SampleID
    set.seed(as.integer(seed))
    train_IDs <- sample(case_IDs, size=length(case_IDs)/2, replace=FALSE)
    cases_tr <- cases[which(cases$SampleID %in% train_IDs),]
    cases_te <- cases[-which(cases$SampleID %in% train_IDs),]

    record[i,4] <- dim(cases_tr)[1]
    record[i,5] <- dim(cases_te)[1]

    # Randomly sample 50% controls
    control_IDs <- controls$SampleID
    set.seed(as.integer(seed))
    control_train_IDs <- sample(control_IDs, size=length(control_IDs)/2, replace=FALSE)
    controls_tr <- controls[which(controls$SampleID %in% control_train_IDs),]
    controls_te <- controls[-which(controls$SampleID %in% control_train_IDs),]

    record[i,6] <- dim(controls_tr)[1]
    record[i,7] <- dim(controls_te)[1]

    # Exclude covariates with missing values from training data (will cause errors in MethypipeR if not)
    cases_tr <- cases_tr[-c(2:6)]
    controls_tr <- controls_tr[-c(2:6)]

    # Subset test cases to those within 10 year onset
    cases_te_under <- cases_te[which(cases_te$time_to_event <= thr),]

    # Subset cases with tte >10 years onset
    cases_te_over <- cases_te[which(cases_te$time_to_event > thr),]

    record[i,8] <- dim(cases_te_under)[1]
    record[i,9] <- dim(cases_te_over)[1]

    # Randomly sample 1:3 case:control ratio in train and test sets once
    Ncases1 <- dim(cases_tr)[1]
    Ncases2 <- dim(cases_te_under)[1]

    # Join cases beyond 10 year limit back into test control population
    controls_te <- rbind(cases_te_over, controls_te)

    length1 <- dim(controls_tr)[1]
    length2 <- dim(controls_te)[1]

    Ncontrols1 <- ratio * Ncases1
    Ncontrols2 <- ratio * Ncases2

    set.seed(as.integer(seed))
    sample1 <- sample(length1, size=Ncontrols1, replace=FALSE)
    set.seed(as.integer(seed))
    sample2 <- sample(length2, size=Ncontrols2, replace=FALSE)

    controls_train <- controls_tr[sample1,]
    controls_test <- controls_te[sample2,]

    record[i,10] <- dim(controls_train)[1]
    record[i,11] <- dim(controls_test)[1]

    y_train <- rbind(cases_tr, controls_train)
    y_test <- rbind(cases_te_under, controls_test)

    record[i,12] <- dim(y_train)[1]
    record[i,13] <- dim(y_test)[1]

    # Check how many cases > 10 years got selected as controls in test set
    record[i,14] <- length(which(controls_test$Event == '1' & controls_test$time_to_event > thr))[1]

    # Subset x files to train and test as inputs, ensuring order of IDs are matched to y files
    x_train <- olink[which(olink$pseudo_ind_id %in% y_train$pseudo_ind_id),]

    x_test <- olink[which(olink$pseudo_ind_id %in% y_test$pseudo_ind_id),]

    x_train <- x_train[match(y_train$pseudo_ind_id, x_train$pseudo_ind_id),]

    x_test <- x_test[match(y_test$pseudo_ind_id,  x_test$pseudo_ind_id),]

    # Sanity check on matching - printed to terminal
    identical(as.character(x_train$pseudo_ind_id), as.character(y_train$pseudo_ind_id))
    identical(as.character(x_test$pseudo_ind_id), as.character(y_test$pseudo_ind_id))

    x_train <- as.matrix(x_train[3:N_prot])
    x_test <- as.matrix(x_test[3:N_prot])

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
    saveRDS(x_train, paste0(location_out, name, '/files/x_train_', i, '.rds'))
    saveRDS(x_test, paste0(location_out, name, '/files/x_test_', i, '.rds'))
    saveRDS(y_train, paste0(location_out, name, '/files/y_train_', i, '.rds'))
    saveRDS(y_test, paste0(location_out, name, '/files/y_test_', i, '.rds'))

  }, error = function(e) cat("skipped"))
}

write.csv(record, paste0(location_out, name, '/record_sampling.csv'), row.names = F)


















