###############################################################################################

## Select ProteinScores from 50 iterations 

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

# Load in additional covariates required in rerun
merge <- read.csv('/path_to_file.../imputed_covs_transformed.csv')

###############################################################################################

### Define diseases

i <- taskid

location_out <- '/path_to_file.../ProteinScores/Run_210723/'

# Get disease names
location <- '/path_to_file.../tables_cut/'
files <- list.files(location, '.csv')
diseases <- sub(".csv", "", files)

### Set the disease for this array 
name <- as.character(diseases[[i]])
print(paste0('The disease trait for this iteration is ', name))

###############################################################################################

### Load 50 seeds randomly sampled between 1 and 5000
# seed <- 4237
# set.seed(seed)
# seed_list <- floor(runif(50, min=1, max=5000))
# write.csv(seed_list, '/path_to_file.../seed_list_50_sampled_with_seed_4237.csv')
seed_list <- read.csv('/path_to_file.../seed_list_50_sampled_with_seed_4237.csv')
seed_list <- seed_list$x

###############################################################################################

### LOAD TEST MODEL FUNCTION
test <- function(mprModel, testData, testTarget, res_list, seed_iteration, iter, name_dis, output, thr){
  tryCatch({
    
    set.seed(seed_iteration)

    # Generate model predictions in the test data
    mprModelTestPredictions <- predictMPRModel(mprModel,
                                               data = testData,
                                               s = 'lambda.min')
    
    # Add test scores to the target y file in the test sample
    testTarget$dScore <- mprModelTestPredictions
    testTarget <- na.omit(testTarget)
    test_count <- dim(testTarget)[1]
    cases_total <- length(which(testTarget$Event == 1))
    cases_recoded <- length(which(testTarget$time_to_event > thr & testTarget$Event == 1))
  
    # Compare null and scores added cox models
    if(name_dis %in% c('Prostate_FO', 'CYS_FO', 'GYN_FO', 'Breast_FO', 'ENDO_FO')){
      riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment, testTarget)
      fullCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + dScore, testTarget)
      ProteinOnly <- coxph(Surv(time_to_event, Event) ~ dScore, testTarget)
      print('Model has been sex stratified in cox.')
      print('Cases and controls available:')
      print(table(testTarget$Event))
    } else {
      riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex), testTarget)
      fullCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + dScore, testTarget)
      ProteinOnly <- coxph(Surv(time_to_event, Event) ~ dScore, testTarget)
      print('Model has not been sex stratified in cox.')
      print('Cases and controls available:')
      print(table(testTarget$Event))
    }
    
    # List models
    models <- list(r = riskFactorsOnlyCoxPH, f = fullCoxPH, t = ProteinOnly)
    
    ### For every Cox PH run in the test sample, calculate 10-year performance
    predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = thr) {
      
      # Extract cumulative baseline hazard at 10 year threshold
      uniqueTimes <- uniqueTimes <- sort(unique(c(dataDF$time_to_event, threshold)))
      thresholdIndex <- match(threshold, uniqueTimes)
      cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$time_to_event, dataDF$Event, predict(coxPHModel), uniqueTimes)
      
      # Calculate onset predictions
      survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
      onsetPredictions <- 1 - survivalPredictions
      
      # Track whether cumulative baseline hazard could be extracted in model record
      print('Cumulative baseline hazard presence NA - print onset predictions:')
      print(table(is.na(onsetPredictions)))
      
      # Recode events for AUC calculaion- event should be 0 if tte is > 10
      dataDF$Event <- sapply(1:nrow(dataDF), function(i) {
        if (dataDF$time_to_event[[i]] > threshold) {
          0
        } else {
          dataDF$Event[[i]]
        }
      })
      
      # Extract performance metrics for 10-year case classification for models
      auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$Event)
      prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$Event)
      roc <- pROC::roc(response = dataDF$Event, predictor = onsetPredictions)
      list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
    }
    
    # Extract test results
    testResults <- lapply(models, function(m) {predictCoxPHOnset(testTarget, m)})
    
    # Extract metrics table
    aucs <- sapply(testResults, function(r) {r$auc})
    praucs <- sapply(testResults, function(r) {r$prauc})
    metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)
    metricsTable[4,1] <- metricsTable[2,1] - metricsTable[1,1]
    metricsTable[4,2] <- metricsTable[2,2] - metricsTable[1,2]
    row.names(metricsTable) <- c('Null', 'Null + ProteinScore', 'ProteinScore Only', 'Difference')
    metricsTable$seed <- seed_iteration
    metricsTable$iteration <- iter
    metricsTable$test_total <- test_count
    metricsTable$test_cases <- cases_total
    metricsTable$test_cases_recoded <- cases_recoded
    
    # Save cox test results
    saveRDS(testTarget, paste0(output, name_dis, '/files/testResults_', iter, '.rds'))
    
    # save metrics table performances
    saveRDS(metricsTable, paste0(output, name_dis, '/files/metricsTable_', iter, '.rds'))
    
  }, error = function(e) cat("skipped"))
  return(metricsTable)
}


###############################################################################################

### Define threshold per trait 

location <- '/path_to_file.../tables_cut/'

files <- list.files(location, '.csv')
format <- sub(".csv", "", files)

### Define the traits that will be assessed over 5yr and 10yr follow up 
list_5yr <- c('ALS_FO', 'CYS_FO', 'ENDO_FO')
list_10yr <- format[-which(format %in% list_5yr)]

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

thr_input <- threshold 
print(paste0('Threshold set to ', thr_input, ' for test sampling.'))

###############################################################################################

### Run testing

list <- list()

set <- 1:50

for(i in set){
  tryCatch({
    seed <- seed_list[i]
    iteration <- i
    
    x_test <- readRDS(paste0(location_out, name, '/files/x_test_', i, '.rds'))
    y_test <- readRDS(paste0(location_out, name, '/files/y_test_', i, '.rds'))
    
    # Merge in added covariates
    y_test <- y_test[-c(2:6)]
    y_test <- left_join(y_test, merge, by = 'SampleID')
    
    model_load <- readRDS(paste0(location_out, name, '/files/model_', i, '.rds'))
    
    mod <- test(model_load, x_test, y_test, list, seed, iteration, name, location_out, thr_input)
    list[[i]] <- mod
  }, error = function(e) cat("skipped"))
}

result <- do.call(rbind, list)
write.csv(result, paste0(location_out, name, '/metricsTables_joint.csv'))


###############################################################################################

# Select proteinscores for each trait

# Get disease names
location <- '/path_to_file.../tables_cut/'
files <- list.files(location, '.csv')
diseases <- sub(".csv", "", files)

# Set location where results are collated
location_out <- "/path_to_file.../Run_210723/"

# Check available models for each that can be used
for(i in 1:length(diseases)){
  tryCatch({
    name <- diseases[i]
    print(name)
    path <- paste0(location_out, name, '/')

    files <- read.csv(paste0(path, 'metricsTables_joint.csv'))
    print(length(unique(files$iteration)))
  }, error = function(e) cat("skipped"))
}

diseases <- c('AL_FO', 'ALS_FO', 'Breast_FO', 'CYS_FO', 'Colorectal_FO', 'COPD_FO', 'ENDO_FO','DEATH', 'Diab_FO', 'GYN_FO',
              'IBD_FO', 'IHD_FO', 'LIV_FO', 'LUNG_FO', 'PD_FO', 'Prostate_FO', 'RA_FO', 'ST_FO', 'VD_FO')


# Find number of features selected for each weights file for each trait
list <- list()

for(i in 1:length(diseases)){
  tryCatch({
    name <- diseases[i]
    print(name)
    path <- paste0(location_out, name, '/files/')

    for(m in 1:50){
      iteration <- m
      file <- read.csv(paste0(path, 'weights_', iteration, '.csv'))
      print(dim(file))
      list[[m]] <- as.numeric(dim(file)[1])
    }

    trait <- do.call(rbind, list)
    write.csv(trait, paste0('/path_to_file.../Run_210723_processing/Weights_summaries/', name, '_V2.csv'), row.names = F)
  }, error = function(e) cat("skipped"))
}

# Identify instances where no features were selected (i.e. dimensions of weights = 1)

for(i in 1:length(diseases)){
  tryCatch({
    name <- diseases[i]
    print(name)
    path <- '/path_to_file.../Run_210723_processing/Weights_summaries/'

    file <- read.csv(paste0(path, name, '.csv'))
    print(length(which(file$V1 == 1)))
  }, error = function(e) cat("skipped"))
}


# Calculate median difference for models
list_chosen <- list()
diseases <- c('AL_FO', 'ALS_FO', 'Breast_FO', 'CYS_FO', 'Colorectal_FO', 'COPD_FO', 'ENDO_FO','DEATH', 'Diab_FO', 'GYN_FO',
              'IBD_FO', 'IHD_FO', 'LIV_FO', 'LUNG_FO', 'PD_FO', 'Prostate_FO', 'RA_FO', 'ST_FO', 'VD_FO')
extra <- c(0,4,1,0,1,0,0,0,0,3,0,0,0,0,0,0,0,0,0)
for(i in 1:length(diseases)){
  tryCatch({
    name <- diseases[i]
    print(name)
    path <- paste0(location_out, name, '/')

    files <- read.csv(paste0(path, 'metricsTables_joint.csv'))
    print(length(unique(files$iteration)))

    # Get differences for models
    res <- files
    diff <- res[grep('Diff', res$X),]
    diff <- diff[order(diff$AUC),]
    lengthd <- dim(diff)[1]

    # Add models with no features selected if applicable to weight median calc
    ex <- extra[[i]]
    blank <- diff[1,]
    blank[1,1:8] <- 0
    if(ex > 0){
      for(j in 1:ex){
        diff <- rbind(diff, blank)
      }
    }

    # Select median model as proteinscore
    diff <- diff[order(diff$AUC),]
    length <- dim(diff)[1]
    cut <- as.integer(length/2)
    median <- diff[cut,]
    median$Outcome <- name
    median$available <- length(unique(files$iteration))

    # Get range of differences
    t1 <- res[which(res$iteration %in% median$iteration),]
    min <- min(diff$AUC)
    max <- max(diff$AUC)

    median$min_diff <- min
    median$max_diff <- max

    list_chosen[[i]] <- median
  }, error = function(e) cat("skipped"))
}

res <- do.call(rbind, list_chosen)
res <- res[order(-res$AUC),]
write.csv(res, '/path_to_file.../Run_210723_processing/models_chosen_accounted_features.csv', row.names = F)

# #########################################################################################################

### Collate weights for chosen models

chosen <- read.csv('/path_to_file.../Run_210723_processing/models_chosen_accounted_features.csv')

location_out <- "/path_to_file.../ProteinScores/Run_210723/"

chosen_it <- chosen$iteration
chosen_dis <- chosen$Outcome

res_weights <- list()

traits <- chosen_dis
table <- data.frame(Trait = 1:length(traits), Features = 1:length(traits), Proteins = 1:length(traits))

for(i in 1:length(traits)){
  tryCatch({
    print(traits[i])
    name <- traits[i]
    loc <- chosen_it[i]
    file <- read.csv(paste0(location_out, name, '/files/weights_', loc, '.csv'))
    file$Trait <- traits[i]
    res_weights[[i]] <- file
    dim <- dim(file)[1]
    table[i,1] <- traits[i]
    table[i,2] <- dim
    file$Protein <- gsub("\\..*","",file$Protein)
    prot <- file$Protein
    table[i,3] <- paste(prot, sep = ",", collapse = ", ")
  }, error = function(e) cat("skipped"))
}

bind1 <- do.call(rbind, res_weights)

table1 <- table[order(table$Features),]

library(tidyverse)
library(readxl)
naming <- read_excel("/path_to_file.../naming_index.xlsx")
naming <- as.data.frame(naming)
names(naming)[1] <- 'Trait'
naming <- naming[c(1,3)]

bind1 <- left_join(bind1, naming, by = 'Trait')

table1 <- left_join(table1, naming, by = 'Trait')

write.csv(bind1, '/path_to_file.../Run_210723_processing/features_selected_table_weights.csv', row.names = F)

write.csv(table1, '/path_to_file.../Run_210723_processing/features_summary_table.csv', row.names = F)


################################################v

# Create a plot showing number of features per trait

dat <- read.csv("/path_to_file.../Run_210723_processing/features_summary_table.csv")

# Plot traits by number of features

dat$onset <- 10
dat$onset[1] <- 5
dat$onset[5] <- 5
dat$onset[7] <- 5


library(ggplot2)

pdf("/path_to_file.../Plot_num_features.pdf", width = 9, height = 4)
ggplot(dat, aes(x = reorder(Naming, -Features), y = Features, fill=as.factor(onset))) +
  geom_bar(stat="identity", width = 0.7) +
  # geom_text(nudge_x = 1) +
  scale_fill_manual(values = c('lightblue2',"midnightblue")) +
  # geom_text(aes(label = Count), hjust = -1) +
  xlab("Number of associations with protein levels") +
  ylab("") +
  coord_flip() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  theme(panel.grid.major.y = element_blank()) + theme(legend.title = element_blank()) + theme(legend.position = 'None') +
  labs(y = "Number of features",
       x = "") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12, angle=0), axis.title = element_text(size=12))
dev.off()





















