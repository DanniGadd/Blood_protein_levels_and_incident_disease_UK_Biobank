###############################################################################################

## Test scores - with covariate combination models 

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

### Define diseases from chosen models

location_out <- "/path_to_file.../ProteinScores/Run_210723/"
core_models_output <- '/path_to_file.../ProteinScores/Run_210723_processing/covariate_assessment/'

chosen <- read.csv("/path_to_file.../ProteinScores/Run_210723_processing/models_chosen_accounted_features.csv")

location_out <- "/path_to_file.../ProteinScores/Run_210723/"

chosen_it <- chosen$iteration
chosen_dis <- chosen$Outcome
diseases <- chosen_dis


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
      allrisk <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Ed) + BMI + Dep + as.factor(Smo) + as.factor(Alc) + as.factor(PA), testTarget)
      allriskprotein <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Ed) + BMI + Dep + as.factor(Smo) + as.factor(Alc)+ as.factor(PA) + dScore, testTarget)
      
      
      clinicalonly <- coxph(Surv(time_to_event, Event) ~ Age_assessment + sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                              albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                              creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L, testTarget)
      clinicalrisk <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Ed) + BMI + Dep + as.factor(Smo) + as.factor(Alc) + as.factor(PA) +
                              sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                              albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                              creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L, testTarget)
      clinicalriskprotein <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                                     as.factor(Smo) + as.factor(Alc) + 
                                     sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                                     albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                                     creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L + dScore, testTarget)
      

      
      print('Model has been sex stratified in cox.')
      print('Cases and controls available:')
      print(table(testTarget$Event))
    } else {
      riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex), testTarget)
      fullCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + dScore, testTarget)
      ProteinOnly <- coxph(Surv(time_to_event, Event) ~ dScore, testTarget)
      allrisk <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(Smo) + as.factor(Alc) + as.factor(PA), testTarget)
      allriskprotein <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(Smo) + as.factor(Alc) + as.factor(PA) + dScore, testTarget)
      
      clinicalonly <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                              albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                              creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L, testTarget)
      clinicalrisk <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(Smo) + as.factor(Alc) + as.factor(PA) +
                              sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                              albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                              creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L, testTarget)
      clinicalriskprotein <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                                     as.factor(Smo) + as.factor(Alc) + 
                                     sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                                     albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                                     creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L + dScore, testTarget)
      
      print('Cases and controls available:')
      print(table(testTarget$Event))
    }
    
    # List models
    models <- list(r = riskFactorsOnlyCoxPH, f = fullCoxPH, t = ProteinOnly, s = allrisk, b = allriskprotein,
                   x = clinicalonly , y = clinicalrisk, z = clinicalriskprotein)

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
    
    # Extract onset predictions
    predictions <- lapply(testResults, function(r) {r$onsetPredictions})
    
    # Calculate p values 
    p_null <- pROC::roc.test(response = testTarget$Event, predictor1 = predictions$r, predictor2 = predictions$f)$p.value

    p_second <- pROC::roc.test(response = testTarget$Event, predictor1 = predictions$s, predictor2 = predictions$b)$p.value
    
    p_third <- pROC::roc.test(response = testTarget$Event, predictor1 = predictions$y, predictor2 = predictions$z)$p.value
    
    # Extract metrics table
    aucs <- sapply(testResults, function(r) {r$auc})
    praucs <- sapply(testResults, function(r) {r$prauc})
    metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)
    metricsTable[9,1] <- metricsTable[2,1] - metricsTable[1,1]
    metricsTable[9,2] <- metricsTable[2,2] - metricsTable[1,2]
    metricsTable[10,1] <- metricsTable[5,1] - metricsTable[4,1]
    metricsTable[10,2] <- metricsTable[5,2] - metricsTable[4,2]
    metricsTable[11,1] <- metricsTable[8,1] - metricsTable[6,1]
    metricsTable[11,2] <- metricsTable[8,2] - metricsTable[6,2]
    row.names(metricsTable) <- c('AgeSex', 'AgeSex_ProteinScore', 'ProteinScore Only', 'AgeSexCovs', 'AgeSexCovs_ProteinScore',
                                 'ClinicalOnly', 'Clinicalrisk', 'Clinicalriskprotein',
                                 'Diff_AgeSex', 'Diff_AgeSexCovs', 'Diff_clinical_protein')
    metricsTable$Outcome <- trait
    metricsTable$seed <- seed_iteration
    metricsTable$iteration <- iter
    
    # Add P values for AUCs to table
    metricsTable$P_ROCs <- 'NA'
    metricsTable[9,6] <- p_null
    metricsTable[10,6] <- p_second
    metricsTable[11,6] <- p_third
    
    metricsTable$test_total <- test_count
    metricsTable$test_cases <- cases_total
    metricsTable$test_cases_recoded <- cases_recoded

    saveRDS(testTarget, paste0(output, name_dis, '/files/testResults_', iter, '.rds'))

  }, error = function(e) cat("skipped"))
  return(metricsTable)
}


###############################################################################################

### Run testing

list <- list()

set <- 1:length(diseases)

for(i in set){
  tryCatch({
    name <- as.character(diseases[i])
    trait <- as.character(diseases[i])
    print(i)
    print(trait)
    print(name)
    
    # Identify median model iteration and batch
    chosen_sub <- chosen[which(chosen$Outcome %in% trait),]
    it <- chosen_sub$iteration
    bat <- chosen_sub$Batch
    see <- chosen_sub$seed
    set.seed(see)
    seed <- see
    
    iteration <- it
    
    
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
    
    # Load test data for specific disease iteration
    x_test <- readRDS(paste0(location_out, name, '/files/x_test_', it, '.rds'))
    y_test <- readRDS(paste0(location_out, name, '/files/y_test_', it, '.rds'))
    
    # Merge in added covariates
    y_test <- y_test[-c(2:6)]
    y_test <- left_join(y_test, merge, by = 'SampleID')
    model_load <- readRDS(paste0(location_out, name, '/files/model_', it, '.rds'))
    
    # Test models
    mod <- test(model_load, x_test, y_test, list, seed, iteration, name, location_out, thr_input)
    
    # Record results
    mod$Outcome <- name
    mod$Dis_num <- i
    mod$Onset_threshold <- thr_input
    list[[i]] <- mod
    print(mod)
  }, error = function(e) cat("skipped"))
}

result <- do.call(rbind, list)
write.csv(result, paste0(core_models_output, '/metricsTables_chosen_joint.csv'))

