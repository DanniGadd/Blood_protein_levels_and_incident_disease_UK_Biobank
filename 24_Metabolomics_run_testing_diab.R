###############################################################################################

## Metabolomics UKB - test scores

###############################################################################################

# ### Load sbatch information when run as array per disease
args <- commandArgs(trailingOnly=TRUE)
taskid <- as.numeric(args[1])
print(paste0('taskid for this iteration is ', taskid))

### Interactive session and R required for script troubleshooting
# srun -p interactive --pty bash
# module load R/4.2.0-foss-2021b

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

### Define diseases
i <- 10

location_out1 <- 'Diab/prot/'
location_out2 <- 'Diab/metab/'
location_out3 <- 'Diab/joint/'

# Get disease names
location <- 'Cox/tables_cut/'
files <- list.files(location, '.csv')
diseases <- sub(".csv", "", files)

### Set the disease for this array 
name <- as.character(diseases[[i]])
print(paste0('The disease trait for this iteration is ', name))

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
    if(name_dis %in% c('Prostate', 'CYS', 'GYN', 'Breast', 'ENDO')){
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
    
    # Generate performance statistics - for years of follow up, comparing base HR to survival HR at threshold
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
    saveRDS(testTarget, paste0(output, '/files/testResults_', iter, '.rds'))
    
    # save metrics table performances
    saveRDS(metricsTable, paste0(output, '/files/metricsTable_', iter, '.rds'))
    
  }, error = function(e) cat("skipped"))
  return(metricsTable)
}


###############################################################################################

### Define threshold per trait 

thr_input <- 10 
print(paste0('Threshold set to ', thr_input, ' for test sampling.'))

###############################################################################################

### Run testing proteins

    seed <- 3104
    iteration <- 18
    i <- 18
    it <- 18
    
    x_test <- readRDS(paste0(location_out1, '/files/x_test_', i, '.rds'))
    y_test <- readRDS(paste0(location_out1, '/files/y_test_', i, '.rds'))
    model_load <- readRDS(paste0(location_out1, '/files/model_', it, '.rds'))
    
    
    
    mod <- test(model_load, x_test, y_test, list, seed, iteration, name, location_out1, thr_input)

result <- mod
write.csv(result, paste0(location_out1, '/metricsTables_joint.csv'))


###############################################################################################


### Run testing metab

seed <- 3104
iteration <- 18
i <- 18
it <- 18

x_test <- readRDS(paste0(location_out2, '/files/x_test_', i, '.rds'))
y_test <- readRDS(paste0(location_out2, '/files/y_test_', i, '.rds'))
model_load <- readRDS(paste0(location_out2, '/files/model_', it, '.rds'))

mod <- test(model_load, x_test, y_test, list, seed, iteration, name, location_out2, thr_input)

result <- mod
write.csv(result, paste0(location_out2, '/metricsTables_joint.csv'))


###############################################################################################

### Run testing joint

seed <- 3104
iteration <- 18
i <- 18
it <- 18

x_test <- readRDS(paste0(location_out3, '/files/x_test_', i, '.rds'))
y_test <- readRDS(paste0(location_out3, '/files/y_test_', i, '.rds'))
model_load <- readRDS(paste0(location_out3, '/files/model_', it, '.rds'))

mod <- test(model_load, x_test, y_test, list, seed, iteration, name, location_out3, thr_input)

result <- mod
write.csv(result, paste0(location_out3, '/metricsTables_joint.csv'))


###############################################################################################

# Load the three scores
prot <- readRDS("Diab/prot/files/testResults_18.rds")
metab <- readRDS("Diab/metab/files/testResults_18.rds")
joint <- readRDS("Diab/joint/files/testResults_18.rds")

# Extract the scores and join for same y_test evaluation set
y_test <- prot[-c(2:6)]
y_test <- y_test[,-12]
y_test$prot <- prot$dScore
y_test$metab <- metab$dScore
y_test$joint <- joint$dScore

# Load in additional covariates required in rerun
merge <- read.csv('imputed_covs_transformed.csv')
y_test <- left_join(y_test, merge, by = 'SampleID')
y_sub <- na.omit(y_test)

case <- y_sub[which(y_sub$Event %in% 1),]
mean(case$time_to_event, na.rm = T)
sd(case$time_to_event, na.rm = T)

library(tidyverse)
library(survival)
library(gbm)
library(precrec)
library(ggplot2)

testTarget <- y_sub
riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex), testTarget)
lifestyleCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                          as.factor(Smo) + as.factor(Alc), testTarget)
Metonly <- coxph(Surv(time_to_event, Event) ~  metab, testTarget)
ProtOnly <- coxph(Surv(time_to_event, Event) ~ prot, testTarget)
jointOnly <- coxph(Surv(time_to_event, Event) ~ joint, testTarget)
manualjoint <- coxph(Surv(time_to_event, Event) ~  prot + metab, testTarget)

clinicalrisk <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                               as.factor(Smo) + as.factor(Alc) +
                               sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                               albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                               creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L, testTarget)
clinicalriskjoint <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                               as.factor(Smo) + as.factor(Alc) +
                               sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                               albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                               creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L + prot + metab, testTarget)

# List models
models <- list(z = riskFactorsOnlyCoxPH,
               k = lifestyleCoxPH,
               o = Metonly,
               t = ProtOnly,
               j = jointOnly,
               q = manualjoint,
               a = clinicalrisk,
               b = clinicalriskjoint)

# Generate performance statistics - for 10 years of follow up, comparing base HR to survival HR at 10 years
predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
  uniqueTimes <- sort(unique(dataDF$time_to_event))
  thresholdIndex <- match(threshold, uniqueTimes)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$time_to_event, dataDF$Event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
  onsetPredictions <- 1 - survivalPredictions

  # Event should be 0 if tte is > 10
  dataDF$Event <- sapply(1:nrow(dataDF), function(i) {
    if (dataDF$time_to_event[[i]] > threshold) {
      0
    } else {
      dataDF$Event[[i]]
    }
  })

  auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  roc <- pROC::roc(response = dataDF$Event, predictor = onsetPredictions)
  list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
}

# Extract test results
testResults <- lapply(models, function(m) {predictCoxPHOnset(testTarget, m)})

# Extract onset predictions
predictions <- lapply(testResults, function(r) {r$onsetPredictions})

# Extract metrics table
aucs <- sapply(testResults, function(r) {r$auc})
praucs <- sapply(testResults, function(r) {r$prauc})
metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)

row.names(metricsTable) <- c('AgeSex',
                             'AgeSexLifestyle',
                             'MetabScore',
                             'ProteinScore',
                             'MetabProteinScore',
                             'MetabScore + ProteinScore',
                             'Extended set (24 covariates)',
                             'Extended set (24 covariates) + MetabProteinScore')

metricsTable <- metricsTable[order(metricsTable$AUC),]
write.csv(metricsTable, 'Cox/metab_analysis/metab_comp_diab.csv', row.names = T)

# Assign target test file to variable and subset such that no event greater than 10 years is present (event = 0, if tte > 10)
w1Target <- testTarget

w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
  if (w1Target$time_to_event[[i]] > 10) {
    0
  } else {
    w1Target$Event[[i]]
  }
})

w1Target$Event <- as.factor(w1Target$Event)

# Assign cox test results from calculation above
coxTestResults <- testResults


# quickly convert list to a dataframe
listToDataframe = function(list) {
  models = colnames(list)
  data = data.frame("Sensitivity" =  numeric(), "Specificity" = numeric(), "Model" = character())

  for (model in models) {
    tmp_data = data.frame("Sensitivity" =  list[['sensitivities', model]], "Specificity" = list[['specificities', model]], "Model" = model)
    data = rbind(data, tmp_data)
  }

  return(data)
}

rocs <- sapply(testResults, function(r) {r$roc})
rocs_df = listToDataframe(rocs)
rocs_df[which(rocs_df$Model =='z'), "Model"] = "Age+Sex"
rocs_df[which(rocs_df$Model =='k'), "Model"] = "Age+Sex+Lifestyle"
rocs_df[which(rocs_df$Model =='o'), "Model"] = "MetabScore"
rocs_df[which(rocs_df$Model =='t'), "Model"] = "ProteinScore"
rocs_df[which(rocs_df$Model =='j'), "Model"] = "MetabProteinScore"
rocs_df[which(rocs_df$Model =='q'), "Model"] = "ProteinScore+MetabScore"
rocs_df[which(rocs_df$Model =='a'), "Model"] = "Age+Sex+Lifestyle+Extended set"
rocs_df[which(rocs_df$Model =='b'), "Model"] = "Age+Sex+Lifestyle+Extended set+MetabProteinScore"

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size =12),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 12, face = "bold"),
  legend.text=element_text(size=16),
  legend.title=element_text(size=16, face = "bold"))
# , legend.position = "none"

rocs_df$Model <-  factor(rocs_df$Model, levels = c(
  "Age+Sex+Lifestyle+Extended set+MetabProteinScore",
  "ProteinScore+MetabScore",
  "Age+Sex+Lifestyle+Extended set",
  'MetabProteinScore',
  "ProteinScore",
  "MetabScore",
  'Age+Sex+Lifestyle',
  'Age+Sex'
  ))

cbPalette <- c("red", "darkorchid", "blue", "tan1", 'turquoise', 'pink', 'grey',
               "black")

# "gold1", "violetred1", "darkturquoise",
# draw a couple of ROC curves
pdf('metab_rocs_diab.pdf', width = 13.5, height = 6)
rocs_df %>%
  ggplot( aes(x=1-Specificity, y=Sensitivity, group=Model, color=Model)) +
  geom_line() +
  theme_light() +
  # scale_color_viridis(discrete=TRUE)
  scale_colour_manual(values=cbPalette)   +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 8)) +
  xlab("False Positive Rate (1 - Specificity)") +
  ylab("True Postive Rate (Sensitivity)")  + theme_classic() + My_Theme
dev.off()


# Collate weights table to suppl

metab <- read.csv("Diab/metab/files/weights_18.csv")
prot <- read.csv("Diab/prot/files/weights_18.csv")
joint <- read.csv("Diab/joint/files/weights_18.csv")
metab$Score <- 'MetabScore (metabolomics only)'
prot$Score <- 'ProteinScore (proteomics only)'
joint$Score <- 'MetabProteinScore (metabolomics and proteomics)'

join <- rbind(joint, metab)
join <- rbind(join, prot)
write.csv(join, 'metab_weights_diab.csv', row.names = F)
