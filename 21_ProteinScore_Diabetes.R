###############################################################################################

### Profiling phenotypes in UKB - additional clinical profiling and value exploration

###############################################################################################

# srun -p interactive --pty bash

library(tidyverse)
library(survival)
library(gbm)
library(precrec)
library(ggplot2)

# Load in additional covariates required in rerun
merge <- read.csv('Covariate_preps/imputed_covs_transformed.csv')

# Load chosen model locations
location_out <- "01_paper/Results/Cox/Proteinscores/"
core_models_output <- '01_paper/Results/Cox/Score_processing/covariate_assessment/'

chosen <- read.csv("01_paper/Results/Cox/Score_processing/models_chosen_accounted_features.csv")

location_out <- "00_Run_260823/01_paper/Results/Cox/Proteinscores/"

chosen_it <- chosen$iteration
chosen_dis <- chosen$Outcome
diseases <- chosen_dis

# Load diabetes info
i <- 1
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
threshold <- 10

# Load test data for specific disease iteration
x_test <- readRDS(paste0(location_out, name, '/files/x_test_', it, '.rds'))
y_test <- readRDS(paste0(location_out, name, '/files/y_test_', it, '.rds'))

# Merge in added covariates
y_test <- y_test[-c(2:6)]
y_test <- left_join(y_test, merge, by = 'SampleID')
# print(head(y_test))

# Load model
model_load <- readRDS(paste0(location_out, name, '/files/model_', it, '.rds'))

# Load PGS
PGS <- read.table("/home/dgadd/ukb671220.tab", header = T)
PGS <- PGS[c(1,35)]
names(PGS) <- c('SampleID', 'T2D_PRS')

y_test <- left_join(y_test, PGS, by = 'SampleID')

# Run assessment of hba1c, PRS and ProteinScore 

library(readxl)
library(tidyverse)
library(MethylPipeR)
library(glmnet)
library(survival)
library(gbm)
library(data.table)
library(pacman)

set.seed(seed_iteration)
mprModel <- model_load
testData <- x_test
testTarget <- y_test
seed_iteration <- see
iter <- iteration
name_dis <- name

# Generate model predictions in the test data
mprModelTestPredictions <- predictMPRModel(mprModel,
                                           data = testData,
                                           s = 'lambda.min')

# Add test scores to the target y file in the test sample
testTarget$dScore <- mprModelTestPredictions
testTarget <- na.omit(testTarget)

test_count <- dim(testTarget)[1]
cases_total <- length(which(testTarget$Event == 1))

testTarget <- na.omit(testTarget)

###############################################################################################

# Rank-based inverse normalisation and scaling of each marker

untransformed <- testTarget

testTarget$hba1c <- qnorm((rank(testTarget$hba1c, na.last='keep')-0.5)/sum(!is.na(testTarget$hba1c)))

testTarget$hba1c <- scale(testTarget$hba1c)

mean(testTarget$hba1c, na.rm = T)
sd(testTarget$hba1c, na.rm = T)

testTarget$dScore <- qnorm((rank(testTarget$dScore, na.last='keep')-0.5)/sum(!is.na(testTarget$dScore)))

testTarget$dScore <- scale(testTarget$dScore)

mean(testTarget$dScore, na.rm = T)
sd(testTarget$dScore, na.rm = T)

transformed <- testTarget



# riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex), testTarget)
riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex), testTarget)
lifestyleCoxPH <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                          as.factor(Smo) + as.factor(Alc), testTarget)
PGSonly <- coxph(Surv(time_to_event, Event) ~  T2D_PRS, testTarget)
ProteinOnly <- coxph(Surv(time_to_event, Event) ~ dScore, testTarget)
hbOnly <- coxph(Surv(time_to_event, Event) ~ hba1c, testTarget)
hbprot <- coxph(Surv(time_to_event, Event) ~  hba1c + dScore, testTarget)
hbprotpgs <- coxph(Surv(time_to_event, Event) ~ hba1c + dScore + T2D_PRS, testTarget)
exset <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(Smo) + as.factor(Alc) + as.factor(PA) +
                 sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                 albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                 creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L, testTarget)
exsetprot <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                     as.factor(Smo) + as.factor(Alc) + 
                     sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                     albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                     creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L + dScore, testTarget)
exsetprotPRS <- coxph(Surv(time_to_event, Event) ~ Age_assessment + as.factor(Sex) + as.factor(Ed) + BMI + Dep + as.factor(PA) +
                        as.factor(Smo) + as.factor(Alc) + 
                        sys_bp_mmHg + leukocyte + erythrocyte + Hb_g_dec + mean_corp_vol + platelet_count +
                        albumin + alanine_aminotransferase + aspartate.aminotransferase + urea + cholesterol_mmol_L +
                        creatinine_umol_L + CRP + cystatin_C + glucose + hba1c + LDL_mmol_L + triglycerides_mmol_L + dScore + T2D_PRS, testTarget)
# print('Model has not been sex stratified in cox.')

# List models
models <- list(p = riskFactorsOnlyCoxPH,
               l = lifestyleCoxPH,
               o = PGSonly, 
               t = ProteinOnly,
               j = hbOnly, 
               q = hbprot, 
               z = hbprotpgs,
               a = exset,
               b = exsetprot,
               c = exsetprotPRS)


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

# Calculate p values for AUC comparisons of interest (age/sex + protein)
p_second <- pROC::roc.test(response = testTarget$Event, predictor1 = predictions$t, predictor2 = predictions$q)$p.value

p_third <- pROC::roc.test(response = testTarget$Event, predictor1 = predictions$q, predictor2 = predictions$z)$p.value


# Extract metrics table
aucs <- sapply(testResults, function(r) {r$auc})
praucs <- sapply(testResults, function(r) {r$prauc})
metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)

row.names(metricsTable) <- c('Minimally-adjusted',
                             'Lifestyle-adjusted',
                             'PRS',
                             'ProteinScore',
                             'HbA1c',
                             'ProteinScore + HbA1c',
                             'ProteinScore +HbA1c + PRS',
                             'Extended set',
                             'Extended set + ProteinScore',
                             'Extended set + ProteinScore + PRS')

metricsTable <- metricsTable[order(metricsTable$AUC),]


metricsTable[11,1] <- metricsTable[8,1] - metricsTable[4,1]
metricsTable[11,2] <- metricsTable[8,2] - metricsTable[4,2]
metricsTable[12,1] <- metricsTable[10,1] - metricsTable[8,1]
metricsTable[12,2] <- metricsTable[10,2] - metricsTable[8,2]


# Add P values for AUCs to table
metricsTable$P_ROCs <- 'NA'
metricsTable[11,3] <- p_second
metricsTable[12,3] <- p_third

write.csv(metricsTable, '/Results/Cox/metab_analysis/diabetes_comparison_ROC_p_PGS.csv', row.names = T)


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
rocs_df[which(rocs_df$Model =='p'), "Model"] = "Minimally-adjusted"
rocs_df[which(rocs_df$Model =='l'), "Model"] = "Lifestyle-adjusted"
rocs_df[which(rocs_df$Model =='o'), "Model"] = "PRS"
rocs_df[which(rocs_df$Model =='t'), "Model"] = "ProteinScore"
rocs_df[which(rocs_df$Model =='j'), "Model"] = "HbA1c"
rocs_df[which(rocs_df$Model =='q'), "Model"] = "ProteinScore + HbA1c"
rocs_df[which(rocs_df$Model =='z'), "Model"] = "ProteinScore + HbA1c + PRS"
rocs_df[which(rocs_df$Model =='a'), "Model"] = "Extended set"
rocs_df[which(rocs_df$Model =='b'), "Model"] = "Extended set + ProteinScore"
rocs_df[which(rocs_df$Model =='c'), "Model"] = "Extended set + ProteinScore + PRS"


My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size =12),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 12, face = "bold"),
  legend.text=element_text(size=16),
  legend.title=element_text(size=16, face = "bold"))


rocs_df <- rocs_df[which(rocs_df$Model %in% c("ProteinScore + HbA1c + PRS",
                                              "ProteinScore + HbA1c",
                                              "HbA1c",
                                              "ProteinScore",
                                              "PRS",
                                              "Lifestyle-adjusted",
                                              "Minimally-adjusted")),]


# > unique(rocs_df$Model)
# [1] "Minimally-adjusted"                "Lifestyle-adjusted"
# [3] "PRS"                               "ProteinScore"
# [5] "HbA1c"                             "ProteinScore + HbA1c"
# [7] "ProteinScore + HbA1c + PRS"        "Extended set"
# [9] "Extended set + ProteinScore"       "Extended set + ProteinScore + PRS"


rocs_df$Model <-  factor(rocs_df$Model, levels = c(
                                                   "ProteinScore + HbA1c + PRS",
                                                   "ProteinScore + HbA1c",
                                                   
                                                   "ProteinScore",
                                                   "HbA1c",
                                                   "Lifestyle-adjusted",
                                                   "PRS",
                                                   "Minimally-adjusted"))

cbPalette <- c("red",   'yellow', "tan1", "darkturquoise", 'grey', "blue",
                  "black")

# "gold1", "violetred1", "darkturquoise",
# draw a couple of ROC curves
pdf('Results/Cox/metab_analysis/DIAB_rocs_diabetes_hba1c_PGS.pdf', width = 11, height = 6)
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

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size =12),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 12, face = "bold"),
  legend.text=element_blank(),
  legend.title=element_blank())

pdf('Results/Cox/metab_analysis/DIAB_rocs_diabetes_hba1c_PGS_unlabelled.pdf', width = 9, height = 6)
rocs_df %>%
  ggplot( aes(x=1-Specificity, y=Sensitivity, group=Model, color=Model)) +
  geom_line() +
  theme_light() +
  # scale_color_viridis(discrete=TRUE) 
  scale_colour_manual(values=cbPalette)   +
  theme(legend.position = "none", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 8)) +
  xlab("False Positive Rate (1 - Specificity)") +
  ylab("True Postive Rate (Sensitivity)")  + theme_classic() + My_Theme
dev.off()


###############################################################################################

# Plot with scaling 
diab <- transformed
# diab <- testTarget
diab <- diab[order(diab$dScore),]
diab <- diab %>% mutate(decile = ntile(diab$dScore, 10))
diab$decile <- as.factor(diab$decile)


diab$Event <- as.factor(diab$Event)
library(ggpubr)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size =12),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 12, face = "bold"),
  legend.text=element_text(size=16),
  legend.title=element_text(size=16, face = "bold"), legend.position = "none")


pdf('plot_scatter_hba1c_diab_no_legend_loes.pdf', width = 6, height = 5)
ggplot(diab, aes(x=dScore, y=hba1c, group = Event)) +
  geom_point(size = 0.2)+
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Event, bins = 3)) + theme_minimal() + theme_classic() +
  scale_fill_manual(values = c("dodgerblue", "tomato")) +
  geom_density_2d(aes(color = Event)) +
  scale_color_manual(values = c("blue4", "firebrick1")) + My_Theme + xlab('Type 2 diabetes ProteinScore') + ylab('HbA1c')
dev.off()

# geom_smooth(method = "loess") +

###############################################################################################

# Plot without scaling 

diab <- untransformed

# raw hba1c and join
phen <- read.table(gzfile("all_quant_20200904.pheno.v2.tsv.gz"), header = T)
hb <- phen[which(colnames(phen) %in% c('IND', 'f_30750_0_0'))]
names(hb)[1] <- 'SampleID'

diab <- left_join(diab, hb, by = 'SampleID')

diab <- diab[order(diab$dScore),]
diab <- diab %>% mutate(decile = ntile(diab$dScore, 10))
diab$decile <- as.factor(diab$decile)

diab$Event <- as.factor(diab$Event)
library(ggpubr)

# Do an anova one way on deciles to get a p value for inclusion
# https://www.r-bloggers.com/2022/05/one-way-anova-example-in-r-quick-guide/
res.aov <- aov(f_30750_0_0 ~ decile, data = diab)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size =12),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 12, face = "bold"),
  legend.text=element_text(size=16),
  legend.title=element_text(size=16, face = "bold"), legend.position = "none")



pdf('plot_deciles_rect_hba1c_diab_no_scaling.pdf', width = 6, height = 5)
ggplot(diab, aes(x=decile, y=f_30750_0_0)) + 
  geom_violin(trim=FALSE, fill='lightcyan2', color="cadetblue4")+
  geom_boxplot(width=0.1, color="darkslategrey", fill="darkslategrey", alpha=0.2) + theme_classic() + 
  labs(x="Type 2 diabetes ProteinScore decile", y = 'HbA1c (mmol/mol)')+ 
  annotate("rect", ymin = 42, ymax = 47, xmin = 0, xmax = 10.5,
           alpha = .1,fill = "blue") +
  My_Theme  + ylim(0,130)
dev.off()


# Get stats for test set 
ev <- testTarget[which(testTarget$Event %in% 1),]
mean(ev$time_to_event, na.rm = T)
sd(ev$time_to_event, na.rm = T)

# Write out for metabolomics test set analyses
write.csv(testTarget, 'Censor_test/diab_test_set_PRS.csv', row.names = F)
