###############################################################################################

##### Run mortality basic and full models

###############################################################################################

# args <- commandArgs(trailingOnly=TRUE)
# taskid <- as.numeric(args[1])
# print(paste0('taskid for this iteration is ', taskid))

taskid <- 1

library(survival)
library(survminer)
library(tidyverse)

print('Packages loaded - now loading d1 file.')

# Read in phenotypes
d1 <- readRDS("/path_to_files.../d1_mortality_202112.rds")

# Add PA as covariate 
PA <- read.csv('/path_to_files.../covariates_additional.csv')
PA <- PA[which(colnames(PA) %in% c('SampleID', 'PA'))]
d1 <- left_join(d1, PA, by = 'SampleID')

# Set protein clock for loops
clock <- names(d1)[56:1523] 

# Set disease and read in disease codes for the iteration (array)
iteration <- taskid

diseases <- c('DEATH')

name <- as.character(diseases[iteration])
print(paste0('disease for this iteration is ', name))
codes <- read.csv(paste0("/path_to_files.../prepped_death/DEATH.csv"))


location_tables <- '/path_to_files.../Death/Tables/'
location_results <- '/path_to_files.../Death/Results/'

now <- Sys.time()
print(paste0('Models initiating. Time start stamp: ', now, '.'))

print('Commencing runs.')

# Basic models per protein

d1_data <- d1

mat_hazard <- matrix(nrow=length(clock),ncol=10)
output_hazard <- as.data.frame(mat_hazard)
for(j in 1:length(clock)){ 
  tryCatch({ 
    dat1= d1_data
    tmp1 = codes[which(codes$SampleID %in% dat1$SampleID),] 
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$SampleID %in% tmp1$SampleID),] 
    age_onset = codes[,c("first", "SampleID")]
    affected = merge(age_onset, affected, by= "SampleID")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1,4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$MOB))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$YOB)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$SampleID %in% codes$SampleID),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$SampleID.y <- NULL
    names(affected)[names(affected)=="SampleID"] <- "SampleID"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age_assessment
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < 0, "NA", cox$tte)
    # cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
      mod = coxph(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Sex) + cox$Age_assessment, data = cox)
      print(j)
      output_hazard[j,1] <- as.character(clock[[j]])
      output_hazard[j,2] <- as.character(name)
      output_hazard[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
      output_hazard[j,6] <- summary(mod)$coefficients[1,5]
      output_hazard[j,8] <- mod$n[1] - mod$nevent[1]
      output_hazard[j,7] <- mod$nevent[1]
      table <- cox.zph(mod)
      p1 <- table$table[,"p"]
      output_hazard[j,9] <-p1[1]
      output_hazard[j,10] <-p1[4]
      cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event)
      # write.csv(cox, paste0(location_tables, name ,'.csv'), row.names = F)

  }, error = function(e) cat("skipped"))
} 

# Save results
comb <- output_hazard
names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph_protein", "cox.zph_global")
comb <- na.omit(comb)
comb <- comb[order(comb$P.Value),]
write.csv(comb, paste0(location_results, name, '.csv'), row.names = F)

end <- Sys.time()
print(paste0('Models complete. Time end stamp: ', end, '.'))



### FULL MODEL DEATH  

d1_data <- d1

mat_hazard <- matrix(nrow=length(clock),ncol=10)
output_hazard <- as.data.frame(mat_hazard)

for(j in 1:length(clock)){ 
  tryCatch({ 
    dat1= d1_data
    tmp1 = codes[which(codes$SampleID %in% dat1$SampleID),] 
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$SampleID %in% tmp1$SampleID),] 
    age_onset = codes[,c("first", "SampleID")]
    affected = merge(age_onset, affected, by= "SampleID")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1,4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$MOB))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$YOB)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$SampleID %in% codes$SampleID),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$SampleID.y <- NULL
    names(affected)[names(affected)=="SampleID"] <- "SampleID"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age_assessment
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < 0, "NA", cox$tte)
    # cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
    mod = coxph(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Sex) + cox$Age_assessment + cox$BMI + cox$Dep + as.factor(cox$Alc) + as.factor(cox$Smo) + as.factor(cox$Edu) + as.factor(cox$PA), data = cox)
    print(j)
    output_hazard[j,1] <- as.character(clock[[j]])
    output_hazard[j,2] <- as.character(name)
    output_hazard[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard[j,6] <- summary(mod)$coefficients[1,5]
    output_hazard[j,8] <- mod$n[1] - mod$nevent[1]
    output_hazard[j,7] <- mod$nevent[1]
    table <- cox.zph(mod)
    p1 <- table$table[,"p"]
    output_hazard[j,9] <-p1[1]
    output_hazard[j,10] <-p1[10]
    cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event)
    # write.csv(cox, paste0(location_tables, name ,'_FULL.csv'), row.names = F)
    
  }, error = function(e) cat("skipped"))
} 

# Save results
comb <- output_hazard
names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph_protein", "cox.zph_global")
comb <- na.omit(comb)
comb <- comb[order(comb$P.Value),]
write.csv(comb, paste0(location_results, name, '_FULL.csv'), row.names = F)

end <- Sys.time()
print(paste0('Models complete. Time end stamp: ', end, '.'))




