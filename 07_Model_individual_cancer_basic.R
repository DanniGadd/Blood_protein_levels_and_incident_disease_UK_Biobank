###############################################################################################

##### Run basic Cox models

###############################################################################################

args <- commandArgs(trailingOnly=TRUE)
taskid <- as.numeric(args[1])
print(paste0('taskid for this iteration is ', taskid))

library(survival)
library(survminer)

print('Packages loaded - now loading d1 file.')

# Set locations
location_codes <- '/path_to_files.../prepped_traits_cancer_reg/'
location_self <- '/path_to_files.../self_report_preps/'

location_tables <- '/path_to_files.../Cancers/Tables/'
location_results <- '/path_to_files.../Cancers/Results/'

# Read in phenotypes
d1 <- readRDS("/path_to_files.../d1_cancer_reg_202012.rds")
# d1 <- as.data.frame(d1)
# Set protein clock for loops
clock <- names(d1)[56:1523] 

# Set disease and read in disease codes for the iteration (array)
iteration <- taskid

diseases <- c('BRAIN_FO', 'Breast_FO', 'Colorectal_FO', 'GYN_FO', 'LUNG_FO', 'Prostate_FO')


restricted_list <- c('AL', 'VD')
sex_list <- c('CYS', 'ENDO', 'Breast_FO', 'GYN_FO')
male_list <- c('Prostate_FO')

name <- as.character(diseases[iteration])
print(paste0('disease for this iteration is ', name))

codes <- read.csv(paste0(location_codes, name, '.csv'))

self <- read.csv(paste0(location_self, name, '.csv'))

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
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline 
    if(dim(self)[1] < 2){
      print('Self less than 2 - skipping')
    } else {
      tmp1 = tmp1[-which(tmp1$SampleID %in% self$SampleID),]
    }
    
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
    
    if(name %in% restricted_list){
      print('Trait in restricted list - subsetting to >=65 years.')
      cox = cox[cox$age_at_event >=65,]
      
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
    } else {
      if(name %in% sex_list){
        print('Trait in sex list - running model without sex as covariate in females.')
        cox <- cox[which(cox$Sex %in% 0),]
        mod = coxph(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age_assessment, data = cox)
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
        output_hazard[j,10] <-p1[3]
        cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event)
        # write.csv(cox, paste0(location_tables, name ,'.csv'), row.names = F)
      } else {
        if(name %in% male_list){
          print('Trait in male list - running model without sex as covariate in males.')
          cox <- cox[which(cox$Sex %in% 1),]
          mod = coxph(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age_assessment, data = cox)
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
          output_hazard[j,10] <-p1[3]
          cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event)
          # write.csv(cox, paste0(location_tables, name ,'.csv'), row.names = F)
        } else {
            print('Trait in neither list - running as standard.')  
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
        }
        
      }
    }
  }, error = function(e) cat("skipped"))
} 

print(dim(cox))

# Save results
comb <- output_hazard
names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph_protein", "cox.zph_global")
comb <- na.omit(comb)
comb <- comb[order(comb$P.Value),]
write.csv(comb, paste0(location_results, name, '.csv'), row.names = F)

end <- Sys.time()
print(paste0('Models complete. Time end stamp: ', end, '.'))



