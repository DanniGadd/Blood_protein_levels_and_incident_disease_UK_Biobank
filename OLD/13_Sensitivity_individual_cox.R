###############################################################################################

##### Sensitivity Cox 

###############################################################################################

args <- commandArgs(trailingOnly=TRUE)
taskid <- as.numeric(args[1])
print(paste0('taskid for this iteration is ', taskid))

library(survival)
library(survminer)
library(tidyverse)

# Set disease and read in disease codes for the iteration (array)
iteration <- taskid
location <- '/path_to_files.../cox_tables/'
files <- list.files(location, '.csv')
diseases <- sub(".csv", "", files)
tables <- list.files('/path_to_files.../cox_tables/', '.csv')
disease <- as.character(diseases[iteration])
table <- read.csv(paste0('/path_to_files.../cox_tables/', disease ,'.csv'))

# Specify sex traits
sex_list <- c('Prostate_FO', 'GYN_FO', 'ENDO_FO', 'Breast_FO', 'CYS_FO')
clock <- names(table)[56:1523] 

print(disease)
print(dim(table))

# Add PA
PA <- read.csv('/path_to_files.../PA.csv')
table <- left_join(table, PA, by = 'SampleID')

# RUN LOOPS OVER YEAR OF FOLLOW UP FOR EACH ASSOC

tot <- 1468 * 17
results <- data.frame(A = 1:tot, B = 1:tot, C = 1:tot, D = 1:tot, E = 1:tot, Z = 1:tot, G = 1:tot, H = 1:tot, I = 1:tot, J = 1:tot)
timer <- seq(0,tot, by = 17)

for(i in 1:length(clock)){
  tryCatch({ 
    cox <- table
    prot <- as.character(clock[[i]])
      for(j in 1:17){
        tryCatch({ 
          k <- timer[[i]]
          
          # cox2 <- cox[which(cox$tte <= j),]

          # Isolate cases and controls
          cases <- cox[which(cox$Event %in% 1),]
          controls <- cox[-which(cox$Event %in% 1),]

          # Restrict to those in follow up time and assign those outside of it as controls
          cases_inc <- cases[which(cases$tte <= j),]
          cases_exc <- cases[which(cases$tte > j),]
          
          # Join to original controls
          cox2 <- rbind(cases_inc, controls)
          
          # If recoded cases present, add in
          if(dim(cases_exc)[1] > 0){
            cases_exc$Event <- 0
            cox2 <- rbind(cox2, cases_exc)
          }

          if(disease %in% sex_list){
            print('Disease in sex list - removing sex from cox.')
            mod <- coxph(Surv(cox2$tte, cox2$Event) ~ scale(cox2[,prot]) + cox2$Age_assessment + cox2$BMI + cox2$Dep + as.factor(cox2$Alc) + as.factor(cox2$Smo) + as.factor(cox2$Edu) + as.factor(cox2$PA), data = cox2)
            results[j+k,1] <- prot
            results[j+k,2] <- disease
            results[j+k,3] <- j
            results[j+k,4:6]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
            results[j+k,7] <- summary(mod)$coefficients[1,5]
            results[j+k,8] <- mod$n[1] - mod$nevent[1]
            results[j+k,9] <- mod$nevent[1]
            table1 <- cox.zph(mod)
            p1 <- table1$table[,"p"]
            results[j+k,10] <-p1[1]
            results[j+k,11] <-p1[9]
          } else {
            mod <- coxph(Surv(cox2$tte, cox2$Event) ~ scale(cox2[,prot]) + factor(cox2$Sex) + cox2$Age_assessment + cox2$BMI + cox2$Dep + as.factor(cox2$Alc) + as.factor(cox2$Smo) + as.factor(cox2$Edu) + as.factor(cox2$PA), data = cox2)
            results[j+k,1] <- prot
            results[j+k,2] <- disease
            results[j+k,3] <- j
            results[j+k,4:6]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
            results[j+k,7] <- summary(mod)$coefficients[1,5]
            results[j+k,8] <- mod$n[1] - mod$nevent[1]
            results[j+k,9] <- mod$nevent[1]
            table1 <- cox.zph(mod)
            p1 <- table1$table[,"p"]
            results[j+k,10] <-p1[1]
            results[j+k,11] <-p1[10]
          }
        }, error = function(e) cat("skipped"))
      }
    print(i)
  }, error = function(e) cat("skipped"))
}

write.csv(results, paste0('/path_to_files.../sensitivity_cox/', disease, '_run_PA_corrected_allocations_controls.csv'), row.names = F)

end <- Sys.time()
print(paste0('Models complete. Time end stamp: ', end, '.'))


####################################################################################

### Collate sensitivity cox results - restricted by both year of follow up

# srun -p interactive --pty bash

library(tidyverse)

location <- '/path_to_files.../sensitivity_cox/'

files <- list.files(location, '_run_PA_corrected_allocations_controls.csv')

res <- list()

for(i in 1:length(files)){
  name <- files[i]
  name <- gsub("\\..*","", name)
  data <- read.csv(paste0(location, name, '.csv'))
  print(name)
  print(dim(data))
  colnames(data) <- c("Predictor", "Outcome", 'Iteration',  "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph_protein", "cox.zph_global")
  data <- as.data.frame(data)
  res[[i]] <- data
}

bind <- res[[1]]

for(i in 2:length(res)){
  bind <- rbind(bind, res[[i]])
  print(i)
  print(unique(res[[i]]$Outcome))
}
# > dim(bind)
# [1] 598944

# # Write off joint file with all associations
write.csv(bind, '/path_to_files.../iteration_sens_joint_assocs_allocation_controls_PA_corrected.csv')

# # ####################################################################################

# Isolate original associations that failed assumptions

library(tidyverse)

bind <- read.csv("/path_to_files.../iteration_sens_joint_assocs_allocation_controls_PA_corrected.csv")
table <- data.frame(Disease = 1:3, Protein = 1:3, Global = 1:3, Both = 1:3)
keep <- read.csv("/path_to_files.../Associations_retained_Bon.csv")
bind$retain <- paste0(bind$Predictor, bind$Outcome)
index <- bind[which(bind$retain %in% keep$retain),]

sixteen <- index[which(index$Iteration %in% 16),]
ten <- index[which(index$Iteration %in% 10),]
five <- index[which(index$Iteration %in% 5),]

sixteen <- sixteen %>% mutate(col = case_when(
  sixteen$cox.zph_protein < 0.05 & sixteen$cox.zph_global > 0.05 ~ 'fail_protein',
  sixteen$cox.zph_global >= 0.05 & sixteen$cox.zph_global >= 0.05 ~ 'not_fail'
))

fail <- sixteen[which(sixteen$cox.zph_protein < 0.05),]

ass1 <- as.data.frame(table(sixteen$Outcome))
tp1 <- as.data.frame(table(fail$Outcome))

# Look to see how many still fail at 10-years
fail <- ten[which(ten$cox.zph_protein < 0.05),]
ten_sig <- ten[which(ten$P.Value < 0.0000031),]
ass2 <- as.data.frame(table(ten_sig$Outcome))
tp2 <- as.data.frame(table(fail$Outcome))

fail <- five[which(five$cox.zph_protein < 0.05),]
five_sig <- five[which(five$P.Value < 0.0000031),]

ass3 <- as.data.frame(table(five_sig$Outcome))
tp3 <- as.data.frame(table(fail$Outcome))

t <- left_join(ass1, tp1, by = 'Var1')
t <- left_join(t, ass2, by = 'Var1')
t <- left_join(t, tp2, by = 'Var1')
t <- left_join(t, ass3, by = 'Var1')
t <- left_join(t, tp3, by = 'Var1')
t <- t[order(t$Freq.x),]
write.csv(t, '/path_to_files.../iteration_sens_summary_table_Bon.csv', row.names = F)
