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
location <- 'Results/Cox/Tables/'
files <- list.files(location, '.csv')
diseases <- sub(".csv", "", files)
tables <- list.files('Results/Cox/Tables/', '.csv')
disease <- as.character(diseases[iteration])
table <- read.csv(paste0('Results/Cox/Tables/', disease ,'.csv'))

# Specify sex traits
sex_list <- c('Prostate_FO', 'GYN_FO', 'ENDO_FO', 'Breast_FO', 'CYS_FO')
clock <- names(table)[56:1523] 
print(disease)
print(dim(table))

# RUN LOOPS OVER YEAR OF FOLLOW UP FOR EACH ASSOC
tot <- 1468 * 16
results <- data.frame(A = 1:tot, B = 1:tot, C = 1:tot, D = 1:tot, E = 1:tot, Z = 1:tot, G = 1:tot, H = 1:tot, I = 1:tot, J = 1:tot)
timer <- seq(0,tot, by = 16)

for(i in 1:length(clock)){
  tryCatch({ 
    cox <- table
    prot <- as.character(clock[[i]])
      for(j in 1:16){
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

write.csv(results, paste0('Sens_successive/', disease, '_run_allocations_controls.csv'), row.names = F)

end <- Sys.time()
print(paste0('Models complete. Time end stamp: ', end, '.'))


####################################################################################

### Collate sensitivity cox results - restricted by both year of follow up

# srun -p interactive --pty bash
library(tidyverse)
location <-
'Results/Cox/Sens_successive/'

files <- list.files(location, '_run_allocations_controls.csv')


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

# Write off joint file with all associations

write.csv(bind,
'iteration_sens_joint_assocs_allocation_controls_PA_corrected.csv')

####################################################################################

# Isolate original associations that failed assumptions

library(tidyverse)

bind <-
read.csv("/iteration_sens_joint_assocs_allocation_controls_PA_corrected.csv")
bind <- bind[-which(bind$Iteration %in% 16),]
table <- data.frame(Disease = 1:3, Protein = 1:3, Global = 1:3, Both = 1:3)
keep <-
read.csv("Cox/Associations_retained_Bon.csv")
bind$retain <- paste0(bind$Predictor, bind$Outcome)
index <- bind[which(bind$retain %in% keep$retain),]

# Now index the associations that failed in the follow-up analyses at 10 years
# and 5 years

sixteen <- index[which(index$Iteration %in% 15),]
ten <- index[which(index$Iteration %in% 10),]
five <- index[which(index$Iteration %in% 5),]
sixteen <- sixteen %>% mutate(col = case_when( sixteen$cox.zph_protein < 0.05
& sixteen$cox.zph_global > 0.05 ~ 'fail_protein', sixteen$cox.zph_global <
0.05 & sixteen$cox.zph_protein > 0.05 ~ 'fail_global', sixteen$cox.zph_protein
< 0.05 & sixteen$cox.zph_global < 0.05~ 'fail_global_and_protein',
sixteen$cox.zph_global >= 0.05 & sixteen$cox.zph_global >= 0.05 ~ 'not_fail'
))

fail <- sixteen[which(sixteen$cox.zph_protein < 0.05),] 
fail2 <- sixteen[which(sixteen$cox.zph_global < 0.05),]
fail4 <- sixteen[which(sixteen$cox.zph_global < 0.05 | sixteen$cox.zph_protein
< 0.05),]

ass1 <- as.data.frame(table(sixteen$Outcome)) 
tp1 <- as.data.frame(table(fail$Outcome)) 
tg1 <- as.data.frame(table(fail2$Outcome))
tgp1 <- as.data.frame(table(fail4$Outcome))

# Look to see how may still fail at 10-years

fail <- ten[which(ten$cox.zph_protein < 0.05),]
fail2 <- ten[which(ten$cox.zph_global < 0.05),]
fail4 <- ten[which(ten$cox.zph_global < 0.05 | ten$cox.zph_protein < 0.05),]
ten_sig <- ten[which(ten$P.Value < 0.0000031),]
ass2 <- as.data.frame(table(ten_sig$Outcome)) 
tp2 <- as.data.frame(table(fail$Outcome)) 
tg2 <- as.data.frame(table(fail2$Outcome))
tgp2 <- as.data.frame(table(fail4$Outcome))

# Look to see how many still fail at 5-years
fail <- five[which(five$cox.zph_protein < 0.05),]  # 409
fail2 <- five[which(five$cox.zph_global < 0.05),]
fail4 <- five[which(five$cox.zph_global < 0.05 | five$cox.zph_protein <
0.05),]
five_sig <- five[which(five$P.Value < 0.0000031),]
ass3 <- as.data.frame(table(five_sig$Outcome)) 
tp3 <- as.data.frame(table(fail$Outcome))
tg3 <- as.data.frame(table(fail2$Outcome))
tgp3 <- as.data.frame(table(fail4$Outcome))

t <- left_join(ass1, tp1, by = 'Var1')
t <- left_join(t, ass2, by = 'Var1') 
t <- left_join(t, tp2, by = 'Var1')
t <- left_join(t, ass3, by = 'Var1') 
t <- left_join(t, tp3, by = 'Var1')

write.csv(t, 'iteration_sens_summary_table_Bon.csv',
row.names = F)

# Subset overall assocs to those kept in original 3209 - 51344
write.csv(index,
'Results/Cox/index_sens_cox.csv',
row.names = F)

# # Annotate naming on main table 
assocs <- read.csv("iteration_sens_joint_assocs_allocation_controls_PA_corrected.csv")
library(tidyverse) 
library(readxl) 
naming <- read_excel("PPP_core_input_files/naming_index.xlsx")
naming <- as.data.frame(naming) 
names(naming)[1] <- 'Outcome' 
naming <- naming[c(1,3)]
assocs <- left_join(assocs, naming, by = 'Outcome')
write.csv(assocs,
"iteration_sens_joint_assocs_allocation_controls_PA_corrected.csv",
row.names = F)

index <- left_join(index, naming, by = 'Outcome')
index <- na.omit(index)
write.csv(index, 'index_6744_assocs_loops_table_Bon.csv', row.names = F)



