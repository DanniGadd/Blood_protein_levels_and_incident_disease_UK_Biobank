###############################################################################################

##### Process models - incident disease and mortality together 

###############################################################################################

# srun -p interactive --pty bash

library(tidyverse)
library(readxl)
naming <- read_excel("PPP_core_input_files/naming_index.xlsx")
naming <- as.data.frame(naming)
names(naming)[1] <- 'Outcome'
naming <- naming[c(1,3)]

# COLLATE RESULTS 
# Set disease and read in disease codes for the iteration (array)
diseases <- c('AL_FO', 'ALS_FO', 'BRAIN_FO', 'Breast_FO', 'Colorectal_FO', 'COPD_FO', 'CYS_FO', 'DEP_FO', 'Diab_FO', 'ENDO_FO', 'GYN_FO',
              'IBD_FO', 'IHD_FO', 'LIV_FO', 'LUNG_FO', 'LUP_FO', 'MS_FO', 'PD_FO', 'Prostate_FO', 'RA_FO', 'SCZ_FO', 'ST_FO', 'VD_FO', 'DEATH_FO')
location <- 'Results/Cox/Results/'
files <- list.files(location, '_FO.csv')

# Basic results collation
results <- list()
for(i in 1:length(diseases)){
  name <- as.character(diseases[i])
  print(name)
  data <- read.csv(paste0(location, name, '.csv'))
  results[[i]] <- data
}

res <- do.call(rbind, results)
res <- res[order(res$P.Value),]
res <- left_join(res, naming, by = 'Outcome')
write.csv(res, 'Cox/basic_Bon.csv', row.names = F)
top <- res[which(res$P.Value < 0.0000031),] # 5273 passing

# Full results collation 
location <- 'Results/Cox/Results/'
files <- list.files(location, '_FO_FULL.csv')
results <- list()

for(i in 1:length(files)){
  name <- as.character(diseases[i])
  print(name)
  data <- read.csv(paste0(location, name, '_FULL.csv'))
  results[[i]] <- data
}

res2 <- do.call(rbind, results)
res2 <- res2[order(res2$P.Value),]
res2 <- left_join(res2, naming, by = 'Outcome')
write.csv(res2, 'Results/Cox/full_Bon.csv', row.names = F)

# Associations passing adjustment
top$retain <- paste0(top$Predictor, top$Outcome)
res2$retain <- paste0(res2$Predictor, res2$Outcome)
keep <- res2[which(res2$retain %in% top$retain),]
keep <- keep[which(keep$P.Value < 0.0000031),] # 3209
keep <- left_join(keep, naming, by = 'Outcome')
write.csv(keep, '/Associations_retained_Bon.csv', row.names = F)
length(unique(keep$Predictor)) 
table(keep$Outcome)

# Failure rate across top models
assocs <- as.data.frame(table(keep$Outcome))
assocs <- assocs[order(assocs$Freq),]
names(assocs) <- c('Incident disease', 'Associations')

fail <- keep[which(keep$cox.zph_protein < 0.05),]
fail2 <- keep[which(keep$cox.zph_global < 0.05),]
fails <- as.data.frame(table(fail$Outcome))
names(fails) <- c('Incident disease', 'Protein_fail')
assocs <- left_join(assocs, fails, by = 'Incident disease')
write.csv(assocs, 'Cox/Protein_failure_full_summary_table_Bon.csv', row.names = F)

###############################################################################################

### Create a summary table showing proteins associated with each disease and direction of effect 

# srun -p interactive --pty bash
library(tidyverse)
# Set disease and read in disease codes for the iteration (array)
diseases <- c('AL_FO', 'ALS_FO', 'BRAIN_FO', 'Breast_FO', 'Colorectal_FO', 'COPD_FO', 'CYS_FO', 'DEP_FO', 'Diab_FO', 'ENDO_FO', 'GYN_FO',
              'IBD_FO', 'IHD_FO', 'LIV_FO', 'LUNG_FO', 'LUP_FO', 'MS_FO', 'PD_FO', 'Prostate_FO', 'RA_FO', 'SCZ_FO', 'ST_FO', 'VD_FO', 'DEATH')

res <- data.frame(Outcome = 1:24, Associations = 1:24, Negative = 1:24, Positive = 1:24, Names_neg = 1:24, Names_pos = 1:24)

for(i in 1:length(diseases)){
  dis <- as.character(diseases[i])
  table <- keep[which(keep$Outcome %in% dis),]
  num <- dim(table)[1]
  
  neg <- table[which(table$Hazard.Ratio < 1),]
  pos <- table[which(table$Hazard.Ratio >= 1),]
 
  neg_list <- neg[,1]
  pos_list <- pos[,1]
  
  neg_list <- gsub("\\..*", "", neg_list)
  pos_list <- gsub("\\..*", "", pos_list)
  
  neg_list <- str_c(neg_list, collapse = ", ")
  pos_list <- str_c(pos_list, collapse = ", ")
  
  res[i,1] <- dis
  res[i,2] <- num
  res[i,3] <- dim(neg)[1] 
  res[i,4] <- dim(pos)[1]
  res[i,5] <- neg_list
  res[i,6] <- pos_list
}

res <- res[order(res[,2]),]
write.csv(res, 'Results/Cox/summary_table_associations_Bon.csv', row.names = F)
