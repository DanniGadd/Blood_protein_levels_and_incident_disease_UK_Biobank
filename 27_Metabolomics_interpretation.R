###############################################################################################

## Metabolomics UKB - interpretation script

###############################################################################################
# 
# # ### Load sbatch information when run as array per disease
# args <- commandArgs(trailingOnly=TRUE)
# taskid <- as.numeric(args[1])
# print(paste0('taskid for this iteration is ', taskid))

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

### Identify median model for each score
prot <- read.csv("metricsTables_joint.csv")
metab <- read.csv("metricsTables_joint.csv")
prot_metab <- read.csv("metricsTables_joint.csv")

data <- list(prot, metab, prot_metab)
names <- c('Proteomics', 'Metabolomics', 'Proteomics_and_metabolomics')

list_chosen <- list()

for(i in 1:3){
  res <- data[[i]]
  name <- names[i]
  
  print(length(unique(res$iteration)))
  
  # Get differences for models
  diff <- res[grep('Diff', res$X),]
  diff <- diff[order(diff$AUC),]
  lengthd <- dim(diff)[1]
  
  # Select median model as proteinscore
  length <- dim(diff)[1]
  cut <- as.integer(length/2)
  median <- diff[cut,]
  median$Outcome <- "Type 2 Diabetes"
  median$Omics <- name
  median$available <- length(unique(res$iteration))
  
  # Get range of differences
  t1 <- res[which(res$iteration %in% median$iteration),]
  min <- min(diff$AUC)
  max <- max(diff$AUC)
  
  median$min_diff <- min
  median$max_diff <- max
  
  list_chosen[[i]] <- median
}

res <- do.call(rbind, list_chosen)

res <- res[order(-res$AUC),]

write.csv(top, '00_Chosen/Models_chosen.csv', row.names = F)

write.csv(res, '00_Chosen/Models_chosen_all.csv', row.names = F)








