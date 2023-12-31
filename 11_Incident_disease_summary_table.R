
############################################################################################

### Incident disease summary table 

############################################################################################

# srun -p interactive --pty bash
# Collate from tables central directory for all traits 
# module load R
# R

library(tidyverse)
files <- list.files("Results/Cox/Tables/", '.csv')
location <- 'Results/Cox/Tables/'
location_cut <- 'Results/Cox/tables_cut/'

# Read in and extract relevant info (i.e. smaller files minus protein data)
for(i in 1:length(files)){
  filename <- files[i]
  name <- sub( '.csv', '', filename)
  print(name)
  cox <- read.csv(paste0(location, name, '.csv'))
  data <- cox[,which(colnames(cox) %in% c('SampleID', 'pseudo_ind_id', "Age_assessment", "Sex", "dead", "aged", "Event", "age_death", "age_at_event", "tte", 'BMI', 'Dep', 'Alc', 'Smo', 'Edu'))]
  data$Outcome <- name
  write.csv(data, paste0(location_cut, name, '.csv'), row.names = F)
  print(dim(cox))
  cox = NULL
  data = NULL
}


# Read subset files back in as list and check max tte
results <- list()
files <- list.files(location_cut, '.csv')
for(i in 1:length(files)){
  name <- files[i]
  name <- gsub("\\..*", "", name)  
  results[[i]] <- read.csv(paste0(location_cut, name, '.csv'))
  print(i)
}

bind <- do.call(rbind, results)
max <- max(bind$tte, na.rm = T)

# Cases 
case <- bind[which(bind$Event == '1'),]
max <- max(case$tte, na.rm = T)

### Get the cases, controls and mean tte for each trait - basic model
output <- matrix(nrow = 1*length(files), ncol = 4)
output <- as.data.frame(output)
j=c(1:length(files))

for(i in 1:length(files)){
  df <- results[[i]]
  df <- as.data.frame(df)
  join <- df
  df <- join
  # generate metrics required
  cases <- df[which(df$Event == "1"),]
  case <- cases %>% filter(!tte == 'NA')
  case_count <- nrow(case)
  control <- df[which(df$Event == "0"),]
  control_count <- nrow(control)
  mean_tte <- mean(case$tte, na.rm = T) 
  mean_tte <- round(mean_tte, digits = 1)
  sd_tte <- sd(case$tte, na.rm = T)
  sd_tte <- round(sd_tte, digits = 1) 
  mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
  name_dis <- unique(df$Outcome)
  # output metrics 
  output[j[i],1] <- name_dis
  output[j[i],2] <- case_count
  output[j[i],3] <- control_count
  output[j[i],4] <- mean_sd
  names(output) <- c("Disease", "Cases", "Controls", "Mean time to event")
}

# order by cases low to high
sort1 <- output[order(output$Cases),]  


### Get the cases, controls and mean tte for each trait - full model

output <- matrix(nrow = 1*length(files), ncol = 4)
output <- as.data.frame(output)
j=c(1:length(files))

for(i in 1:length(files)){
  df <- results[[i]]
  df <- as.data.frame(df)
  join <- df
  
  PA <- read.csv('Censor_test/PA.csv')
  PA <- PA[which(colnames(PA) %in% c('SampleID', 'PA'))]
  join <- left_join(join, PA, by = 'SampleID')
  
  # remove missing covariates 
  join <- join[!is.na(join$Alc), ]
  join <- join[!is.na(join$Smo), ]
  join <- join[!is.na(join$Dep), ]
  join <- join[!is.na(join$Edu), ]
  join <- join[!is.na(join$BMI), ]
  join <- join[!is.na(join$PA), ]
  df <- join
  # generate metrics required
  cases <- df[which(df$Event == "1"),]
  case <- cases %>% filter(!tte == 'NA')
  case_count <- nrow(case)
  control <- df[which(df$Event == "0"),]
  control_count <- nrow(control)
  mean_tte <- mean(case$tte, na.rm = T) 
  mean_tte <- round(mean_tte, digits = 1)
  sd_tte <- sd(case$tte, na.rm = T)
  sd_tte <- round(sd_tte, digits = 1) 
  mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
  name_dis <- unique(df$Outcome)
  # output metrics 
  output[j[i],1] <- name_dis
  output[j[i],2] <- case_count
  output[j[i],3] <- control_count
  output[j[i],4] <- mean_sd
  names(output) <- c("Disease", "Cases", "Controls", "Mean time to event")
}

# order by cases low to high
sort <- output[order(output$Cases),]  


### JOIN THEM TOGETHER

basic <- sort1

full <- sort

names(basic) <- c("Disease", "Cases basic", "Controls basic", "Mean time to event basic")
names(full) <- c("Disease", "Cases full", "Controls full", "Mean time to event full")

join <- left_join(basic, full, by = "Disease")

write.csv(join, "Results/Cox/joint_Ns.csv", row.names = F)

###############################################################

### Plot tte distribution for each disease outcome 

library(ggpubr)
library(tidyverse)

# Read subset files back in as list and check max tte

location_new <- 'Results/Cox/tables_cut/'
results <- list()
files <- list.files(location_new, '.csv')

for(i in 1:length(files)){
  name <- files[i]
  name <- gsub("\\..*", "", name)  
  results[[i]] <- read.csv(paste0(location_new, name, '.csv'))
  print(i)
}

###############################################################

### Plot cumulative tte distribution for each disease outcome 

library(ggpubr)
library(tidyverse)

# Read subset files back in as list and check max tte
results <- list()
files <- list.files(location_new, '.csv')

for(i in 1:length(files)){
  name <- files[i]
  name <- gsub("\\..*", "", name)  
  results[[i]] <- read.csv(paste0(location_new, name, '.csv'))
  print(i)
}

bind <- do.call(rbind, results)
max <- max(bind$tte, na.rm = T)
plot_list <- list()
tte_list <- data.frame(Year = 1:15, Cases = 1:15, trait = 1:15)
tte_all <- list()

naming <- c("Alzheimer's dementia", 'Amyotrophic lateral sclerosis',
            'Brain/CNS cancer', 'Breast cancer', 'Colorectal cancer',
            'COPD', 'Cystitis', 'Death', 'Major depression', 'Type 2 diabetes', 'Endometriosis',
            'Gynaecological cancer', 'Inflammatory bowel disease', 'Ischaemic heart disease',
            'Liver fibrosis/cirrhosis', 'Lung cancer', 'Systemic lupus erythematosus',
            'Multiple sclerosis', "Parkinson's disease", 'Prostate cancer', 'Rheumatoid arthritis',
            'Schizophrenia', 'Ischaemic stroke', 'Vascular dementia')

format <- sub(".csv", "", files)

for(i in 1:length(files)){
  trait_format <- format[i]
  name_plot <- naming[i]
  # trait_format <- as.character(naming[[i]])
  df <- results[[i]]
  trait <- unique(df$Outcome)
  df <- as.data.frame(df)
  case <- df[which(df$Event %in% '1'),]
  # case <- cases %>% filter(tte > 0)
  N <- dim(case)[1]
  for(j in 1:15){
    sub <- case[which(case$tte <= j),]
    N_tte <- dim(sub)[1]
    tte_list[j,1] <- j
    tte_list[j,2] <- N_tte
    tte_list[j,3] <- trait
    
    prop <- (N_tte / N ) * 100
    tte_list[j,4] <- prop
    
    tte_list$V5 <- ifelse(tte_list$V4 < 80, 'blue', ifelse(tte_list$V4 > 90, 'grey', 'red'))
  }
  
  tte_all[[i]] <- tte_list
  
  p <- ggplot(tte_list, aes(x=Year, y=Cases)) + 
    geom_line() + geom_point() + theme_bw() + ggtitle(paste0(trait_format, ' - N = ', N)) +
    theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),
          axis.text = element_text(color='black',face='bold', size = 20),
          axis.title = element_text(color='black',face='bold', size = 20),
          plot.title = element_text(size = 20))+
    ylab('Cumulative cases')
  
  r <- ggplot(tte_list, aes(Year, Cases, fill = V5)) +     # Manually specifying colors
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("blue" = "#1b98e0",
                                 "red" = "salmon",
                                 'grey' = 'grey')) +
    theme_bw() + ggtitle(paste0(name_plot, ' - N = ', N)) +
    theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),
          axis.text = element_text(color='black',face='bold', size = 20),
          axis.title = element_text(color='black',face='bold', size = 20),
          plot.title = element_text(size = 20), legend.position = 'none')+
    ylab('Cumulative cases')
  # dev.off()
  
  plot_list[[i]] <- r
}

pdf(file = paste0("Cox/tte_by_traits_cum_V2.pdf"))
for (i in 1:length(files)) {
  print(plot_list[[i]])
}
dev.off()
