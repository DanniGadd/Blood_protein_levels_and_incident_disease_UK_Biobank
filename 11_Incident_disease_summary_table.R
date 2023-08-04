
############################################################################################

### Incident disease summary table 

############################################################################################

# srun -p interactive --pty bash

# module load R

# R

library(tidyverse)

files <- list.files("/path_to_file.../cox_tables/", '.csv')

location <- '/path_to_file.../cox_tables/'

location_cut <- '/path_to_file.../tables_cut/'

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

write.csv(sort, file = '/path_to_file.../Collated_tables/summary_N_table.csv', row.names =F)


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
tte_list <- data.frame(Year = 1:16, Cases = 1:16, trait = 1:16)
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
  for(j in 1:16){
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
  
  # pdf('/path_to_file.../bar_test.pdf')
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

pdf(file = paste0("/path_to_file.../tte_by_traits_cum.pdf"))
for (i in 1:length(files)) {
  print(plot_list[[i]])
}
dev.off()
