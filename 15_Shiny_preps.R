
### Prep data for shiny app integrations


##################################################

### FORMAT INPUT DATA: - protein network

### Format for network plot - 6744 assocs full model as input
assocs <- read.csv("/path_to_file.../Associations_retained_Bon.csv")

# Select same cols as before
assocs <- assocs[which(colnames(assocs) %in% c('Predictor', 'Naming.x', 'Hazard.Ratio', 'LCI', 'UCI','P.Value', 'No..of.Cases',
                                               'No..of.Controls', 'cox.zph_protein'))]
assocs <- assocs[c(1,9,2,3,4,5,6,7,8)]

names(assocs) <- c('Predictor', 'Outcome', 'HR', 'LCI', 'UCI', 'P.Value', 'N Cases', 'N Controls', 'Protein zph')

assocs$Protein <- sub("\\..*", "", assocs$Predictor)


# Set edge colours based on HR < 1 or > 1
assocs$edge_col <- ifelse(assocs[,3] < 1, 'blue', 'red')

# Set thickness of edges to be deviation from HR value of 1
assocs[,3] <- ifelse(assocs[,3] < 1, 1-assocs[,3], 1-(2-assocs[,3]))

# Save out for manju
write.csv(assocs, '/path_to_file.../assocs.csv', row.names = F)

##################################################

### FORMAT INPUT DATA: PANEL ONE - COX PH 16 year

library(tidyverse)
bind <- read.csv("/path_to_file.../iteration_sens_joint_assocs_allocation_controls_PA_corrected.csv")

names(bind) <- c('X', 'Predictor', 'Outcome', 'Iteration', 'Hazard.Ratio', 'LCI', 'UCI', 'P.Value', 'Controls', 'Cases',
                 'cox.zph_protein', 'cox.zph_global', 'Dis_name')

bind$Prot_name <- sub("\\..*", "", bind$Predictor)

# Create violations and significance variables
bind <- bind %>% mutate(Schoenfeld = case_when(
    bind$cox.zph_protein < 0.05 ~ 'Local P < 0.05',
    bind$cox.zph_protein >= 0.05 ~ 'Local P >= 0.05'
))

bind <- bind %>% mutate(Association = case_when(
    bind$P.Value < 0.0000031 ~ 'P sig',
    bind$P.Value >= 0.0000031 ~ 'P nonsig'
))

# Format disease names
bind <- bind %>% mutate(Dis_name = case_when(
    bind$Outcome == 'AL_FO' ~ 'Alzheimers dementia',
    bind$Outcome == 'COPD_FO' ~ 'COPD',
    bind$Outcome == 'VD_FO' ~ 'Vascular dementia',
    bind$Outcome == 'PD_FO' ~ 'Parkinsons disease',
    bind$Outcome == 'IBD_FO' ~ 'Inflammatory bowel disease',
    bind$Outcome == 'IHD_FO' ~ 'Ischaemic heart disease',
    bind$Outcome == 'ALS_FO' ~ 'Amyotrophic lateral sclerosis',
    bind$Outcome == 'DEATH' ~ 'Death',
    bind$Outcome == 'RA_FO' ~ 'Rheumatiod arthritis',
    bind$Outcome == 'DEP_FO' ~ 'Major depression',
    bind$Outcome == 'Diab_FO' ~ 'Type 2 diabetes',
    bind$Outcome == 'GYN_FO' ~ 'Gynaecological cancer',
    bind$Outcome == 'LUNG_FO' ~ 'Lung cancer',
    bind$Outcome == 'Colorectal_FO' ~ 'Colorectal cancer',
    bind$Outcome == 'LIV_FO' ~ 'Liver disease',
    bind$Outcome == 'BRAIN_FO' ~ 'Brain CNS cancer',
    bind$Outcome == 'ST_FO' ~ 'Ischaemic stroke',
    bind$Outcome == 'Prostate_FO' ~ 'Prostate cancer',
    bind$Outcome == 'SCZ_FO' ~ 'Schizophrenia',
    bind$Outcome == 'LUP_FO' ~ 'Systemic lupus erythematosus',
    bind$Outcome == 'CYS_FO' ~ 'Cystitis',
    bind$Outcome == 'Breast_FO' ~ 'Breast cancer',
    bind$Outcome == 'ENDO_FO' ~ 'Endometriosis',
    bind$Outcome == 'MS_FO' ~ 'Multiple sclerosis'
))

# Set variable to dictate colour based on failures of local PH
bind$col_variable <- ifelse(bind$Schoenfeld == 'Local P < 0.05', 'coral1', 'azure4')

# Set variable to dictate shape based on whether P full < 0.01
bind$shape_variable <- ifelse(bind$Association == 'P sig', 'circle', 'triangle')

bind <- bind[c(1:12,14:16,13,17,18)]

# Format iterations
bind <- bind[-which(bind$Iteration %in% 16),]

bind$Iteration <- ifelse(bind$Iteration == 17, 16, bind$Iteration)
write.csv(bind, '/path_to_file.../bind.csv', row.names = F)
bind <- na.omit(bind)
write.csv(bind, '/path_to_file.../bind_narm.csv', row.names = F)
