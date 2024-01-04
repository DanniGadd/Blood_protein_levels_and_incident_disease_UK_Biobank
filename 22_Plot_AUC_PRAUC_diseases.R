
### CREATE PLOTFOR FIG 3

# srun -p interactive --pty bash

library(MetBrewer)
library(readxl)
library(ggplot2)
library(ggthemes)
library(readxl)
library(tidyverse)
library(dplyr)
library(ggalt)

# install.packages("viridis")  
library("viridis")       
# https://waldyrious.net/viridis-palette-generator/

table <- read.csv("metricsTables_chosen_joint_cov_imputed_V4_results.csv")

##############################################################v

# Plot individuals diseases that had P for extended set 
My_Theme = theme(
  axis.title.x = element_text(size = 11),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11),
  axis.title.y = element_text(size = 11),
  strip.text = element_text(size = 11, face = "bold"),
  legend.text=element_text(size=11),
  legend.title=element_text(size=11, face = "bold"), legend.position = "none")

# Create matching index for labels
match <- data.frame(Labels = c('Age+Sex', 'Age+Sex+Lifestyle', 'Age+Sex+Expanded set', 'Age+Sex+Lifestyle+Expanded set', 'ProteinScore',
                               'Age+Sex+ProteinScore', 
                               'Age+Sex+Lifestyle+ProteinScore',
                               'Age+Sex+Lifestyle+Expanded set+ProteinScore'),
                    original = c('AgeSex','AgeSexCovs', 'ClinicalOnly', 'Clinicalrisk', 
                                 'ProteinScore Only', 'AgeSex_ProteinScore', 'AgeSexCovs_ProteinScore',
                                 'Clinicalriskprotein'))

# Diabetes
diseases <- c('Diab_FO')
res <- table[which(table$Outcome %in% diseases),]
AUC <- res[c(1:2)]
AUC <- AUC[-c(11,10,9),]
AUC <- AUC[match(match$original, AUC$X),]
AUC$X <- match$Labels
AUC$order <- 1:8
AUC$Naming = factor(AUC$X, levels=unique(AUC$X[order(AUC$order)]))

pdf('covariate_assessment/diab.pdf', width = 7.0, height = 2)
ggplot(AUC, aes(x=AUC, y=Naming)) +
  # geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color= c('slateblue4','slateblue4','slateblue4', 'slateblue4', 'salmon1', 'salmon1', 'salmon1', 'salmon1'), size=2, alpha=1) +theme_light() +
  # coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +ylab("") + My_Theme + xlim(0.5, 1.0)
dev.off()

# COPD
diseases <- c('COPD_FO')
res <- table[which(table$Outcome %in% diseases),]
AUC <- res[c(1:2)]
AUC <- AUC[-c(11,10,9),]
AUC$X <- sub('2', '', AUC$X)
AUC <- AUC[match(match$original, AUC$X),]
AUC$X <- match$Labels
AUC$order <- 1:8
AUC$Naming = factor(AUC$X, levels=unique(AUC$X[order(AUC$order)]))

pdf('covariate_assessment/copd.pdf', width = 7.0, height = 2)
ggplot(AUC, aes(x=AUC, y=Naming)) +
  # geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color= c('slateblue4','slateblue4','slateblue4', 'slateblue4', 'salmon1', 'salmon1', 'salmon1', 'salmon1'), size=2, alpha=1) +theme_light() +
  # coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +ylab("") + My_Theme + xlim(0.5, 1.0)
dev.off()

# AL
diseases <- c('AL_FO')
res <- table[which(table$Outcome %in% diseases),]
AUC <- res[c(1:2)]
AUC <- AUC[-c(11,10,9),]
AUC$X <- sub('10', '', AUC$X)
AUC <- AUC[match(match$original, AUC$X),]
AUC$X <- match$Labels
AUC$order <- 1:8
AUC$Naming = factor(AUC$X, levels=unique(AUC$X[order(AUC$order)]))

pdf('covariate_assessment/al.pdf', width = 7.0, height = 2)
ggplot(AUC, aes(x=AUC, y=Naming)) +
  # geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color= c('slateblue4','slateblue4','slateblue4', 'slateblue4', 'salmon1', 'salmon1', 'salmon1', 'salmon1'), size=2, alpha=1) +theme_light() +
  # coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +ylab("") + My_Theme + xlim(0.5, 1.0)
dev.off()

# IHD
diseases <- c('IHD_FO')
res <- table[which(table$Outcome %in% diseases),]
AUC <- res[c(1:2)]
AUC <- AUC[-c(11,10,9),]
AUC$X <- sub('11', '', AUC$X)
AUC <- AUC[match(match$original, AUC$X),]
AUC$X <- match$Labels
AUC$order <- 1:8
AUC$Naming = factor(AUC$X, levels=unique(AUC$X[order(AUC$order)]))

pdf('covariate_assessment/ihd.pdf', width = 7.0, height = 2)
ggplot(AUC, aes(x=AUC, y=Naming)) +
  # geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color= c('slateblue4','slateblue4','slateblue4', 'slateblue4', 'salmon1', 'salmon1', 'salmon1', 'salmon1'), size=2, alpha=1) +theme_light() +
  # coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +ylab("") + My_Theme + xlim(0.5, 1.0)
dev.off()

# DEATH
diseases <- c('DEATH')
res <- table[which(table$Outcome %in% diseases),]
AUC <- res[c(1:2)]
AUC <- AUC[-c(11,10,9),]
AUC$X <- sub('5', '', AUC$X)
AUC <- AUC[match(match$original, AUC$X),]
AUC$X <- match$Labels
AUC$order <- 1:8
AUC$Naming = factor(AUC$X, levels=unique(AUC$X[order(AUC$order)]))

pdf('covariate_assessment/death.pdf', width = 7.0, height = 2)
ggplot(AUC, aes(x=AUC, y=Naming)) +
  # geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color= c('slateblue4','slateblue4','slateblue4', 'slateblue4', 'salmon1', 'salmon1', 'salmon1', 'salmon1'), size=2, alpha=1) +
  theme_light() +
  # coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +ylab("") + My_Theme + xlim(0.5, 1.0)
dev.off()

# PD
diseases <- c('PD_FO')
res <- table[which(table$Outcome %in% diseases),]
AUC <- res[c(1:2)]
AUC <- AUC[-c(11,10,9),]
AUC$X <- sub('12', '', AUC$X)
AUC <- AUC[match(match$original, AUC$X),]
AUC$X <- match$Labels
AUC$order <- 1:8
AUC$Naming = factor(AUC$X, levels=unique(AUC$X[order(AUC$order)]))

pdf('covariate_assessment/pd.pdf', width = 7.0, height = 2)
ggplot(AUC, aes(x=AUC, y=Naming)) +
  # geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color= c('slateblue4','slateblue4','slateblue4', 'slateblue4', 'salmon1', 'salmon1', 'salmon1', 'salmon1'), size=2, alpha=1) +
  theme_light() +
  # coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +ylab("") + My_Theme + xlim(0.5, 1.0)
dev.off()


##############################################################v

# Plot staggered for each disease with difference
table <- read.csv("covariate_assessment/metricsTables_chosen_joint_cov_imputed_V4_results.csv")
table[,1] <- rep(c('m1','m6','m5','m2','m7', 'm3', 'm4', 'm8', 'd1', 'd2', 'd3'), 19)
match2 <- match
names(match2)[2] <- 'original'
match2$X <- c('m1','m2','m3','m4','m5', 'm6', 'm7', 'm8')
table <- left_join(table, match2, by = 'X')

# Join disease naming
library(readxl)
naming <- read_excel("naming_index.xlsx")
naming <- as.data.frame(naming)
names(naming)[1] <- 'Outcome'
table <- left_join(table, naming, by = 'Outcome')
AUC <- table[-3]
PRAUC <- table[-2]
AUC <- AUC[c(1,2,15)]
AUC <- reshape(AUC, idvar = "Naming", timevar = "X", direction = "wide")
PRAUC <- PRAUC[c(1,2,15)]
PRAUC <- reshape(PRAUC, idvar = "Naming", timevar = "X", direction = "wide")

My_Theme = theme(
  axis.title.x = element_text(size = 11),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11),
  axis.title.y = element_text(size = 11),
  strip.text = element_text(size = 11, face = "bold"),
  legend.text=element_text(size=11),
  legend.title=element_text(size=11, face = "bold"), legend.position = "none")

# DO NULL AUC 
AUC$Naming = factor(AUC$Naming, levels=unique(AUC$Naming[order(AUC$AUC.d1)]))
pdf('plot_AUC_basic.pdf', width = 5, height = 4)
scales::hue_pal()(2)
ggplot(AUC, aes(x = AUC.m1, xend = AUC.m6, y = Naming)) + 
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "#5b8124", colour_xend = "#bad744",
                dot_guide=TRUE, dot_guide_size=0.25) +
  scale_y_discrete(limits = rev, name = NULL) +
  xlab('AUC')+
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.09)) + xlim(0.5,1.0)
dev.off()


# DO FULL AUC 
AUC$Naming = factor(AUC$Naming, levels=unique(AUC$Naming[order(AUC$AUC.d1)]))
pdf('plot_AUC_full.pdf', width = 5, height = 4)
scales::hue_pal()(2)
ggplot(AUC, aes(x = AUC.m2, xend = AUC.m7, y = Naming)) + 
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "cadetblue4", colour_xend = "cadetblue3",
                dot_guide=TRUE, dot_guide_size=0.25) +
  scale_y_discrete(limits = rev, name = NULL) +
  xlab('AUC')+
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.09)) + xlim(0.5,1.0)
dev.off()

# DO EXTENDED AUC 
AUC$Naming = factor(AUC$Naming, levels=unique(AUC$Naming[order(AUC$AUC.d1)]))
pdf('plot_AUC_extended.pdf', width = 5, height = 4)
scales::hue_pal()(2)
ggplot(AUC, aes(x = AUC.m4, xend = AUC.m8, y = Naming)) + 
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "coral2", colour_xend = "lightsalmon1",
                dot_guide=TRUE, dot_guide_size=0.25) +
  scale_y_discrete(limits = rev, name = NULL) +
  xlab('AUC')+
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.09)) + xlim(0.5,1.0)
dev.off()


######################################################################################## 

death <- read.csv("death_results_full_10yr.csv")
disease <- read.csv("disease_results_full_10yr.csv")
disease_5yr <- read.csv("disease_results_full_5yr.csv")

table <- rbind(disease, death)
table$Onset <- '10'
disease_5yr$Onset <- '5'
table <- rbind(table, disease_5yr)
table[,1] <- rep(c('m1','m2','m3','m4','m5', 'd1', 'd2'), 20)
table <- table[which(table$Outcome %in% c('DEATH', 'ST', 'PD', 'RA', 'LUNG', 'LIV', 'IHD', 'Diab', 'COPD', 'AL', 'ALS')),]
AUC <- table[-3]
PRAUC <- table[-2]
AUC <- AUC[c(1,2,7)]
AUC <- reshape(AUC, idvar = "Naming", timevar = "X", direction = "wide")
PRAUC <- PRAUC[c(1,2,7)]
PRAUC <- reshape(PRAUC, idvar = "Naming", timevar = "X", direction = "wide")

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 24),
  strip.text = element_text(size = 20, face = "bold"),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24, face = "bold"), legend.position = "none")

# DO NULL AUC 
AUC$Naming = factor(AUC$Naming, levels=unique(AUC$Naming[order(AUC$AUC.d1)]))
pdf('plot_AUC_basic.pdf', width = 5, height = 4)
scales::hue_pal()(2)
ggplot(AUC, aes(x = AUC.m1, xend = AUC.m2, y = Naming)) + 
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "#5b8124", colour_xend = "#bad744",
                dot_guide=TRUE, dot_guide_size=0.25) +
  scale_y_discrete(limits = rev, name = NULL) +
  xlab('AUC')+
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.09)) + xlim(0.5,0.9)
dev.off()

# DO NULL PRAUC 
join <- AUC[c(1,7)]
PRAUC <- left_join(PRAUC, join, by = 'Naming')
PRAUC$Naming = factor(PRAUC$Naming, levels=unique(PRAUC$Naming[order(PRAUC$AUC.d1)]))

pdf('plot_PRAUC_null_updated_ROC_P.pdf', width = 5, height = 4)
scales::hue_pal()(2)
ggplot(PRAUC, aes(x = PRAUC.m1, xend = PRAUC.m2, y = Naming)) + 
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "#5b8124", colour_xend = "#bad744",
                dot_guide=TRUE, dot_guide_size=0.25) +
  scale_y_discrete(limits = rev, name = NULL) +
  xlab('PRAUC')+
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.09)) + xlim(0.2,0.7)
dev.off()

# DO FULL AUC 
AUC$Naming = factor(AUC$Naming, levels=unique(AUC$Naming[order(AUC$AUC.d1)]))
pdf('plot_AUC_full_updated_ROC_P.pdf', width = 5, height = 4)
scales::hue_pal()(2)
ggplot(AUC, aes(x = AUC.m4, xend = AUC.m5, y = Naming)) + 
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "cadetblue4", colour_xend = "cadetblue3",
                dot_guide=TRUE, dot_guide_size=0.25) +
  scale_y_discrete(limits = rev, name = NULL) +
  xlab('AUC')+
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.09)) + xlim(0.5,0.9)
dev.off()

# DO FULL PRAUC 
join <- AUC[c(1,8)]
PRAUC <- left_join(PRAUC, join, by = 'Naming')
PRAUC$Naming = factor(PRAUC$Naming, levels=unique(PRAUC$Naming[order(PRAUC$AUC.d1)]))

pdf('plot_PRAUC_full_updated_ROC_P.pdf', width = 5, height = 4)
scales::hue_pal()(2)
ggplot(PRAUC, aes(x = PRAUC.m4, xend = PRAUC.m5, y = Naming)) + 
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "cadetblue4", colour_xend = "cadetblue3",
                dot_guide=TRUE, dot_guide_size=0.25) +
  scale_y_discrete(limits = rev, name = NULL) +
  xlab('PRAUC')+
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.09))+ xlim(0.2,0.7)
dev.off()






