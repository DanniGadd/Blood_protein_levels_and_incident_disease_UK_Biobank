##################################################################################################

### MAKE A PLOT SHOWING TOTAL NUMBER OF ASSOCIATIONS PER TRAIT FOR ALL TRAITS

##################################################################################################

# srun -p interactive --pty bash

library(MetBrewer)
library(readxl)
library(ggplot2)
library(ggthemes)
library(readxl)
library(tidyverse)
library(dplyr)

keep <- read.csv("/path_to_file.../Associations_retained_Bon.csv")

table <- as.data.frame(table(keep$Naming.x))

names(table) <- c("Phenotype", "Count")

plot_data <- table[order(-table$Count),]

pdf("/path_to_file.../Plot_showing_number_assocs_Bon.pdf", width = 18, height = 5)
ggplot(plot_data, aes(x = reorder(Phenotype, -Count), y = Count)) +
  geom_segment(aes(x = reorder(Phenotype, -Count),
                   xend = reorder(Phenotype, -Count),
                   y = 0, yend = Count),
               color = "#21918c", lwd = 1) +
  geom_point(size = 4, color="#3b528b", fill="#3b528b", pch = 21, bg = 4) +
  xlab("Number of associations with protein levels") +
  ylab("") +
  coord_flip() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  theme(panel.grid.major.y = element_blank()) + theme(legend.title = element_blank()) + theme(legend.position = 'None') +
  labs(y = "Number of associations with protein levels",
       x = "") +
  theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12, angle=0), axis.title = element_text(size=12)) 
dev.off()


##################################################################################################

### MAKE A PLOT SHOWING PRTOEIN ASSOCIATIONS WITH MULTIPLE MORBIDITIES

##################################################################################################

library(MetBrewer)
library(readxl)
library(ggplot2)
library(ggthemes)
library(readxl)
library(tidyverse)
library(dplyr) 

table <- read.csv("/path_to_file.../Associations_retained_Bon.csv")

table  <- table[c(1,11,3)]
names(table)[1] <- 'Protein'
table$Assoc <- ifelse(table[,3] < 1, 'neg', 'pos')
# > table(table$Assoc)
# 
# neg  pos
# 303 2898

freq <- as.data.frame(table(table$Protein))
length(unique(freq$Var1)) # 
freq <- freq[order(-freq$Freq),]
names(freq)[1] <- 'Protein'
table <- left_join(table, freq, by = 'Protein')
table <- table %>% group_by(Protein)
table <- table[order(-table$Freq),]
table <- as.data.frame(table)

top <- freq[which(freq$Freq >= 8),] 
list <- top$Protein
table <- table[which(table$Protein %in% list),]

table$Protein  <- ifelse(table$Protein %in% "IL6.P05231.OID21276.v1", "IL6_1.P05231.OID21276.v1", table$Protein)
table$Protein  <- ifelse(table$Protein %in% "IL6.P05231.OID20911.v1", "IL6_2.P05231.OID20911.v1", table$Protein)
table$Protein  <- ifelse(table$Protein %in% "IL6.P05231.OID20101.v1", "IL6_3.P05231.OID20101.v1", table$Protein)
table$Protein  <- ifelse(table$Protein %in% "IL6.P05231.OID20563.v1", "IL6_4.P05231.OID20563.v1", table$Protein)

table$Protein  <- ifelse(table$Protein %in% "TNF.P01375.OID21237.v1", "TNF_1.P01375.OID21237.v1", table$Protein)
table$Protein  <- ifelse(table$Protein %in% "TNF.P01375.OID20074.v1", "TNF_2.P01375.OID20074.v1", table$Protein)
table$Protein  <- ifelse(table$Protein %in% "TNF.P01375.OID20473.v1", "TNF_3.P01375.OID20473.v1", table$Protein)
table$Protein  <- ifelse(table$Protein %in% "TNF.P01375.OID20848.v1", "TNF_4.P01375.OID20848.v1", table$Protein)

table$Protein <- sub("\\..*", "", table$Protein)

library(ggplot2)
table$Freq <- as.factor(table$Freq)


table$Protein <- factor(table$Protein, levels = c("GDF15", "IL6_1", "IL6_2", "PLAUR", 'NEFL', 'ASGR1', 'CHI3L1', 'IL6_3', 'IL6_4',
                                                  "TNFRSF1A", 'CSF1', 'TNFRSF1B', "LGALS9", "CD74", "HAVCR2","CD300E" , "TNFRSF4",
                                                  "CD274", "CD27", "TNF_2", 'TNF_1', 'TNF_3', 'TNF_4', 'CCL7', 'ST6GAL1', 
                                                  'WFDC2', 'PRSS8', "IGFBP4", "BST2", 'HGF', 
                                                  "TNFRSF10A", "TIMP1",     "VSIG4",     "MMP12" ,    "MSR1" ,     "CDCP1",
                                                  "PGF" ,      "IL2RA" ,    "IL18BP",    "LAIR1"  ,   "LAMP3" ,    "CST3",
                                                  "CCL3" ,     "CXCL9" ,    "TNFRSF9" ,  "LILRB4"  ,  "CXCL13" ,   "MDK",
                                                  "TNFSF13",   "MZB1" ,     "TNFRSF14",  "ZBTB17" ,   "ITGA11" ,   "ITGAV"
                                                  ))
# 
table$Naming.x <- factor(table$Naming.x, levels = c("Death",
                                                    "COPD",
                                                    "Ischaemic heart disease" ,
                                                    "Liver disease",
                                                    "Ischaemic stroke",
                                                    "Type 2 diabetes",
                                                    "Systemic lupus erythematosus",
                                                    "Inflammatory bowel disease" ,
                                                    "Rheumatoid arthritis",
                                                    "Lung cancer" ,
                                                    "Vascular dementia",
                                                    "Alzheimer's dementia",
                                                    "Parkinson's disease",
                                                    "Amyotrophic lateral sclerosis",
                                                    "Multiple sclerosis",
                                                    "Cystitis"
                                                    ))


colour.scale <-  c('pos' = 'red', 'neg' = 'blue')

or <- table[order(table[,3]),]

pdf("/path_to_file.../plot_multimorbidity_associations_updated_heatmap_Bon_legend.pdf", width = 18, height =7)
ggplot(data = table, aes(x=Protein, y=Naming.x, fill = Hazard.Ratio)) +  theme_minimal() +
  geom_tile(colour="black", size = 0.25) +
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") +
  scale_alpha(range = c(0.9, 3)) +
  xlab("") + 
  ylab("") +
  theme(legend.title=element_blank(),axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size=12, angle = 90, hjust = 1),
        axis.title = element_text(size=12), legend.position = "right")
dev.off()

