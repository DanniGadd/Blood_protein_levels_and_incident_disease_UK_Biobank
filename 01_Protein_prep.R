###############################################################################################

## Protein data preparation 1.5k

###############################################################################################

# srun -p interactive --pty bash
# module load R
# R

library(tidyverse)
library(data.table)
library(pacman)
library(tidyverse)

olink_internal = readRDS('correct_1.5k_proteins_withdrawals.rds')
olink_internal <- as.data.frame(olink_internal)
names(olink_internal)[1] <- 'SampleID' 

## Remove related individuals
rel <- read.table("ukb_rel_a26041_s488168.dat", header = T)
dim(rel)
rel <- subset(rel, ID1 %in% olink_internal$SampleID & ID2 %in% olink_internal$SampleID) # 1276 pairs
write.csv(rel, "rel_index.csv")

list1 <- as.data.frame(rel$ID1)
list2 <- as.data.frame(rel$ID2)
names(list1) <- "X"
names(list2) <- "X"
list <- rbind(list1, list2)
list <- list$X
length(unique(list)) 
t <- as.data.frame(table(list))
t <- t[order(-t$Freq),]
m <- t[which(t$Freq > 1),]
n <- t[which(t$Freq < 2),]

 for(i in 1:1276){
   pair <- rel[i,1:2]
   pair_list <- c(pair$ID1, pair$ID2)
   if(pair$ID1 %in% m$list == TRUE){
     keep1 <- 'ID1'
   } else {
     keep1 <- 'No'
   }
   rel[i,6] <- keep1
 }

 for(i in 1:1276){
   pair <- rel[i,1:2]
   pair_list <- c(pair$ID1, pair$ID2)
   if(pair$ID2 %in% m$list == TRUE){
     keep1 <- 'ID2'
   } else {
     keep1 <- 'No'
   }
   rel[i,7] <- keep1
 }

 rel$V8 <- ifelse(rel$V6 == 'ID1' & rel$V7 == 'ID2', 'BOTH', 'No')

both <- rel[which(rel$V8 %in% 'BOTH'),]

# Exclude pairs with both IDs in multiple families
olink_ex <- olink_internal[-which(olink_internal$SampleID %in% both$ID1),]
olink_ex <- olink_ex[-which(olink_ex$SampleID %in% both$ID2),]
dim(olink_ex) 

# Find instances where individuals were part of multiple family structures and in pairs
# with others that were not - then select the one in single pairing

sub <- rel[-which(rel$V8 %in% 'BOTH'),]
sub1 <- sub[which(sub$V6 %in% 'ID1'),]
olink_ex <- olink_ex[-which(olink_ex$SampleID %in% sub1$ID1),]
dim(olink_ex)

sub2 <- sub[which(sub$V7 %in% 'ID2'),]
olink_ex <- olink_ex[-which(olink_ex$SampleID %in% sub2$ID2),]
dim(olink_ex)

# So individuals in multiple family structures have been removed
# We can then sample the unique pairs and remove one person at random
# We can be sure that the remaining person is not part of any other pairs, as
# we've filtered out those with multiple pairs.

# Now look at those in single pairs and select one at random
single <- rel[which(rel$V6 %in% 'No'),]
single <- single[which(single$V7 %in% 'No'),]
single <- single[which(single$V8 %in% 'No'),] 

for(i in 1:length(single$ID1)){
  pair <- single[i,1:2]
  pair_list <- c(pair$ID1, pair$ID2)
  keep <- sample(pair_list, size=1, replace=FALSE)
  single[i,9] <- keep
}

# Remove the final set of 1/2 pairs in single pair matches
olink_ex <- olink_ex[-which(olink_ex$SampleID %in% single$V9),]
which(olink_ex$SampleID %in% m$list) 
dim(olink_ex)
write.csv(olink_ex, 'excluolink_internal_rel_adj.csv', row.names = F)

## Check people's missingness for the protein measurements
olink_internal <- olink_ex
length <- length(rownames(olink_internal))
res <- data.frame(SampleID = 1:length, Complete = 1:length, Missing = 1:length)

for (i in 1:length){
  variable <- as.character(olink_internal$SampleID[i])
  individual <- olink_internal[which(olink_internal$SampleID %in% variable),]
  individual <- individual[1:1474]
  missing <- sum(is.na(individual))
  complete <- 1472 - missing
  res[i,1] <- variable
  res[i,2] <-  complete
  res[i,3] <- missing
  print(i)
}

# Order by those that have the most missing data
res <- res[order(res$Complete),]

# Index how many people have >10% missingness as a proportion
res$prop <- res$Missing / 1472
res$exclude <- ifelse(res$prop > 0.1, '1', '0')
table(res$exclude)

# Save off missingness summary
write.csv(res, "missingness_olink_internal_people_v2.csv", row.names = F)
keep <- res[which(res$exclude %in% 0),]
exclude <- res[which(res$exclude %in% 1),]
olink_internal <- olink_internal[which(olink_internal$SampleID %in% keep$SampleID),]

## Check protein missingness in the subset of individuals that have <10% missing data
olink_internal_sub <- olink_internal
length <- length(colnames(olink_internal_sub))
res2 <- data.frame(Variable = 1:length, Complete = 1:length, Missing = 1:length)

for (i in 3:length){
  variable <- as.character(colnames(olink_internal_sub)[i])
  individual <- olink_internal[,which(colnames(olink_internal) %in% variable)]
  missing <- sum(is.na(individual))
  complete <- 47600 - missing
  res2[i,1] <- variable
  res2[i,2] <-  complete
  res2[i,3] <- missing
  print(i)
}

# Save off missingness summary
res2 <- res2[-c(1:2),]
res2 <- res2[order(res2$Complete),]

# Index how many proteins have >10% missingness
res2$prop <- res2$Missing / 47600
res2$exclude <- ifelse(res2$prop > 0.1, '1', '0')
keep <- res2[which(res2$exclude %in% 0),]
IDs <- olink_internal_sub[,c(1:2)]
overlap <- which(colnames(olink_internal_sub) %in% keep$Variable)
test <- olink_internal_sub[, overlap]
test <- cbind(IDs, test)

# Save out file with exclusion summaries for proteins
write.csv(res2, 'missingness_olink_internal_proteins_v2.csv', row.names = F)
# Save out protein data at this stage
write.csv(test, 'olink_proteins_rel_excluded_missingness_excluded.csv', row.names = F)
# There are now 47600 people and 1468 protein measurements

###################################################################################

## Imputations of missing data
# KNN imputation

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("impute")
#https://www.rdocumentation.org/packages/impute/versions/1.46.0/topics/impute.knn

library(tidyverse)
library(impute)
olink_internal_sub2 <- read.csv('olink_proteins_rel_excluded_missingness_excluded.csv')
IDs <- olink_internal_sub2[c(1:2)]
data <- olink_internal_sub2[c(3:1470)]
data <- as.matrix(data)
rownames(data) <- IDs$pseudo_ind_id
data <- t(data)
print(dim(data))
imputed <- impute.knn(data)
imputed_data <- as.data.frame(t(imputed$data))
identical(as.character(rownames(imputed_data)), as.character(IDs$pseudo_ind_id))
imputed_data <- cbind(IDs, imputed_data)
write.csv(imputed_data, 'knn_imputed_processed_olink_internal_proteins.csv', row.names = F)
print('done')

#####################################################################

olink_internal_sub2 <- read.csv('knn_imputed_processed_olink_internal_proteins.csv')

# Save as RDS file for predictor score models
saveRDS(olink_internal_sub2, 'knn_imputed_processed_olink_internal_proteins.rds')

# Rank-based inverse normalisation and scaling of each protein
for(i in colnames(olink_internal_sub2)[c(3:1470)]){
  olink_internal_sub2[,i] <- qnorm((rank(olink_internal_sub2[,i], na.last='keep')-0.5)/sum(!is.na(olink_internal_sub2[,i])))
}

## Scale protein data
olink_internal_sub2[,3:1470] <- apply(olink_internal_sub2[,3:1470], 2, scale)
mean(olink_internal_sub2[,200], na.rm = T)
sd(olink_internal_sub2[,200], na.rm = T)
write.csv(olink_internal_sub2, 'knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv', row.names = F)

# Save as RDS file
olink <- read.csv('knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')
saveRDS(olink, 'knn_imputed_olink_proteins_normalised_scaled_related_excluded.rds')

#####################################################################

### PR COMP - PCA on proteins

#####################################################################

olink <- read.csv('knn_imputed_olink_proteins_normalised_scaled_related_excluded.csv')

library(tidyverse)
library(ggplot2)
library(readxl)
library(psych)
library(ggcorrplot)
library(cowplot)
library(tidyverse)

# Read in the protein file
prot <- olink
# Isolate just the protein columns of interest for PCA
prot <- prot[c(3:1470)]
# Try prcomp
library(factoextra)
res.pca <- prcomp(prot, scale = TRUE)
names(res.pca)
# [1] "sdev"     "rotation" "center"   "scale"    "x"
loadings <- res.pca$rotation # this provides the loadings
loadings[1:5,1:4]

# the matrix x has the principal component score vectors
dim(res.pca$x)

# Compute variance explained by componenets using standard deviations
std_dev <- res.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
head(prop_varex)
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val
write.csv(eig.val, "eig_values_prcomp.csv")

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation

### PLOT EIGENVALUES AND CUMULATIVE VARIANCE WITH CORR PLOT
var <- eig.val
var$num <- 1:1468
names(var)[1] <- "mes"
# var[1,1] <- 300
var$col <- ifelse(var$mes >= 1, "darkgrey", "orange")
dim(var[which(var$mes >= 1),]) # 159 eigenvalues greater than or equal to 1

pdf("proteins_eigen_vals.pdf")
ggplot(var, aes(num, mes, col = col)) +
  geom_point(size = 1) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
  xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) +
  scale_color_manual(values=c("#E69F00", "#999999"))
dev.off()

names(var)[3] <- "cum"

pdf("proteins_cum_var.pdf")
ggplot(var, aes(num, cum)) +
  geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
  xlab("Principal component") + ylab("Cumulative proportion")
dev.off()

# ggcorrplot
library(ggcorrplot)
corr <- cor(prot)
pdf("proteins_corr_plot.pdf", width = 50, height = 50)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           outline.col = "white") + theme(text = element_text(size = 0.2))
dev.off()

