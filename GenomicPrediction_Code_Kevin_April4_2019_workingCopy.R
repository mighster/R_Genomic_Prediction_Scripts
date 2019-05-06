library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(4,type="SOCK")
registerDoSNOW(cl)
library(rrBLUP)
library(dplyr)

Pheno_ALL <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Data/KGF_AdjustedBLUPsAllDays_thinned_Oct24_tall_TRL_GR.csv")
GWAS_GD <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data Stuff/GWAS_GD.txt",header = T)
Pheno_ALL[1:10,31:35]
GWAS_GD[1:5,1:10]
SNP_Markers <- (GWAS_GD[,4:ncol(GWAS_GD)])

Day6_data <- subset(Pheno_ALL, Pheno_ALL$Day == "6")
Day9_data <- subset(Pheno_ALL, Pheno_ALL$Day == "9")
Day12_data <- subset(Pheno_ALL, Pheno_ALL$Day == "12")

cbind(Day6_data[,2],GWAS_GD[1:292,3])

#impute missing SNP_Markers with A.mat
impute = (A.mat(SNP_Markers,max.missing=0.5,impute.method="mean",return.imputed=TRUE))
SNP_Markers_1_impute = impute$imputed

########################################################################
########################################################################
########################################################################
Pheno <- Day9_data %>% select(VOL)
########################################################################
training_entries <- as.matrix(sample(1:292, 146))
testing_entries <-setdiff(1:292,training_entries)
########################################################################
Pheno_training_data=as.matrix(Pheno[training_entries,])
SNP_training_data=as.matrix(SNP_Markers_1_impute[training_entries,])
########################################################################
Pheno_testing_data=as.matrix(Pheno[testing_entries,])
SNP_testing_data=as.matrix(SNP_Markers_1_impute[testing_entries,], K=NULL)
########################################################################
trained_model <- mixed.solve(y=Pheno_training_data, Z=SNP_training_data)
########################################################################
marker_effects <- as.matrix(trained_model$u)
BLUE <- as.vector(trained_model$beta)
########################################################################
predicted_test =  as.matrix(SNP_testing_data) %*% marker_effects
predicted_train =  as.matrix(SNP_training_data) %*% marker_effects
predicted_test_result <- as.vector((predicted_test[,1])+BLUE)
predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
########################################################################
summary(as.vector(predicted_train_result))
summary(predicted_test_result)
########################################################################
cor(as.vector(Pheno_testing_data), predicted_test_result, use="complete")
cor.test(as.vector(Pheno_training_data),predicted_train_result, method="spearman")
#cor.test(Pheno_training_data,predicted_train_result, method="spearman")
########################################################################
plot(Pheno_testing_data,predicted_test_result)
plot(Pheno_training_data,predicted_train_result)
