#Parallel computing
library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(20,type="SOCK")
registerDoSNOW(cl)
library(tidyverse)
library(rrBLUP)

Big_SNP_file <- read.table("C:/Users/jmshook/Desktop/BLINK/20110geno42080markers08262018.txt",head=T)
MG3list <- read.csv("C:/Users/jmshook/Desktop/MG3PIs.csv",head=T)
GD_292 <- read.csv("C:/Users/jmshook/Desktop/Kevin_GD.csv",head=T)
GD_292 <- GD_292[1:292,3:35550]
GD_292[GD_292 == 2] <- 0.5
GD_292[GD_292 == 1] <- 8
GD_292[GD_292 == 0] <- 1
GD_292[GD_292 == 8] <- 0
GD_292[1:30,1:30]
GDsub2[1:30,1:30]
t_Big_SNP_file <- t(Big_SNP_file)

GDsub <- t_Big_SNP_file[which(rownames(t_Big_SNP_file) %in% MG3list$PI),]


dim(MG3list)
GDsub2[1:2,1:2]
dim(GDsub)
GDsub2 <- GDsub[,which(colnames(GDsub) %in% colnames(GD_292_new))]
dim(GDsub2)
testing_names <- rownames(GDsub2)
GD_testing <- GDsub2[1:30,1:30]


write.csv(GDsub2,"C:/Users/jmshook/Desktop/MG3_GD.csv")
write.csv(GD_292,"C:/Users/jmshook/Desktop/Kevin_BEASTLY_292.csv")
#######John Edits#####
GD_292_names <- read.csv('Kevin_GD_names_only.csv')
GD_292_new <- t_Big_SNP_file[which(rownames(t_Big_SNP_file) %in% GD_292_names$PI),]
Pheno_ALL <- read.csv("C:/Users/jmshook/Desktop/KGF_AdjustedBLUPsAllDays_thinned_Mar25_tall_TRL_GR_CLchanged.csv")
######End John Edits#####
GD_292[1:30,1:30]
GDsub2[1:30,1:30]

Pheno_ALL <- read.csv("C:/Users/jmshook/Desktop/KGF_AdjustedBLUPsAllDays_thinned_Oct24_tall_TRL_GR.csv")
Day6_data <- subset(Pheno_ALL, Pheno_ALL$Day == "6")
Day9_data <- subset(Pheno_ALL, Pheno_ALL$Day == "9")
Day12_data <- subset(Pheno_ALL, Pheno_ALL$Day == "12")
colnames(Day9_data)

Pheno <- Day9_data[,34]
plot(Day9_data[,34])

Pheno
PhenoName= colnames(Pheno)
rownames(trait_accuracy) <- PhenoName

testing_df <- GDsub2
testing_df[1:10,1:10]
training_df <- GD_292_new # John Edit to add _new
training_df[1:10,1:10]

mean(Pheno)

#use all 292 for training in training matrix
summary(Pheno)
plot(Pheno)
m_train <- as.matrix(training_df)
str(m_train)
dim(m_train)
m_train[1:10,1:10]
#train the model
TL_answer<-mixed.solve(Pheno, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)

#set BLUP(u) as the marker_effects
marker_effect = TL_answer$u
marker_effects = as.matrix(marker_effect)
marker_effects[1:10,1]
BLUE <- as.vector(TL_answer$beta)

#set the new genotypes to testing matrix
m_valid <- as.matrix(testing_df)
m_valid <- as.matrix(training_df)

dim(testing_df)
dim(m_valid)
dim(marker_effects)
marker_effects[1:10,1]
m_valid[1:10,1:10]
#matrix multiply marker validation matrix times the marker effects
predicted_valid =  m_valid %*% marker_effects
predicted_train =  m_train %*% marker_effects

predicted_result <- as.vector((predicted_valid[,1])+BLUE)
summary(predicted_result)
summary(Pheno)
kyle <- (m_train %*% marker_effects)+BLUE
Pheno <- as.matrix(Pheno)
cor.test(Pheno,kyle, method="kendall")
cor.test(Pheno,kyle, method="spearman")
cor.test(Pheno,kyle, method="pearson")

cor(predicted_train+BLUE, Pheno, use="complete")
summary(predicted_train+BLUE)
summary(Pheno)

