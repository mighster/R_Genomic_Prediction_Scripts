#Parallel computing
library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(4,type="SOCK")
registerDoSNOW(cl)
library(tidyverse)
library(rrBLUP)


Pheno_ALL <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Data/KGF_AdjustedBLUPsAllDays_thinned_Oct24_tall_TRL_GR.csv")
GM <- read.table("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data Stuff/298_nmis0.1_maf.05_GM.txt",header = T)
GD_292 <- read.table("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data Stuff/GD_292.txt",header = T)
GWAS_GD <- read.table("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data Stuff/GWAS_GD.txt",header = T)

Day6_data <- subset(Pheno_ALL, Pheno_ALL$Day == "6")
Day9_data <- subset(Pheno_ALL, Pheno_ALL$Day == "9")
Day12_data <- subset(Pheno_ALL, Pheno_ALL$Day == "12")

SNP_Markers <- (GWAS_GD[,3:ncol(GWAS_GD)])
SNP_Markers[1:292,1:10]
#SNP_Markers_1 <- t(SNP_Markers)
training_df <- SNP_Markers
m_train <- as.matrix(training_df)
#impute missing SNP_Markers with A.mat
impute=(A.mat(SNP_Markers,max.missing=0.5,impute.method="mean",return.imputed=TRUE))

A <- A.mat(as.matrix(snp[,5:35580]),impute.method="EM", #Make relatedness matrix to be used in GBLUP
           n.core=4,max.missing=0.5,shrink = T)

dim(impute$A)
colnames(impute)
SNP_Markers_1_impute=impute$imputed
dim(SNP_Markers_1_impute)
#############

colnames(Day9_data)
day_list <- list(Day6_data[,c(32:80)],Day9_data[,c(32:80)],Day12_data[,c(32:80)])
column_names <- colnames(Day6_data[,c(32:80)])
DayList <- c(6,9,12)
test_list <- c(146,60,90,120,150,180,210,240,270)
81-34
j=1
i=1
k=1
r=1
ratios=1
traits=2
cycles=1
days=1
accuracy = matrix(nrow=cycles, ncol=ratios)
all_accuracy <- data.frame()
colnames(accuracy) <- test_list
rownames(accuracy) <- c(1:cycles)

   for (j in c(1:days)){ #this  loop runs through each DAY, one at a time
      Day_data <- day_list[j]
      Day_data <- as.data.frame(Day_data)
      Day <- DayList[j]
      day_accuracy <- data.frame()
      
     for (i in 1:traits){  #this  loop runs through each TRAIT, one at a time
       Pheno=as.data.frame(Day_data[i])
       PhenoName=(colnames(Day_data[i]))
       trait_accuracy <- matrix()
       rownames(trait_accuracy) <- PhenoName
       #colnames(accuracy)[i] <- PhenoName
       #rownames(accuracy) <- column_names[i]
      
       for (k in c(1:ratios)){ #this  loop runs through each RATIO, one at a time
          #define the training and test populations
          #training-60% validation-40%
         b <- test_list[k]
         ratio_number <- b
         accuracy = matrix(nrow=cycles)
         #accuracy <- as.data.frame(accuracy)
        for(r in 1:cycles) {#this loop runs through each CYCLE, one at a time
          train= as.matrix(sample(1:292, b))
          test<-setdiff(1:292,train)
          Pheno_train=as.data.frame(Pheno[train,])
          m_train=as.data.frame(SNP_Markers_1_impute[train,])
          Pheno_valid=as.data.frame(Pheno[test,])
          m_valid=(SNP_Markers_1_impute[test,])
          dim(SNP_Markers_1_impute)
          TL=(Pheno_train[,1])
          TL_answer<-mixed.solve(TL, Z=m_train, K=impute$A, SE = FALSE, return.Hinv=FALSE)
          dim(impute)
          
          #set BLUP(u) as the marker_effects
          marker_effect = TL_answer$u
          dim(marker_effect)
          marker_effects = as.matrix(marker_effect)
          e = as.matrix(TL_answer$u)
          pred_TRL_valid =  m_valid %*% e
          pred_TRL_train =  as.matrix(m_train) %*% e
          
          TRL_valid = Pheno_valid[,1]
          summary(TRL_valid)
          accuracy[r,1] <- cor(pred_TRL_valid, TRL_valid, use="complete")
          cor.test(pred_TRL_valid,TRL_valid, method="spearman")
        }
         mean_accuracy <- as.data.frame(mean(accuracy))
         row.names(mean_accuracy) <- PhenoName
         colnames(mean_accuracy) <- b
         
         #trait_accuracy <- as.data.frame(trait_accuracy)
         trait_accuracy <- cbind(mean_accuracy,trait_accuracy)
         day6_accuracy <- rbind(day_accuracy,trait_accuracy)
       }
       all_accuracy <- rbind(all_accuracy,day6_accuracy)
       
      # rownames(all_accuracy[1,]) <- PhenoName
     }
  }
#write.csv(all_accuracy,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#3/Genomic Prediction/all_accuracy_feb.csv")
i=1
b=146
Pheno=as.data.frame(Day_data[i])
PhenoName=(colnames(Day_data[i]))
train<-as.matrix(sample(1:292, b))
test<-setdiff(1:292,train)
Pheno_train=as.data.frame(Pheno[train,])
m_train=as.data.frame(SNP_Markers_1_impute[train,])
dim(m_train)
Pheno_valid=as.data.frame(Pheno[test,])
m_valid=(SNP_Markers_1_impute[test,])
dim(m_valid)
TL=(Pheno_train[,1])
TL_answer<-mixed.solve(TL, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
#set BLUP(u) as the marker_effects
marker_effect = TL_answer$u

marker_effects = as.matrix(marker_effect)
BLUE <- as.vector(TL_answer$beta)

predicted_valid =  as.matrix(m_valid) %*% marker_effects
predicted_train =  as.matrix(m_train) %*% marker_effects
predicted_valid_result <- as.vector((predicted_valid[,1])+BLUE)
predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
PhenoName
summary(predicted_train_result)
summary(predicted_valid_result)
cor.test(predicted_train,predicted_valid, method="spearman")
cor(predicted_train, predicted_valid, use="complete")


summary(Pheno)
kyle <- (m_train %*% marker_effects)+BLUE
Pheno <- as.matrix(Pheno)


cor.test(Pheno,kyle, method="pearson")

cor(predicted_train+BLUE, Pheno, use="complete")
summary(predicted_train+BLUE)
summary(Pheno)

