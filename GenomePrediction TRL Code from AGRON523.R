#Genomic Selection R Code

#This first line of code defines the working directory for R.
#setwd("C:/Users/falk/Google Drive/PhD/Past Courses/AGRON 523/Genomic Selection")
setwd('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data stuff')

GWAS_GD = read.table("GWAS_GD.txt", sep = '\t',header = T)
impute <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data stuff/GWAS_GD_imputed.csv")
SNP_Markers_1_impute <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data stuff/SNP_Markers_1_impute.csv")
#load rrBLUP
impute
library(rrBLUP)
memory.limit()
#load data files
#SNP_Markers <- as.matrix(read.table(file="24_GD.txt"), header=FALSE)
#SNP_Markers <- as.matrix(read.table("GD_292.txt", sep = '\t',header = T))
colnames(GWAS_GD)
SNP_Markers <- (GWAS_GD[,3:ncol(GWAS_GD)])
head(SNP_Markers)
SNP_Markers_1 <- t(SNP_Markers)
#Pheno <-as.matrix(read.delim(file ="GS_Pheno.txt", header=TRUE))
Pheno_ALL <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Data/KGF_AdjustedBLUPsAllDays_thinned_Oct24_tall_TRL_GR.csv")
Day6_data <- subset(Pheno_ALL, Pheno_ALL$Day == "6")
Day9_data <- subset(Pheno_ALL, Pheno_ALL$Day == "9")
Day12_data <- subset(Pheno_ALL, Pheno_ALL$Day == "12")

colnames(Day9_data)
Pheno <- as.matrix(Day9_data[,34])
#dimensions of the matrix
dim(SNP_Markers)
dim(SNP_Markers_1)
dim(Pheno)
#############
str(SNP_Markers_1)
#impute missing SNP_Markers with A.mat
impute=A.mat(SNP_Markers,max.missing=0.5,impute.method="mean",return.imputed=TRUE)
colnames(impute)
SNP_Markers_1_impute=impute$imputed
dim(SNP_Markers_1_impute)
#############

#define the training and test populations
#training-60% validation-40%
train= as.matrix(sample(1:292, 291))
test<-setdiff(1:292,train)
Pheno_train=Pheno[train,]
m_train=SNP_Markers_1_impute[train,]
Pheno_valid=Pheno[test,]

m_valid=SNP_Markers_1_impute[test,]
m_valid<-m_valid[,2:ncol(m_valid)]
str(m_valid)
#############

#This code is for a single run of RR-BLUP for a particular trait (TRL in this case.)
str(Pheno_train)
TL=(Pheno_train)
TL_answer<-mixed.solve(TL, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
TRL = TL_answer$u
e = as.matrix(TRL)
pred_TRL_valid =  m_valid %*% e
pred_TRL=(pred_TRL_valid[,1])+TL_answer$beta
pred_TRL
TRL_valid = Pheno_valid
TRL_accuracy <-cor(pred_TRL_valid, TRL_valid, use="complete" )
TRL_accuracy
#############

#For 500 iterations run the following code. This is set up for TRL. 
traits=1
cycles=100
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:300, 180))
  test<-setdiff(1:300,train)
  Pheno_train=Pheno[train,]
  m_train=Markers_1_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_1_impute[test,]
  
  TL=(Pheno_train[,1])
  TL_answer<-mixed.solve(TL, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  TRL = TL_answer$u
  e = as.matrix(TRL)
  pred_TRL_valid =  m_valid %*% e
  pred_TRL=(pred_TRL_valid[,1])+TL_answer$beta
  pred_TRL
  TRL_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_TRL_valid, TRL_valid, use="complete" )
}

mean(accuracy)
