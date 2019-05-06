#Parallel computing
library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(4,type="SOCK")
registerDoSNOW(cl)
#Single-Trait Genomic Selection
install.packages("dplyr")
library(dplyr)
library(ggplot2)

print('What column is your data on in the pheno file?')

###############
#Pre-processing
#Import df of BLUP values 
setwd('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genomic Selection')
df = read.csv('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genomic Selection/CombinedAdjustedBLUPsNoCentroid.csv') #Use all BLUPs or subset from RFE analysis
df = df[,c(11:ncol(df))]#Remove descriptor vars


#Import genotype file
df_GD = read.table('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data stuff/GWAS_GD.txt',sep = '\t',header = T)
df_GD = df_GD[,c(1,5:35583)] #Remove any descriptor information
###############################################################################################################################
###############################################################################################################################
Pheno_ALL <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Data/KGF_AdjustedBLUPsAllDays_thinned_Oct24_tall_TRL_GR.csv")
Day6_data <- subset(Pheno_ALL, Pheno_ALL$Day == "6")
Day9_data <- subset(Pheno_ALL, Pheno_ALL$Day == "9")
Day12_data <- subset(Pheno_ALL, Pheno_ALL$Day == "12")

SNP_Markers <- (GWAS_GD[,3:ncol(GWAS_GD)])
SNP_Markers[1:292,1:10]
###############################################################################################################################
###############################################################################################################################
#Put snp and pheno df in same order
target <- df$name #Make pheno accession order the target
snp = df_GD[match(target, df_GD$name),] #Puts snp file in the same order as pheno file

snp[1:5,1:5]
#Prepare and install rrBLUP package
library(rrBLUP)
A <- A.mat(as.matrix(snp[,5:35580]),impute.method="EM", #Make relatedness matrix to be used in GBLUP
           n.core=4,max.missing=0.5,shrink = T)
rownames(A) = df[,1]  
  
  
 #Make df to save CV output for A matrix GBLUP
  x <- c("Fold", "Trait", "Accuracy","Vg","Ve","PEV")
  cv_results = data.frame(matrix(ncol=6,nrow=0))
  colnames(cv_results) = x
  
#Make df to save output for snp GBLUP
  x <- c("Fold", "Trait", "Accuracy","Vu","Ve")
  cv_results_mrk = data.frame(matrix(ncol=5,nrow=0))
  colnames(cv_results_mrk) = x
  
#Make initial df to save predictions: (BLUPs)
    gBLUP_predictions = data.frame(matrix(ncol=0,nrow=1490))
    gBLUP_predictions$Fold = rep(1:5, each = 298) #Insert CV fold number 
  

    n <- 10  #number of iterations
    m <- 78 #number of traits
    trait <- 3  #start column

for (i in 1:m){  
  
  
  for (i in 1:n) {   
    
    
#Prepare training/test data set
  test_train = sample(1:5, 298, replace=T,prob=c(.2,.2,.2,.2,.2)) #Make random sample of numbers for k-fold CV
  table(test_train)
  #Attach test_train to pheno and geno files
  pheno = cbind(test_train,df) 
  snp1 = cbind(test_train,snp)

  
###################
#Model training
#rrBLUP with A matrix
  #Fold 1
    #Add what trait and fold working on
      Fold = '1'; Trait = colnames(pheno[trait]) #Change the trait name
      #Partition training/testing and train model
      test <- which(pheno$test_train==1) #Fold number
      #Choose the pheno column in the df and mask observations used in validation set
      pheno_NA <- pheno[,trait];pheno_NA[test] <- NA ;data1 <- data.frame(y=pheno_NA,pheno$name) #Make new pheno df with trait value and genotype name
      #Fit gBLUP model using A matrix
      ans1 <- kin.blup(data1,K=A,geno="pheno.name",pheno="y",PEV = T)
      acc=round(cor(ans1$g[test],pheno[test,trait]),2) #Estimates accuracy of k-fold values
      #Save CV output
      newRow <- data.frame(Fold =Fold,Trait=Trait, Accuracy =acc ,Vg=ans1$Vg,Ve =ans1$Ve,PEV = mean(ans1$PEV))
      cv_results = rbind(cv_results,newRow)
      #Save gBLUP prediction
      new_preds_1 = data.frame(ans1$pred)
  #Fold 2
      Fold = '2'; Trait = colnames(pheno[trait]) 
      test <- which(pheno$test_train==2) 
      pheno_NA <- pheno[,trait];pheno_NA[test] <- NA ;data1 <- data.frame(y=pheno_NA,pheno$name)
      ans1 <- kin.blup(data1,K=A,geno="pheno.name",pheno="y",PEV = T)
      acc=round(cor(ans1$g[test],pheno[test,trait]),2)
      newRow <- data.frame(Fold =Fold,Trait=Trait, Accuracy =acc ,Vg=ans1$Vg,Ve =ans1$Ve,PEV = mean(ans1$PEV))
      cv_results = rbind(cv_results,newRow)
      new_preds_2 = data.frame(ans1$pred)
  #Fold 3
      Fold = '3'; Trait = colnames(pheno[trait])
      test <- which(pheno$test_train==3) 
      pheno_NA <- pheno[,trait];pheno_NA[test] <- NA ;data1 <- data.frame(y=pheno_NA,pheno$name)
      ans1 <- kin.blup(data1,K=A,geno="pheno.name",pheno="y",PEV = T)
      acc=round(cor(ans1$g[test],pheno[test,trait]),2)
      newRow <- data.frame(Fold =Fold,Trait=Trait, Accuracy =acc ,Vg=ans1$Vg,Ve =ans1$Ve,PEV = mean(ans1$PEV))
      cv_results = rbind(cv_results,newRow)
      new_preds_3 = data.frame(ans1$pred)
  #Fold 4
      Fold = '4'; Trait = colnames(pheno[trait])
      test <- which(pheno$test_train==4)
      pheno_NA <- pheno[,trait];pheno_NA[test] <- NA ;data1 <- data.frame(y=pheno_NA,pheno$name)
      ans1 <- kin.blup(data1,K=A,geno="pheno.name",pheno="y",PEV = T)
      acc=round(cor(ans1$g[test],pheno[test,trait]),2)
      newRow <- data.frame(Fold =Fold,Trait=Trait, Accuracy =acc ,Vg=ans1$Vg,Ve =ans1$Ve,PEV = mean(ans1$PEV))
      cv_results = rbind(cv_results,newRow)
      new_preds_4 = data.frame(ans1$pred)
  #Fold 5
      Fold = '5'; Trait = colnames(pheno[trait]) 
      test <- which(pheno$test_train==5) 
      pheno_NA <- pheno[,trait];pheno_NA[test] <- NA ;data1 <- data.frame(y=pheno_NA,pheno$name)
      ans1 <- kin.blup(data1,K=A,geno="pheno.name",pheno="y",PEV = T)
      acc=round(cor(ans1$g[test],pheno[test,trait]),2)
      newRow <- data.frame(Fold =Fold,Trait=Trait, Accuracy =acc ,Vg=ans1$Vg,Ve =ans1$Ve,PEV = mean(ans1$PEV))
      cv_results = rbind(cv_results,newRow)
      new_preds_5 = data.frame(ans1$pred)
  
      rm(pheno)
      rm(snp1)
      
  }
 trait = trait+1
}      
      
cv_results$Iteration = rep(1:n, each = 5)

acc_trait <- cv_results %>% group_by(Trait)
      
    GS_summary =acc_trait %>% summarise(
        acc = mean(Accuracy)
      )
      
      
      write.csv(blah, 'output_retestforposter.csv')
      

      
      
       rm(cv_results) #clean house
 x <- c("Fold", "Trait", "Accuracy","Vg","Ve","PEV")
 cv_results = data.frame(matrix(ncol=6,nrow=0))
 colnames(cv_results) = x      
      
      
#####Output predcition accuracy summary statistics from above
 library(dplyr)
 grouped <- group_by(cv_results, Trait)
 
 
 
 ####Compute heritability estimates
 GS_Heritability = grouped$Vg/(grouped$Vg+grouped$Ve)
 hist(GS_Heritability)
 grouped$Heritability = GS_Heritability
 
 
 GS_Summary = summarise(grouped, Acc_mean=mean(Accuracy), Acc_sd=sd(Accuracy),Vg_mean = mean(Vg),Vg_sd = sd(Vg),
                        Ve_mean = mean(Ve),Ve_sd = sd(Ve),PEV_mean = mean(PEV),PEV_sd = sd(PEV),
                        Heritability_mean = mean(Heritability),Heritability_sd = sd(Heritability))
 
 
 #Make figure of prediction accuracy values
 tiff('GWAS_GS_Prediction_Accuracy.tiff',res = 300, units = 'in',width = 6, height = 8)
 ggplot(data = GS_Summary, aes(x = reorder(Trait,Acc_mean), y = Acc_mean, ymin = Acc_mean-Acc_sd, ymax = Acc_mean+Acc_sd)) +
   geom_point(position = position_dodge(width = 0.15)) +
   geom_errorbar(position = position_dodge(width = 0.15), width = 0.2,size=0.8) +
   coord_flip() +scale_colour_manual(values = c("blue", "red","darkgreen")) +
   theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=0.5))+labs(x = "",y="Genomic Prediction Accuracy (r)")
 dev.off()
 
  #Make figure of heritability
 ggplot(data = GS_Summary, aes(x = reorder(Trait,Heritability_mean), y = Heritability_mean, ymin = Heritability_mean-Heritability_sd, ymax = Heritability_mean+Heritability_sd)) +
   geom_point(position = position_dodge(width = 0.15)) +
   geom_errorbar(position = position_dodge(width = 0.15), width = 0.2,size=0.8) +
   coord_flip() +scale_colour_manual(values = c("blue", "red","darkgreen")) +
   theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=0.5))+labs(x = "",y="Genomic Estimated Heritability")
 
 write.csv(GS_Summary, 'GS_Summary_output_retestforposter.csv')
 
 ############################################################
GS_Summary = read.csv("GS_Summary_output_500iterations_onlyImportantTraits.csv")
 tiff('GWAS_GS_Prediction_Accuracy_NB.tiff',res = 100, units = 'in',width = 6, height = 4)
 ggplot(data = GS_Summary, aes(x = reorder(Trait,Acc_mean), y = Acc_mean, ymin = Acc_mean-Acc_sd, ymax = Acc_mean+Acc_sd)) +
   geom_point(position = position_dodge(width = 0.15)) +
   geom_errorbar(position = position_dodge(width = 0.15), width = 0.2,size=0.8) +
   coord_flip() +scale_colour_manual(values = c("blue", "red","darkgreen")) +
   theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=0.5))+labs(x = "",y="Genomic Prediction Accuracy (R) 500 Interations")
 dev.off()
 