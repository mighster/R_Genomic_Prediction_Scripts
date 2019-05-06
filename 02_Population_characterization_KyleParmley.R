# Phenomic Diversity
#population structure analyses reproducing the results of the widely-used computer program structure

#Install packages
install.packages(c("fields","RColorBrewer","mapplots"))
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
install.packages('adegenet')
install.packages("adegenet", repos="http://cran.rstudio.com/", dependencies=TRUE)
install.packages('ggplot2')
install.packages('stringi')
library(stringi)
library(adegenet)
library(dplyr)
library(ggplot2)
library(dendextend)

#Parallel computing
library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(4,type="SOCK")
registerDoSNOW(cl)
set.seed(8)
set.seed(10)


#Read in csv data
setwd('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data stuff')
GWAS_GD = read.table("GWAS_GD.txt", sep = '\t',header = T)
AllData <-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/GWAS/KGF_AdjustedBLUPsAllDays_thinned_Oct24_tall.csv", header = T)
colnames(AllData)
AllData <- subset(AllData, AllData$Day == "6")
AllData <- AllData[,c(1,5)]
GWAS_GD <- left_join(AllData, GWAS_GD, by="Entry")
GWAS_GD <- GWAS_GD %>% select(Entry,everything()) #move entry to the start of the df
GWAS_GD <- arrange(GWAS_GD, Entry) #sort new df based on master ranking column with no.1 being best
GWAS_GD[1:10,1:10]
GD <- GWAS_GD[1:292,6:ncol(GWAS_GD)]

metadata <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Randomizations Origin Data GWAS Names/Meta_data.csv")
metadata <-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Cluster_Summary/Cluster_metadata.csv", header = T)
meta <- GWAS_GD[1:nrow(GWAS_GD),c(1,3)]
### Population Structure Analysis ###
#################################################################################
######## Find optimal number of clusters using iterative K means approach #######
#################################################################################
library(stringi)
library(adegenet)

#HOLD ON HERE JUST A MINUTE. DOES THIS ITERATIVE K MEANS APPROACH PROVIDE SOMEWHAT DIFFERENT GROUPINGS EACH TIME??????
#THE GROUPINGS FROM THIS METHOD ARE VERY SIMILAR HOWEVER NOT IDENTICAL TO THE GROUPINGS OF THE NEI DISTANCE METHOD BELOW

      #  obj <- df2genind(GD, ploidy=2,sep = '/t') # 1. Make genind object to be used in further analysis
      grp <- find.clusters(obj, max.n=20, n.pca=200, scale=FALSE) # 2. try different values of k (interactive) using kmeans
      #The rule of thumb consists in increasing K until it no longer leads to an appreciable improvement of fit (i.e., to a decrease of BIC)
      phenogrp <- find.clusters(Day9_dataonly[,1:43], max.n=40, n.pca=200, scale=TRUE)
      Day9_pheno.18 <- c(Day9_metadata,phenogrp)
      
      ## number of accessions per group
      table(grp$grp)
      table(phenogrp$grp)
      
      grouping = data.frame(Day6_data$Entry,grp$grp)
      colnames(grouping)[1] = 'Name'
      colnames(grouping)[2] = 'kmeans.9'
      
      phenogrouping = data.frame(Day6_data$Name,phenogrp$grp)
      colnames(phenogrouping)[1] = 'name'
      colnames(phenogrouping)[2] = 'phenocluster'
      
      #Write out grouping of genotype
      write.csv(grouping, "C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/Population_Clustering_kmeans9_Dec12.csv",row.names = F)
      write.csv(phenogrouping, "C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Dendrogram/Pheno_Clustering_6groups.csv",row.names = F)
      
      getwd()
      grouping <- read.csv("SNP.8_grouping.csv")
      grouping <- read.csv("Population_Clustering_6groups.csv")
      
      metadata <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Randomizations Origin Data GWAS Names/Meta_data.csv")
      
#Read in BLUP file with explanatory information
      #setwd('/Volumes/GoogleDrive/My Drive/Phenomic_Diversity/Results')
      dat2 = df<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/GWAS/KGF_AdjustedBLUPsAllDays_thinned_Oct24.csv",header = T)
      table(dat2$Country)
#Merge grouping and explanatory info
      #dat3 = merge(metadata, dat2, by = 'Entry', sort = F)
      dat3 <- dat2
      dat2 <- metadata
#Modify dat3 so only keep; China, Japan, U.S., and make all other 'other'
      dat2$Country <- as.character(dat2$Country)
      str(dat2$Country)
      dat2$Country[dat2$Country == 'Algeria'] <- 'Other';dat2$Country[dat2$Country == 'France'] <- 'Other'
      dat2$Country[dat2$Country == 'Georgia'] <- 'Other';dat2$Country[dat2$Country == 'Korea'] <- 'South Korea'
      dat2$Country[dat2$Country == 'Morocco'] <- 'Other';dat2$Country[dat2$Country == 'Poland'] <- 'Other'
      dat2$Country[dat2$Country == 'Portugal'] <- 'Other';dat2$Country[dat2$Country == 'Taiwan'] <- 'Asia'
      dat2$Country[dat2$Country == 'Turkey'] <- 'Other';dat2$Country[dat2$Country == 'Uzbekistan'] <- 'Other'
      dat2$Country[dat2$Country == 'Vietnam'] <- 'Asia';dat2$Country[dat2$Country == 'Yugoslavia'] <- 'Other'
      #dat2$Country[dat2$Country == 'Russia'] <- 'Other';dat2$Country[dat2$Country == 'South Korea'] <- 'Other'
      #dat2$Country[dat2$Country == 'North Korea'] <- 'Other'
      table(dat2$Country)

      newmetadata<-cbind(metadata,dat2$Country)
      write.csv(newmetadata,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Cluster_Summary/Cluster_metadata1.csv")
#################################################################################            
###### Compute LD, Fst, and genetic distance parameters using NAM package  ######
#################################################################################
      
      install.packages('NAM')
      install.packages('colorspace')
      library(NAM) #nested association mapping - Alencar Xavier, Katy Rainey
      library(phylogram)   
      library(dendextend)
      library(circlize)
      library(colorspace)
      #Convert GD into matrix form 
      geno = as.matrix(GD)
      geno[1:10,1:10]
      grouping[1:10,1:2]
      
      
      #Compute Fst among subpopulations found using the Population Structure Analysis Methods (NAM Package)
      #Make fam vector
      grouping <- Day9_metadata$snp35k.8 
      str(grouping)
      fam = as.vector(grouping)
      fst = Fst(geno,fam) #Genetic variation associated with markers distributed among subpopulations. The function generates a plot for structure diagnosis.
      plot(fst$fst)
      
#####################################################################################################################################                  
      #Compute genetic distance using Nei Distance (however there are a total of 5 methods - see gdist documentation)
      #gdist = You must choose one genetic distance to calculate: choose gdist from:"Nei GST", "Hedrick G'ST", "Jost D", "WC Theta", "PhiST", "Chi2", or "NL dA"
      gdist = Gdist(geno, method = 1) #This function computes measures of genetic distances between populations using a genpop object. 
      
      
      fit<-hclust(gdist,method='ward.D') #Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it - see hclust documentation
      CLUST.9 <-  cbind(Day6_data[,1:2],cutree(fit, k = 9))
      
      write.csv(CLUST.9,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/Population_Clustering_snp35k9_Dec12.csv",row.names = F)

#####################################################################################################################################            
########              https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html                     ########
#####################################################################################################################################
      iris <- datasets::iris
      iris2 <- iris[,-5]
      iris[,5]
      
      cluster_label <- (GWAS_GD$Greek8)[order.dendrogram(dend)]
      cluster_label <- factor(c("1","2","3","4","5","6","7","8"))
      
      
      dend <- as.dendrogram(hclust(gdist,method='ward.D'))
      # order it the closest we can to the order of the observations:
      dend <- rotate(dend, 1:150)
      # Color the branches based on the clusters:
      dend <- color_branches(dend, k=8)
      #dend <- (groupLabels=cluster_label)
      
      
      # Manually match the labels, as much as possible, to the real classification of the flowers:
      labels_colors(dend) <-
        rainbow_hcl(8)[sort_levels_values(
          as.numeric(CLUST.8$`cutree(fit, k = 8)`)[order.dendrogram(dend)]
        )]
      
      
      # We shall add the flower type to the labels:
      labels(dend) <- paste(as.character(GWAS_GD$Country)[order.dendrogram(dend)],
                            "(",labels(dend),")", 
                            sep = "")
      
      labels(dend) <- paste(as.character(GWAS_GD$Country)[order.dendrogram(dend)],
                            "(",labels(dend),")", 
                            sep = "")
      
      # We hang the dendrogram a bit:
      dend <- hang.dendrogram(dend,hang_height=0.1)
      # reduce the size of the labels:
      dend <- set(dend, "labels_cex", 0.5)
      # And plot:
      par(mar = c(3,3,3,7))
      
      tiff("SNP8_composition.tiff", width = 20, height = 20, units = 'in', res = 300)
      plot(dend, 
           main = "8 Genotypic Clusters", 
           horiz =  TRUE,  nodePar = list(cex = .007))
      legend("topleft", legend = "Cluster 2",  fill = rainbow_hcl(8))
      colored_bars(rainbow_hcl(8), dend, rowLabels = cluster_label, "clusterzzz", horiz = FALSE)
      
      dev.off()
      
      tiff("SNP8_composition_Circle.tiff", width = 10, height = 8, units = 'in', res = 300)
      circlize_dendrogram(dend)
      dev.off()
      getwd()
 #####################################################################################################################################     
      
      dend1 <- color_branches(dend, k = 8)
      dend2 <- color_labels(dend, k = 8)
      par(mfrow = c(1,2))
      circlize_dendrogram(dend1, main = "Colored branches")
      plot(dend2, main = "Colored labels")
      plot(fit)
      
      plot(as.phylo(fit),cex = 0.5,show.tip.label = T) # converts an object into a tree of class "phylo". 
      #In an object of class "hclust", the height gives the distance between the two sets that are being agglomerated. So these distances are divided by two when setting the branch lengths of a phylogenetic tree.
      meta <- GDmerged[,c(1:3,6:14)]
      tree = cbind.data.frame(meta,cluster=cutree(fit, k = 8)) #Cuts a tree, e.g., as resulting from hclust, into several groups either by specifying the desired number(s) of groups or the cut height(s).
      #k	= an integer scalar or vector with the desired number of groups
      #h	= numeric scalar or vector with heights where the tree should be cut.
      table(tree)
      library(ape) #Analyses of Phylogenetics and Evolution; provides functions for reading, writing, manipulating, analysing, and simulating phylogenetic trees and DNA sequences
      try1 <- plot(as.phylo(fit),type='fan',label.offset=0.1,no.margin=TRUE,cex=0.5,show.tip.label = T)
      as.dendrogram(fit),type='fan',label.offset=0.1,no.margin=TRUE,cex=0.5,show.tip.label = T

################################################################################           
######################   PCA of Genotypic Information  #########################
################################################################################       
      
# PCA with function prcomp
      #pca1 = prcomp(GD[,2:ncol(GD)], scale. = TRUE)
      pca1 <- prcomp(geno, scale. = TRUE) #Performs a principal components analysis on the given data matrix and returns the results as an object of class prcomp
      summary(pca1)
      
      # loadings
      pca1_loading = pca1$x #the value of the rotated data (the centred data multiplied by the rotation matrix) 
      
      # add cluster info to pca
      #write.csv(pca1_loading,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/PCA.csv", row.names = T)
      
      pca1_loading <- cbind(GWAS_GD[,1:5], pca1_loading)
      names(pca1_loading)[1] <- "Name"
      dat4 <- pca1_loading
      #percent variance explaine
      
      PCA <- as.matrix(pca1)
      #write.csv(PCA,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/PCA_nov7.csv", row.names = T)
      #PCA <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/PCA_nov7.csv")
      
# make figure colored by pca
setwd('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/')
tiff("Subpopulation_pca_292.tiff",compression = "lzw", width = 4, height = 3, units = 'in', res = 300)
ggplot(dat4, aes(x = PC1, y = PC2, color = GenetBackground)) +
  geom_point(alpha=0.5) +
  labs(x = "PC1 (44.8%)", y = "PC2 (31.49%)")+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dev.off()

#make figure of subpopuation assignment by country

  dat5 = data.frame(table(dat2$Country,grouping$subpop))
  tiff("02_Subpopulation_composition.tiff", width = 6, height = 6, units = 'in', res = 300)
  ggplot(dat5, aes(Var2))+
    geom_bar(aes(weight=Freq,fill=Var1), width = 0.5) + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6))+
    theme_classic()+
    scale_color_brewer(palette="Set1")+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8,face="bold"))+
    coord_flip()+
    labs(x = "Cluster", y = 'Count ')+
    guides(fill=guide_legend(title="Country"))
  dev.off()

######## Not In Use ##########
# Not in use:  Alternative method using snapclust implements fast maximum-likelihood (ML) genetic clustering 
#Find optimal number of clusters
#x = snapclust.choose.k(20,obj, IC = BIC, IC.only = T)
#plot(1:20, x, xlab = "Number of clusters (k)",
# ylab = "AIC", type = "b", pch = 20, cex = 3)

## run EM algo with defined number of clusters
res <- snapclust(obj, 5, pop.ini = grp$grp ,hybrids = F)
# names(res)
#res$converged
#res$n.iter
## plot result
#compoplot(res)

#d.dapc <- dapc(obj, n.pca = 20, n.da = 2)
#scatter(d.dapc, clab = 0.85, col = funky(24),
#posi.da="topleft", posi.pca = "bottomleft", scree.pca = TRUE)
  
  
  
  ### Compute LD using synbreed package###
  #Import GM file
  #install.packages("synbreed",repos="http://r-forge.r-project.org")
  library(synbreed)
  #Preprocessing
  #GM = read.table("GWAS_GM.txt", sep = '\t',header = T)
  #colnames(GM)[2] = 'chr'
  #map = GM[,2:3]
  #row.names(map) = GM[,1]
  #geno = as.matrix(GD[,2:35581])
  #rownames(geno) = GD[,1]
  #Convert to gpData format
  #gpData<- create.gpData(pheno=NULL, geno=geno, map=map, covar= NULL, reorderMap=F,map.unit="bp")
  #gData_coded <- codeGeno(gpData, maf= 0.05, nmiss= 0.1, impute=F,
                          #impute.type = "random",  label.heter="alleleCoding",cores = 4) 
  
  #Compute pairwaise LD
  #pLD = pairwiseLD(gData_coded, chr = NULL, type = "matrix",use.plink=FALSE,
                   #ld.threshold=0, ld.window=99999, rm.unmapped = TRUE, cores=4)
  
  #plot(pLD)
  

