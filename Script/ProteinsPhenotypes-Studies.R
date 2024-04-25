
# Libraries
library(tidyverse)
library(moments) #Kurtosis &Skewness
library(dummy)
library(ggplot2)
library(umap)
library(dbscan)
library(ConsensusClusterPlus)


set.seed(789)




#### DATA INGESTION ####

#Working directory
  rm(list = ls())
  getwd()
  path = "/Users/alvarogarciamunoz/Documents/DATOS"
  setwd(path)

  
  
#File reading
  path <- "Data"
  fn_phenotypes<- "phenotypes.tsv"
  fn_proteins <- "Olink.NPX_values.Inflammatory-panels.csv"

  phenotypes <- read_tsv(paste(path,fn_phenotypes, sep="/"))
  proteins <- read.csv(paste(path,fn_proteins,sep="/"))
  
  
  
  
  
#### DATA AGGREGATION ####
  
#EXPLORATORY DATA ANALYSIS#
  
#UNIVARIANT
  
  
  
  #Search NA's
  summary(phenotypes) #Phenotypes  
  
  for(i in 1:ncol(proteins)){ #Boucle for Proteins
    NA.index <- which(is.na(proteins[,i]))
    if(length(NA.index) != 0){
      print(paste('The protein', colnames(proteins)[i], 'has a missing value in index', NA.index, sep = ' '))
    }
  }

#NA correction or elimination?  

  
  
  #Phenotypes plotting
  par(mfrow = c(2,2))
  temp <- phenotypes[, -which(colnames(phenotypes)=='Age')] #Age is the only phenotype with a non-discrete distribution.
  
  for(i in 2:ncol(temp)){
    tempt <- table(temp[,i]) #Table of phenotype
    print(tempt) #Obtain table
    print(prop.table(tempt)) #Obtain prop table
    barplot(tempt, main = paste('Variable', colnames(temp)[i]),
            xlab = 'Levels', ylab = 'Quantity') #Shows a barplot
  }
  
  par(mfrow = c(1,1))
  temp <-phenotypes[, which(colnames(phenotypes) == 'Age')]
  plot(density(temp$Age), main = 'Density of variable Age')
  
  
  
  #Continue distribution analysis
  
  temp <- phenotypes[, which(colnames(phenotypes) == 'Age')]
  index <- sample(2:ncol(proteins), 15, replace = FALSE)
  
  temp <- cbind(proteins[,index], temp)
  par(mfrow=c(4,4))
  
    #Density (sample)
  for (x in 1:ncol(temp)){
    plot(density(temp[,x]), main = paste("Density of protein number", colnames(temp)[x], sep = " "))
  }
  
  
    #Quartiles, Mean, MAD, CV & Outliers (sample)
  for (i in 1:ncol(temp)){ 
    tempt <- fivenum(temp[,i])
    q <- paste('-The quartiles of ', colnames(temp)[i], 
               '. Min:',tempt[1], 
               '. 1Q:' ,tempt[2],
               '. Median:', tempt[3],
               '. 3Q:', tempt[4],
               '. Máx:', tempt[5]
               , sep = '')
    m <- paste('-The mean of', colnames(temp)[i], 'is', mean(temp[,i], sep= ' '))
    mad <- paste("-The median absolute deviation of",colnames(temp)[i], "is", mad(temp[,i]), sep= " ")
    
    tempt <- sd(temp[,i])/mean(temp[,i])
    cv <- paste("-The coefficient of variation of",colnames(temp)[i], "is", tempt, sep= " ")
    title <- paste('Analysis of variable ', colnames(temp)[i], ':', sep = '')
    cat(title ,q, m, mad, cv, sep = "\n")   

    boxplot(temp[,i], main = paste("Boxplot of Variable", colnames(temp)[i], sep = " "))
  } #We don't exclude the outliers, we need to study them...
  
  
    #Kurtosis & Skewness (all)
  kurtosis <- list('Normal' = c(), 'Leptokurtic' = c(), 'Platykurtic' = c())
  skewness <- list('Symmetric' = c(), 'Right' = c(), 'Left' = c())
  for (i in 2:ncol(proteins)){
    temp <- kurtosis(proteins[,i])
    tempt <- skewness(proteins[,i])
    
    if(temp == 3){
      kurtosis[[1]] <- cbind(kurtosis[[1]], colnames(proteins)[i])
    }
    else if(temp >3){
      kurtosis[[2]] <- cbind(kurtosis[[2]], colnames(proteins)[i])
    }
    else{
      kurtosis[[3]] <- cbind(kurtosis[[3]], colnames(proteins)[i])
    }
    
    
    if(temp == 0){
      skewness[[1]] <- cbind(skewness[[1]], colnames(proteins)[i])
    }
    else if(temp >0){
      skewness[[2]] <- cbind(skewness[[2]], colnames(proteins)[i])
    }
    else{
      skewness[[3]] <- cbind(skewness[[3]], colnames(proteins)[i])
    }
  }
  
  N <- paste(kurtosis[[1]], sep = ", ")
  L <- paste(kurtosis[[2]], sep = ", ")
  P <- paste(kurtosis[[3]], sep = ", ")
  
  cat('Kurtosis lengths:', 'Normal:', length(kurtosis[[1]]), 'Leptokurtic:', length(kurtosis[[2]]), 'Platykurtic:', length(kurtosis[[3]]), sep = '\n')
  
  cat('Kurtosis results:', 'Normal: ', N, 'Leptokurtic: ', L, 'Platykurtic: ', P)
  
  N <- paste(skewness[[1]], sep = ", ")
  L <- paste(skewness[[2]], sep = ", ")
  P <- paste(skewness[[3]], sep = ", ")
  
  cat('Skewness lengths:', 'Symmetric:', length(skewness[[1]]), 'Right:', length(skewness[[2]]), 'Left:', length(skewness[[3]]), sep = '\n')
  
  cat('Skewness results:', 'Symmetric: ', N, 'Rigth: ', L, 'Left: ', P)
  
  #BIVARIATE#
  #Proteins Covariance
  temp <- cov(proteins[,-1])
  temp2 <- eigen(temp)
  newbase <- temp2$vectors
  m<-c()
  for (i in 1:nrow(proteins)){
    pr <- newbase * proteins[i,-1]
    m <- rbind(m, pr)
  }
  
  proteins.nb <-cbind("id" = proteins$id, m) #No covariance base proteins
  
  temp <- phenotypes[, which(colnames(phenotypes) == 'Age')]
  index <- sample(2:ncol(proteins.nb), 15, replace = FALSE)
  
  temp <- cbind(proteins.nb[,index], temp)
  par(mfrow=c(4,4))
  
  #Density (sample)
  for (x in 1:ncol(temp)){
    plot(density(temp[,x]), main = paste("Density of protein number", colnames(temp)[x], sep = " "))
  }
  
  tempt <- sample(1, 1:ncol(proteins))
  model <- lm(colnames(proteins)[, tempt]~.)
  
  #Phenotype chr management
  temp <- dummy(phenotypes[,-1], int = TRUE)
  temp2 <- cbind(phenotypes[,-ncol(phenotypes)], temp)[,-1]
  
  
  #The structure stand, lets use normalized
  phenotypes.int<- cbind(phenotypes[,1:ncol(phenotypes)-1], temp)

#Data grouping - Labels, em, Data. (Labels = dependent variable, em = extended matrix, Data = independent variables )

  labels <- subset(phenotypes.int, select = c(id, Diagnosis))
  labels$Diagnosis <- as.factor(labels$Diagnosis)
  
  em <- merge(proteins.nb, phenotypes.int, by="id") 
  data <- em[,-which(colnames(em)=="Diagnosis")]

#Normalize
  em.norm <-scale(em[,-ncol(proteins)])
  data.norm <-scale(data[,-ncol(proteins)])
  
  
#Normal
  
  
  
#### PIPELINE A ####
  
#UMAP 1
  em.split <- split(em, em$Diagnosis)
  
  em.umap <-umap(em[,-1])
  data.umap <-umap(data.norm)
  proteins.umap <- umap(proteins[,-1])
  
  coordinates.em.umap <- data.frame(em.umap$layout, labels = labels$Diagnosis)
  coordinates.data.umap <- data.frame(data.umap$layout, labels = labels$Diagnosis)
  coordinates.proteins.umap <-data.frame(proteins.umap$layout, labels = labels$Diagnosis)
  
    di1<- umap(em.split[[1]][,-1])
    di2<- umap(em.split[[2]][,-1])
    di3<- umap(em.split[[3]][,-1])
    di4<- umap(em.split[[4]][,-1])
  
  coordinates.di1 <-data.frame(di1$layout, labels = em.split[[1]]$Diagnosis)
  
  ggplot(data = coordinates.di1, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="UMAP Reduction of Extended Matrix")
  
  
  coordinates.di2 <-data.frame(di2$layout, labels = em.split[[2]]$Diagnosis)
  
  ggplot(data = coordinates.di2, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="UMAP Reduction of Extended Matrix")
  
  coordinates.di3 <-data.frame(di3$layout, labels = em.split[[3]]$Diagnosis)
  
  ggplot(data = coordinates.di3, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="UMAP Reduction of Extended Matrix")
  
  coordinates.di4 <-data.frame(di4$layout, labels = em.split[[4]]$Diagnosis)
  
  ggplot(data = coordinates.di4, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="UMAP Reduction of Extended Matrix")
  
#Plotting 1st 3 matrix
  
  ggplot(data = coordinates.em.umap, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="UMAP Reduction of Extended Matrix")
  
  ggplot(data = coordinates.data.umap, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="UMAP Reduction of Data")
  
  ggplot(data= coordinates.proteins.umap, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="UMAP Reduction of Proteins")
  

#Resultados raros, probar con diferentes distancias -  Pearson, Minkowski, Kullback-Leibler y validación cruzada
  
  
#HDBSCAN. Unsupervised - Clustering
  
  #Clusters
  em.hdbscan <- hdbscan(em[,-1], minPts = 5)
  data.hdbscan <- hdbscan(data.norm, minPts = 3)
  proteins.hdbscan <- hdbscan(proteins[,-1], minPts = 5)
  
  coordinates.em.hdbscan <- data.frame(proteins.umap$layout, labels = as.factor(em.hdbscan$cluster))
  coordinates.data.hdbscan <- data.frame(data.umap$layout, labels = as.factor(data.hdbscan$cluster))
  coordinates.proteins.hdbscan <- data.frame(proteins.umap$layout, labels = as.factor(proteins.hdbscan$cluster))
    


  #Plotting
  
  ggplot(data = coordinates.em.hdbscan, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="HDBSCAN cluster in High Dim of Extended Matrix", subtitle = "HDBSCAN EM clustering")
  
  ggplot(data = coordinates.data.hdbscan, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="HDBSCAN cluster in High Dim of Data", "HDBSCAN Data clustering")
  
  ggplot(data= coordinates.proteins.hdbscan, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labs(x="Axis x", y="Axis y", title="HDBSCAN cluster in High Dim of Proteins", "HDBSCAN Proteins clustering")


#### PIPELINE B ####
  
#ConsensusClusters
  
  em.consensus<- ConsensusClusterPlus(em.norm, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1, title = "Extended Matrix Consensus Cluster", 
                                      clusterAlg = "hc", distance = "pearson", plot = "png")
  
  data.consensus <- ConsensusClusterPlus(data.norm, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1, title = "Data Consensus Cluster", 
                                         clusterAlg = "hc", distance = "pearson", plot = "png")
  
  #Selecting a Cluster or an Item from em
  icl.em <- calcICL(em.consensus, title = "Extended Matrix Consensus Cluster", plot = "png" )
  icl.em[["clusterConsensus"]]  
  
  most.consensusItem.em <-icl.em[["itemConsensus"]][which(icl.em[["itemConsensus"]][,ncol(icl.em[["itemConsensus"]])]>0.99),]
  
  #Selecting a Cluster or an Item from Data
  icl.data <- calcICL(data.consensus, title = "Data Consensus Cluster", plot = "png" )
  icl.data[["clusterConsensus"]]  
  
  most.consensusItem.data<- icl.data[["itemConsensus"]][which(icl.data[["itemConsensus"]][,ncol(icl.data[["itemConsensus"]])]>0.8),]

  
  
#Plots of Cluster Consensus
  
  ggplot(data = icl.data[["clusterConsensus"]], aes(x = as.factor(k), y = clusterConsensus, fill = factor(cluster))) + geom_bar(position = "dodge", stat = "identity") + theme_light() +
    labs(x = "Number of clusters", y = "Cluster Consensus", title = "Cluster Consensus of Data Matrix", legend = "Cluster")
  
  ggplot(data = icl.em[["clusterConsensus"]], aes(x = as.factor(k), y = clusterConsensus, fill = factor(cluster))) + geom_bar(position = "dodge", stat = "identity") + theme_light() +
    labs(x = "Number of clusters", y = "Cluster Consensus", title = "Cluster Consensus of Extended Matrix", legend = "Cluster")
  
  
#Plots of Item Consensus
  
  ggplot(data = most.consensusItem.data, aes(x = as.factor(k), y = itemConsensus, fill = factor(item))) + geom_bar(position = "dodge", stat = "identity") + theme_light() +
    labs(x = "Number of clusters", y = "Item Consensus", title = "Item Consensus >0.99 of Data Matrix", legend = "Item")
  
  ggplot(data = most.consensusItem.em, aes(x = as.factor(k), y = itemConsensus, fill = factor(item))) + geom_bar(position = "dodge", stat = "identity") + theme_light() +
    labs(x = "Number of clusters", y = "Item Consensus", title = "Item Consensus >0.99 of Extended Matrix", legend = "Item")
  