
# Libraries
library(tidyverse)
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
  path = "/Users/alvarogarciamunoz/Documents/GitHub/TFM-2.0"
  setwd(path)

#File reading
  path <- "Data"
  fn_phenotypes<- "phenotypes.tsv"
  fn_proteins <- "Olink.NPX_values.Inflammatory-panels.csv"

  phenotypes <- read_tsv(paste(path,fn_phenotypes, sep="/"))
  proteins <- read.csv(paste(path,fn_proteins,sep="/"))
  
  
#### DATA AGGREGATION ####
  
#Phenotype chr management
  temp <- dummy(phenotypes[,-1], int = TRUE)
  temp1 <- scale(cbind(phenotypes[,-ncol(phenotypes)], temp)[,-1])
  temp2 <- cbind(phenotypes[,-ncol(phenotypes)], temp)[,-1]
  

  
  #Phenotype normalization analysis
  par(mfrow = c(3,3))
  for (x in 2:ncol(temp1)){
    plot(density(temp1[,x]), main = paste("Density of phenotype", colnames(temp1)[x], sep = " "))
  }
  
  par(mfrow = c(3,3))
  for (x in 2:ncol(temp2)){
    plot(density(temp2[,x]), main = paste("Density of phenotype", colnames(temp2)[x], sep = " "))
  }
  
  #The structure stand, lets use normalized
  phenotypes.int<- cbind(phenotypes[,1], temp1)
  
  
  
#Proteins Plotting
  index <- sample(1:ncol(proteins), 16, replace = FALSE)

  temp <- proteins[,index]
  par(mfrow=c(4,4))
  for (x in 1:ncol(temp)){
    plot(density(temp[,x]), main = paste("Density of protein number", colnames(temp)[x], sep = " "))
  }
  
  rm(temp)


#Data grouping - Labels, em, Data. (Labels = dependent variable, em = extended matrix, Data = independent variables )

  labels <- subset(phenotypes.int, select = c(id, Diagnosis))
  labels$Diagnosis <- as.factor(labels$Diagnosis)
  
  em <- merge(proteins, phenotypes.int, by="id") 
  data <- em[,-which(colnames(em)=="Diagnosis")]

#Normalize
  em.norm <-scale(em[,-ncol(proteins)])
  data.norm <-scale(data[,-ncol(proteins)])
  
  
#### EXPLORATORING DATA ANALYSIS ####
  
#Normal
  
  
  
#### PIPELINE A ####
  
#UMAP 1
  em.umap <-umap(em.norm)
  data.umap <-umap(data.norm)
  proteins.umap <- umap(proteins[,-1])
  
  coordinates.em.umap <- data.frame(em.umap$layout, labels = labels$Diagnosis)
  coordinates.data.umap <- data.frame(data.umap$layout, labels = labels$Diagnosis)
  coordinates.proteins.umap <-data.frame(proteins.umap$layout, labels = labels$Diagnosis)
  
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
  em.hdbscan <- hdbscan(em.norm, minPts = 3)
  data.hdbscan <- hdbscan(data.norm, minPts = 3)
  proteins.hdbscan <- hdbscan(proteins[,-1], minPts = 3)
  
  coordinates.em.hdbscan <- data.frame(em.umap$layout, labels = as.factor(em.hdbscan$cluster))
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
  