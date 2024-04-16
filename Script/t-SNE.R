#### UMAP Script ####

  #This script is the development of code with a SIMPLE pipeline, extracted from the development of the R-Script.R file.

#Libraries
  library(ggplot2) #Visualisation of datasets
  library(tidyverse) #Handling of datasets
  library(SNPRelate) #Main analysis of genetic data
  library(openxlsx) #Excel file management
  library(e1071) #SVM
  library(stats) #kmeans
  library(Rtsne) #Dimension reduction

#### DATA INGESTION

#Working directory
  rm(list = ls())
  getwd()
  path = "Set-your-wd-path-Data"
  setwd(path)

#Reading raw data
  fn_vcf <- ""
  fn_phenotypes <- ""
  fn_covariates <- ""
  fn_proteins <- ""
  fn_GDS <- ""


  snpgdsVCF2GDS(fn_vcf, fn_GDS)
  phenotypes <-unlist(unname( read.csv(fn_phenotypes)))
  covariates <-unlist(unname( read.csv(fn_covariates)))
  proteins <-unlist(unname( read.csv(fn_proteins)))
  
  GDS <- snpgdsOpen(fn_GDS)
  
#LD Eliminating
  snpset <- snpgdsLDpruning(GDS, ld.threshold = 0.2, method = "r")

  snpset.id <- unlist(unname(snpset)) 

  genome.matrix <- snpgdsGetGeno(GDS, snp.id = snpset.id)



#### DATA AGGREGATION

#Data treatment

  #We transform chr column of covariates in matrix

  hot_encoding <- model.matrix(~, data = covariates, contrasts = "contr.treatment")

  hot_encoding <- cbind(sample = covariates[,1],hot_encoding)

#PCA

  pca <- snpgdsPCA(GDS, snp.id = snpset.id, num.thread = 10)
  pca.matrix <- cbind(sample = pca$sample.id, pca$eigenvect)

  hot_encoding <- merge(hot_encoding, pca.matrix, by = "sample")

# Data grouping

  #We assume that all data collect samples in rows and properties in columns.

  genome.matrix <- cbind(sample = pca$sample.id, genome.matrix) #We suppose they have the same order, i need a better method (pca, maybe)
  data.matrix <- merge(genome.matrix, hot_encodings, proteins, by="sample")
  label.matrix <- phenotypes
  

#Dim reduction (t-SNE)

  tsne<- Rtsne(data.matrix)

  coordinates.tsne <- data.frame(tsne$Y, labels = label.matrix)

  ggplot(data = coordinates.tsne, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labels() + xlabel() + ylabel()

#### DATA INTERPRETATION

# Clusters K-means
  k_valores <- c()
  kmeans_PCA <- list()
  
  set.seed(123)
  
  for (k in k_valores){ #The loop gives the freedom to make various groupings depending on the number of groups we want.
    agrupamiento <- kmeans(data.matrix[,-1], centers = k)
    k<-as.character(k)
    kmeans[[k]] <- agrupamiento  
  }
  
  for (i in kmeans_PCA){
    print(i$size) #Show the size of the groups on the screen
  }
  
#Plotting
  
  coordinates.tsne <- data.frame(tsne$Y, labels = kmeans$"introduce k number"$cluster)
  
  ggplot(data = coordinates.tsne, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labels() + xlabel() + ylabel()
  
  
# Clusters SVM  
  
  #Split data
  index<- 1:nrow(data.matrix)
  N <- trunc(length(index)/3)
  testindex <- sample(index, N)
  testset <- data.matrix[testindex,]
  trainset <- data.matrix[-testindex,]
  trainset$label <- as.factor(label.matrix[-testindex])
  
  #Model
  svm.model <- svm(label~ ., data = trainset[,-1], cost = 100, gamma = 1)
  svm.pred <- predict(svm.model, testset)
  
  conf.matrix <-table(pred = svm.pred, true = label.matrix[testindex]) #compute SVM confusion matrix
  conf.matrix
  classAgreement(conf.matrix)
  
  testset$cluster <- svm.pred
  testset$label <- label.matrix[testindex]
  
  par(mfrow = c(1,2))
  
  ggplot(testset, aes(EV2, EV1, color = cluster)) + geom_point()
  ggplot(testset, aes(EV2, EV1, color = label)) + geom_point()
