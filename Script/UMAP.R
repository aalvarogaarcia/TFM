#### UMAP Script ####

#This script is the development of code with a SIMPLE pipeline, extracted from the development of the R-Script.R file.

#Libraries
library(ggplot2) #Visualisation of datasets
library(tidyverse) #Handling of datasets
library(SNPRelate) #Main analysis of genetic data
library(openxlsx) #Excel file management
library(dbscan) #Hierarchical DB scan
library(umap) #Dimension reduction

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


#Dim reduction (UMAP)

  umap<- umap(data.matrix)

  coordinates.umap <- data.frame(umap$layout, labels = label.matrix)

  ggplot(data = coordinates.umap, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labels() + xlabel() + ylabel()

#### DATA INTERPRETATION

#HDBSCAN Clustering

  clusters <- hdbscan(data.matrix, minPts = 5)

  coordinates.hdbscan <- data.frame(umap$layout, labels = clusters)

  ggplot(data = coordinates.umap, aes(X1, X2, color = labels)) + geom_point() + theme_light() +
    labels() + xlabel() + ylabel()













