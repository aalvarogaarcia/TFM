#Script to process high throughput proteomics from ILD patients

# Libraries
library(tidyverse)
library(moments)                #Kurtosis &Skewness
library(dummy)
library(ggplot2)
library(umap)
library(Rtsne)
library(dbscan)
library(ConsensusClusterPlus)
library(data.table)
library(Cairo)       #Save plots
library(gridExtra)
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(MASS)

set.seed(123)

#### DATA INGESTION ####

#Define base path
getwd()
path <- "/mnt/756b4178-6cf1-4cab-bdfb-8342cc4c885a/jlorsal/ANALISIS_275GB/Alvaro/06_TFM"
setwd(path)

#File reading
path <- "data"
fn_phenotypes<- "phenotypes.tsv"
fn_proteins <- "Olink.NPX_values.Inflammatory-panels.csv"
phenotypes <- read_tsv(paste(path,fn_phenotypes, sep="/"))
proteins <- read.csv(paste(path,fn_proteins,sep="/"))
#head(phenotypes)  
#head(proteins)  

#### DATA AGGREGATION ####
  
#EXPLORATORY DATA ANALYSIS#
  
#UNIVARIANT ANALYSES

#Search NA's
summary(phenotypes) #Phenotypes  

#Iterate over the proteins and find 'NA'
for (i in 1:ncol(proteins)){
  NA.index <- which(is.na(proteins[,i]))
  if(length(NA.index) != 0){
    print(paste('The protein', colnames(proteins)[i], 'has a missing value in index', NA.index, sep = ' '))
  }
}

#'NA' correction or elimination? Consider imputation methods as alternative
#To be done

#Plot phenotypes
#CairoPNG(filename = "Fig.proteomics.EDA1.png", width = 800, height = 500, pointsize = 12)
#CairoPDF(file = "Fig.proteomics.EDA1.pdf", width = 8, height = 5)
#temp <- phenotypes[, -which(colnames(phenotypes)=='Age')] #Age is the only phenotype with a non-discrete distribution.
#pl <- list()
#To be checked
#for(i in 2:ncol(temp)){
#  temp.table <- as.data.frame(table(temp[,i])) #Table of phenotype
#  #temp[,i] <- as.factor(temp[,i])
#  #print(tempt) #Obtain table
#  #print(prop.table(tempt)) #Obtain prop table
#  pl[[(i-1)]] <- ggplot(temp.table, aes(x = temp.table[,1], y = Freq)) + 
#    geom_bar(stat = "identity") + 
#    xlab('Levels') + 
#    ylab('Quantity')
#}
#grid.arrange(pl[[1]], pl[[2]], pl[[3]], pl[[4]], ncol = 2)

#Straty age
phenotypes$AgeGroup <- cut(phenotypes$Age, seq(0, 100, 10))
phenotypes$Sex <- as.character(phenotypes$Sex)
#Plot
CairoPNG(filename = "Fig.proteomics.EDA1.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.EDA1.pdf", width = 8, height = 5)
p1 <- ggplot(phenotypes, aes(x = Site)) + geom_bar(fill = c("blue", "red")) + xlab('Site') + ylab('Quantity') + theme_light()
p2 <- ggplot(phenotypes, aes(x = Sex)) + geom_bar(fill = c("#0080FF", "#80FFFF")) + xlab('Sex') + ylab('Quantity') + theme_light()
p3 <- ggplot(phenotypes, aes(x = Ancestry)) + geom_bar(fill = c("#80FFFF", "#008000", "#0080FF", "#800080", "#C0C0C0")) + xlab('Ancestry') + ylab('Quantity') + theme_light()
p4 <- ggplot(phenotypes, aes(x = Diagnosis)) + geom_bar(fill = c("#158752", "#8000FF", "#FF8000", "#CCCC56")) + xlab('Diagnosis') + ylab('Quantity') + theme_light()
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

#Plot
CairoPNG(filename = "Fig.proteomics.EDA2.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.EDA2.pdf", width = 8, height = 5)
p5 <- ggplot(phenotypes, aes(x = Age)) + geom_bar() + xlab('Age') + ylab('Quantity') + theme_light()
p6 <- ggplot(phenotypes, aes(x = AgeGroup)) + geom_bar() + xlab('Age group') + ylab('Quantity') + theme_light()
grid.arrange(p5, p6, ncol = 2)
dev.off()

#Distribution analysis
index <- sample(2:ncol(proteins), 16, replace = FALSE)
#Random indexes to plot 16 proteins
temp <- cbind(proteins[,index])
#dim(temp)
#par(mfrow=c(4,4))
#Density (sample)
CairoPNG(filename = "Fig.proteomics.EDA3.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.EDA3.pdf", width = 8, height = 5)
par(mfrow=c(2,2))
for (x in 1:4){
  plot(density(temp[,x]), main = paste("Density of protein number", colnames(temp)[x], sep = " "))
}
dev.off()
  
#Quartiles, Mean, MAD, CV & Outliers (sample)
CairoPNG(filename = "Fig.proteomics.EDA4.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.EDA4.pdf", width = 8, height = 5)
par(mfrow=c(2,2))
for (i in 1:4){ 
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

  boxplot(temp[,i], main = paste("Boxplot of protein number", colnames(temp)[i], sep = " "))
} #Outliers are not discarded, just marked
dev.off()  
  
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

#Show high order moments
N <- paste(kurtosis[[1]], sep = ", ")
L <- paste(kurtosis[[2]], sep = ", ")
P <- paste(kurtosis[[3]], sep = ", ")
cat('Kurtosis lengths:', 'Normal:', length(kurtosis[[1]]), 'Leptokurtic:', length(kurtosis[[2]]), 'Platykurtic:', length(kurtosis[[3]]), sep = '\n')
cat('Kurtosis results:', 'Normal: ', N, 'Leptokurtic: ', L, 'Platykurtic: ', P)
#Kurtosis lengths:
#Normal:
#  0
#Leptokurtic:
#  329
#Platykurtic:
#  39

#Show high order moments
N <- paste(skewness[[1]], sep = ", ")
L <- paste(skewness[[2]], sep = ", ")
P <- paste(skewness[[3]], sep = ", ")
cat('Skewness lengths:', 'Symmetric:', length(skewness[[1]]), 'Right:', length(skewness[[2]]), 'Left:', length(skewness[[3]]), sep = '\n')
cat('Skewness results:', 'Symmetric: ', N, 'Rigth: ', L, 'Left: ', P)
#Skewness lengths:
#Symmetric:
#  0
#Right:
#  368
#Left:
#  0

#BIVARIATE ANALYSES

#Proteins Covariance <- CHECK whether we include this in the final script
temp <- cov(proteins[,-1])
dim(temp)
#[1] 368 368

#Define a threshold for strong correlation
threshold1 <- 0.95
cr <- cor(proteins[,-1])
temp1 <- c()
for(i in 1:ncol(cr)){
  for(j in i:nrow(cr)){
    if(abs(cr[i,j]) > threshold1 && i!=j){
      temp1 <- rbind(temp1, c(i,j))
    }
  }
}
#Define a threshold for weak correlation
threshold2 <- 0.01
cr <- cor(proteins[,-1])
temp2 <- c()
for(i in 1:ncol(cr)){
  for(j in i:nrow(cr)){
    if(abs(cr[i,j]) < threshold2 && i!=j){
      temp2 <- rbind(temp2, c(i,j))
    }
  }
}

#Show the number of proteins with high/low correlation for the selected thresholds
dim(temp1)
dim(temp2)
temp1.cor <- temp1
temp2.cor <- temp2

#Linear regression selecting 9 proteins highly and positively correlated
temp <- sample(1:nrow(temp1), 9)
#'temp1' contains the protein indexs
temp1 <- temp1[temp,]
pl <- list()
for(i in 1:nrow(temp1)){
  data <- data.frame("x" =proteins[,temp1[i,1]+1], "y"=proteins[,temp1[i,2]+1])
  pl[[i]] <- ggplot(aes(x, y), data = data) +
    geom_point()+
    geom_smooth(method = "lm", se = TRUE)+
    xlab(paste('Protein',colnames(proteins)[temp1[i,1]+1]))+
    ylab(paste('Protein',colnames(proteins)[temp1[i,2]+1]))+
    theme_light()
}
CairoPNG(filename = "Fig.proteomics.EDA5.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.EDA5.pdf", width = 8, height = 5)
grid.arrange(pl[[1]], pl[[2]], pl[[3]],
             pl[[4]], pl[[5]], pl[[6]],
             pl[[7]], pl[[8]], pl[[9]],
             ncol = 3)
dev.off()

#Linear regression selecting 9 proteins highly correlated
temp <- sample(1:nrow(temp2), 9)
#'temp2' contains the protein indexs
temp2 <- temp2[temp,]
pl <- list()
for(i in 1:nrow(temp2)){
  data <- data.frame("x" =proteins[,temp2[i,1]+1], "y"=proteins[,temp2[i,2]+1])
  pl[[i]] <- ggplot(aes(x, y), data = data) +
    geom_point()+
    geom_smooth(method = "lm", se = TRUE)+
    xlab(paste('Protein',colnames(proteins)[temp2[i,1]+1]))+
    ylab(paste('Protein',colnames(proteins)[temp2[i,2]+1]))+
    theme_light()
}
CairoPNG(filename = "Fig.proteomics.EDA6.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.EDA6.pdf", width = 8, height = 5)
grid.arrange(pl[[1]], pl[[2]], pl[[3]],
             pl[[4]], pl[[5]], pl[[6]],
             pl[[7]], pl[[8]], pl[[9]],
             ncol = 3)
dev.off()

#Histogram of correlations
correlations <- data.frame(x = cr[lower.tri(cr)])
CairoPNG(filename = "Fig.proteomics.EDA7.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.EDA7.pdf", width = 8, height = 5)
ggplot(correlations, aes(x = x)) + 
  geom_histogram(binwidth = 0.001, fill = "skyblue") + 
  labs(x = "Value", y = "Frequency") + 
  theme_minimal()
dev.off()

#Data merging
#Proteins
labels <- subset(phenotypes, select = c(id, Diagnosis))
head(labels)
dim(proteins)
proteins.nb <- proteins[,-temp1.cor[,1]]
dim(proteins.nb)
rownames(proteins.nb) <- proteins.nb$id

#Hot-encoding of Phenotypes
phenotypes$Diagnosis <- as.factor(phenotypes$Diagnosis)
phenotypes$Site <- as.factor(phenotypes$Site)
phenotypes$Sex <- as.factor(phenotypes$Sex)
temp <- dummy(phenotypes[,-1], int = TRUE)
temp1 <- cbind('Age' = phenotypes$Age, temp)
phenotypes.int <- cbind(phenotypes[,1], temp1)

#Create an extended matrix 'em' merging proteins (without high stepwise-correlated proteins) and phenotypes
em <- merge(proteins.nb, phenotypes.int, by='id')
dim(em)

#Scale proteins and em matrix
proteins.nb <- scale(proteins.nb[,-which(colnames(proteins.nb) == "id")])
em <- scale(em[,-1])
#Extract phenoype from em matrix
data <- em[,-grep('Diagnosis',colnames(em))]
rownames(data) <- rownames(em)

#Save intermediate data
write.csv2(proteins.nb, "proteins.nb.csv", quote = FALSE, row.names = FALSE)
write_tsv(phenotypes.int, "phenotypes.int.tsv", quote = "none")

###DIMENSIONAL REDUCTION

# A) UMAP on the three matrix
em.umap <- umap(em)
data.umap <- umap(data)
proteins.umap <- umap(proteins.nb)

#Plot results
CairoPNG(filename = "Fig.proteomics.DimRed.UMAP.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.DimRed.UMAP.pdf", width = 8, height = 5)
par(mfrow=c(1,3))
#Proteins+covariates+phenotype
eqscplot(em.umap$layout, main = "UMAP Extended Matrix", 
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
#Proteins+covariates
eqscplot(data.umap$layout, main = "UMAP Proteins and Covariates", 
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
#Proteins
eqscplot(proteins.umap$layout, main = "UMAP Proteins", 
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
dev.off()

# B) tSNE
em.tsne <- Rtsne(em)
data.tsne <- Rtsne(data)
proteins.tsne <- Rtsne(proteins.nb)

#Plot results
CairoPNG(filename = "Fig.proteomics.DimRed.tSNE.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.DimRed.tSNE.pdf", width = 8, height = 5)
par(mfrow=c(1,3))
eqscplot(em.tsne$Y, main = "t-SNE Extended Matrix", 
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "t-SNE1", 
         ylab = "t-SNE2")
eqscplot(data.tsne$Y, main = "t-SNT Proteins and Covariates", 
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "t-SNE1", 
         ylab = "t-SNE2")
eqscplot(proteins.tsne$Y, main = "t-SNE Proteins", 
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "t-SNE1", 
         ylab = "t-SNE2")
dev.off()



###UNSPERVISED CLUSTERING
for k = 6
em.kmeans <- kmeans(em, centers = k)
data.kmeans <- kmeans(data, centers = k)
proteins.kmeans <- kmeans(proteins.nb, centers = k)

#K-means
k <- 2:7
em.kmeans <- list()
data.kmeans <- list()
proteins.kmeans <- list()

set.seed(123)
for (i in k){ 
  cluster.em <- kmeans(em, centers = i)
  cluster.data <- kmeans(data, centers = i)
  cluster.proteins <- kmeans(proteins[,-1], centers = i)
  i <- as.character(i)
  em.kmeans[[i]] <- cluster.em
  data.kmeans[[i]] <- cluster.data 
  proteins.kmeans[[i]] <- cluster.proteins 
}

v <-c(em.kmeans, data.kmeans, proteins.kmeans)

#Tocheck
#for(kmeans in v){
#  print(kmeans)
#  for (i in kmeans){
#    print(i$size) #Show cluster sizes
#  }
#}

#Show data with UMAP and color by K-means 
CairoPNG(filename = "Fig.proteomics.K-means.UMAP.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.K-means.UMAP.pdf", width = 8, height = 5)
par(mfrow=c(1,3))
eqscplot(em.umap$layout, col = em.kmeans[[1]]$cluster, 
         main = "UMAP Extended Matrix (k-means)", 
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
eqscplot(data.umap$layout, col =data.kmeans[[1]]$cluster, 
         main = "UMAP Data (k-means)",
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
eqscplot(proteins.umap$layout, col =proteins.kmeans[[1]]$cluster, 
         main = "UMAP Proteins (k-means)",
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
dev.off()

#Show data with t-SNE and color by K-means
CairoPNG(filename = "Fig.proteomics.K-means.t-SNE.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.proteomics.K-means.t-SNE.pdf", width = 8, height = 5)
par(mfrow=c(1,3))
eqscplot(em.tsne$Y, col =em.kmeans[[1]]$cluster, 
         main = "t-sne Extended Matrix (k-means)",
         pch = 15,
         xlab = "t-SNE1", 
         ylab = "t-SNE2")
eqscplot(data.tsne$Y, col =data.kmeans[[1]]$cluster, 
         main = "t-sne Data (k-means)",
         pch = 15,
         xlab = "t-SNE1", 
         ylab = "t-SNE2")
eqscplot(proteins.tsne$Y, col = proteins.kmeans[[1]]$cluster, 
         main = "t-sne Proteins (k-means)",
         pch = 15,
         xlab = "t-SNE1", 
         ylab = "t-SNE2")
dev.off()

par(mfrow=c(1,2))
eqscplot(em.umap$layout, main = "UMAP Extended Matrix", 
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
eqscplot(em.umap$layout, col = em.kmeans[[1]]$cluster, 
         main = "UMAP Extended Matrix (k-means)", 
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")

par(mfrow=c(1,2))
eqscplot(data.umap$layout, 
         main = "UMAP of Proteins and Covariates",
         sub = "Colored by Diagnosis",
         col = (phenotypes.int$Diagnosis_1 + 1), 
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")
eqscplot(data.umap$layout, col = em.kmeans[[1]]$cluster, 
         main = "UMAP of Proteins and Covariates",
         sub = "Colored by K-means clusters",
         pch = 15,
         xlab = "UMAP1", 
         ylab = "UMAP2")

points(em.kmeans[[1]]$centers[,1:2])

par(mfrow=c(1,2))
eqscplot(data.umap$layout, 
         main = "UMAP of Proteins and Covariates",
         sub = "Colored by Diagnosis",
         col = as.integer(phenotypes$Diagnosis) + 2, 
         pch = 16,
         cex = .75,
         xlab = "UMAP1", 
         ylab = "UMAP2")
eqscplot(data.umap$layout, col = em.kmeans[[1]]$cluster, 
         main = "UMAP of Proteins and Covariates",
         sub = "Colored by K-means clusters",
         pch = 1,
         xlab = "UMAP1", 
         ylab = "UMAP2")








#PCA
pca.proteins <- prcomp(proteins[,-1])
# 'm' contains the matrix of PC

#Phenotype chr management
temp <- dummy(phenotypes[,-1], int = TRUE)
temp2 <- cbind(phenotypes[,-ncol(phenotypes)], temp)[,-1]

#The structure stands, lets use normalized
phenotypes.int <- cbind(phenotypes[,1:ncol(phenotypes)-1], temp)


#Data grouping - Labels, em, Data. (Labels = dependent variable, em = extended matrix, Data = independent variables )
labels <- subset(phenotypes.int, select = c(id, Diagnosis))
labels$Diagnosis <- as.factor(labels$Diagnosis)
  
em <- merge(proteins.nb, phenotypes.int, by="id") 
data <- em[,-which(colnames(em)=="Diagnosis")]

#Normalize
em.norm <-scale(em[,-ncol(proteins)])
data.norm <-scale(data[,-ncol(proteins)])


#### PIPELINE A ####
#STEPS:
# > Step1. Dimension reduction using UMAP

#UMAP 1
#'data' stands for 'input data' (proteins and covaraites) excluded the phenotype
#'em' stands for 'extended matrix': a matrix with the proteins data, covariates and phenotypes
#Split the extended matrix by 'Diagnosis'
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
  