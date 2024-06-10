#### R - SCRIPT ##### 

#Este documento recoge ordenadamente y con comentarios el código de R utilizado para el informe dinámico

#### INTRODUCCIÓN ####

#install R packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")

#Load R packages
library(ggplot2)     #Visualización de los datasets
library(tidyverse)   #Manejo de los datasets
library(dummy)       #Hot-encoding
library(SNPRelate)   #Librería básica encargada del análisis principal de los datos genéticos
library(reticulate)  #Uso de python
library(RefManageR)  #Bibliografía
library(openxlsx)    #Gestión de archivos tipo excel
library(stats)       #Incluye diferentes estadísticas
library(e1071)       #Supervised Learning
library(glmnet)      #Supervised Learning
library(umap)        #Dimension Reduction
library(Rtsne)       #Dimension Reduction
library(tensorflow)  #Unsupervised Learning
library(Cairo)       #Save plots
library(gridExtra)   #Plots in panels
library(dbscan)      #DBSCAN clustering


#Define base path
getwd()
path <- "/mnt/756b4178-6cf1-4cab-bdfb-8342cc4c885a/jlorsal/ANALISIS_275GB/Alvaro/06_TFM"
setwd(path)

#### DATA MANAGEMENT ####
rm(list = ls())

#Access to file names
fn_vcf <- "data/chr22.CHB.CEU.FIN.IBS.TSI.YRI.715inds.25-35Mb.vcf.gz"  #VCF
fn_cv <- "data/lists.xlsx"                                             #Covariates
fn_red <- "data/list715"                                               #List of individuals
fn_GDS <- "data/chr22"                                                 #GDS file

#Load data
snpgdsVCF2GDS(fn_vcf, fn_GDS)          #CreateGDS file
co.var <- read.xlsx(fn_cv)             #Read covariates
red <- read.csv(fn_red, header= FALSE) #Read list of individuals
red <- unlist(unname(red))             #Unlist individuals

#Explore some of the loaded data
#head(co.var)
#head(red)

#Population labels (taken from '1000 Genomes project', aka '1KGP') used as a filter
filter <- c("CEU", "FIN", "GBR", "IBS", "TSI", "CHB", "YRI")
filter

#Data from 1KGP can contain up to 2,504 individuals from 26 different populations
#Filter 'co.var' using only the listed individuals
dim(co.var.red)
co.var.red <- filter(co.var, co.var$sample %in% red)
co.var.red <- co.var[co.var$sample %in% red & co.var$pop %in% filter, ]
dim(co.var.red)


#### EXPLORATORY ANALYSIS ####
GDS <- snpgdsOpen(fn_GDS)

#Reduction of data set using linkage disequilibrium
snpset <- snpgdsLDpruning(GDS, ld.threshold = 0.2, method = "r")
#1,743 markers are selected in total.
snpset.id <- unlist(unname(snpset)) 
#IDs of remainings SNPs in dataset
head(snpset.id)
#[1]  19  84 217 316 801 804
#Matrix of genotype
genome.matrix <- snpgdsGetGeno(GDS, snp.id = snpset.id)
#Genotype matrix: 338 samples X 1743 SNPs


#### DIMENSION REDUCTION ####

## PCA ##
pca <- snpgdsPCA(GDS, snp.id= snpset.id, num.thread = 10) #PCA GDS file
str(pca)

#PCA as df
df <- data.frame(sample = pca$sample.id,
                 PC1=pca$eigenvect[,1],
                 PC2=pca$eigenvect[,2],
                 PC3=pca$eigenvect[,3],
                 PC4=pca$eigenvect[,4],
                 PC5=pca$eigenvect[,5],
                 PC6=pca$eigenvect[,6],
                 PC7=pca$eigenvect[,7],
                 PC8=pca$eigenvect[,8],
                 PC9=pca$eigenvect[,9],
                 PC10=pca$eigenvect[,10],
                 stringsAsFactors=FALSE)
str(df)

#Variance explained by each principal component (as percentage)
#'pca$eigenval' contains the eigenvalues from the PCA
variance.explained <- ( pca$eigenval / sum(pca$eigenval, na.rm = TRUE) ) * 100
head(variance.explained)
paste("PCA 1-50", round(sum(variance.explained[1:50], na.rm = TRUE), 2), "%")
paste("PCA 1-10", round(sum(variance.explained[1:10], na.rm = TRUE),2), "%")

#Same method but for any matrix
pca <- prcomp(genome.matrix) #PCA any matrix

df <- data.frame(sample = co.var.red[,1], pca$x[,1:10], stringsAsFactors = FALSE)

variance.explained <- ( pca$sdev / sum(pca$sdev) ) * 100

paste("PCA 1-50", round(sum(variance.explained[1:50], na.rm = TRUE), 2), "%")
paste("PCA 1-10", round(sum(variance.explained[1:10], na.rm = TRUE), 2), "%")


#Scree-plot
CairoPNG(filename = "Fig.scree-plots.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.scree-plots.pdf", width = 8, height = 5)
par(mfrow=c(1,2))
screeplot(pca, npcs = 50, type ='lines', 
          main = 'Explained variance', 
          sub = paste('50 PC:', round(sum(variance.explained[1:50], na.rm = TRUE),2), '% Variance'))
screeplot(pca, npcs = 10, type ='lines', 
          main = 'Explained variance', 
          sub = paste('10 PC:', round(sum(variance.explained[1:10], na.rm = TRUE),2), '% Variance'))
dev.off()

## Plot ##

#Coloring individuals from 26 populations of 1KGP
colors <- c("#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF",
             "#00FFFF", "#FF8000", "#80FF00", "#8000FF", "#FFFF80", "#FF80FF",
             "#80FFFF", "#008000", "#0080FF", "#800080", "#808080", "#C0C0C0",
             "#A0A0A0", "#606060", "#404040", "#202020", "#AA12AA", "#158752", 
             "#CCCC56", "#949513") 
colors <- c("#808080", "#FF0000", "#00FFFF", "#0080FF", "#FF8000", "#CCCC56",
             "#008000", "#FF8000") 

#Command for plot
pca_cv <- cbind(df, pop = co.var.red$pop, super_pop = co.var.red$super_pop)
CairoPNG(filename = "Fig.PCA.1KGP.7pops.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.PCA.1KGP.7pops.pdf", width = 8, height = 5)
par(mfrow=c(1,1))
ggplot(pca_cv, aes(PC1, PC2, color = pop, shape = super_pop)) + 
  geom_point(size = 3) + 
  scale_color_manual(values = colors) + 
  theme_light()
dev.off()


#Define phen hot-encoding (General)
df_hot_encoding <- dummy(co.var.red[,-1], int = TRUE)
df_hot_encoding <- cbind(sample = co.var.red[,1], df_hot_encoding)
head(df_hot_encoding)

#Write CSV
pca_cv <- merge(df, df_hot_encoding, by = "sample")
fn_csv <- "PCA10_covariates"
write.csv2(pca_cv, fn_csv, quote = FALSE, row.names = FALSE)


### UMAP ###

# Data and label
label <- NA
df.data <- cbind(df[,1], genome.matrix, df_hot_encoding[,-1])           #df.data contains all information
if(!is.na(label)){df.data <- df.data[,-grep(label, colnames(df.data))]} #Do not use label in the reduction
A <- as.matrix(df.data[,-1])                                            #A will be the matrix of all the samples
v <- co.var.red[,which(colnames(co.var.red) == label)]                  #V will be a vector label


#Algorithm apply

#Set UMAP params
custom.config <- umap.defaults
custom.config$n_neighbors <- 15
custom.config$n_components <- 2
custom.config$metric <- 'euclidean'
custom.config$n_epochs <- 200

#Implementation
umap <- umap(A, custom.config)
umap.comp <- data.frame(umap$layout, labels = co.var.red$pop)

#Plot
CairoPNG(filename = "Fig.UMAP.1KGP.7pops.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.UMAP.1KGP.7pops.pdf", width = 8, height = 5)
par(mfrow=c(1,1))
ggplot(umap.comp, aes(X1, X2, color = labels)) + 
  geom_point() + 
  theme_light() + 
  labs(color = "Populations", x = "UMAP1", y = "UMAP2")
dev.off()

#UMAP learn - Python
df.umap.learn <- umap(df.data, method = "umap-learn") #Needs a python module


### t-SNE ###

#Implementation
df.tsne <- Rtsne(A,
                 dims = 2,
                 perplexity = 30, 
                 max_iter = 1000,
                 check_duplicates = FALSE)

tsne.comp <- data.frame(df.tsne$Y, labels = co.var.red$pop)

#Plot
CairoPNG(filename = "Fig.tSNE.1KGP.7pops.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.tSNE.1KGP.7pops.pdf", width = 8, height = 5)
par(mfrow=c(1,1))
ggplot(tsne.comp, aes(X1, X2, color = labels)) + 
  geom_point() + 
  theme_light() +
  labs(color = "Populations", x = "t-SNE1", y = "t-SNE2")
dev.off()


#### CLUSTERS #####

#### Unsupervised - Learning ####

#K-means
k <- 2:7
kmeans <- list()

set.seed(123)
for (i in k){
  cluster <- kmeans(A, centers = i)
  i <- as.character(i)
  kmeans[[i]] <- cluster  
}

for (i in kmeans){
  print(i$size) #Show cluster sizes
}

#Keep clusters information in a data frame
df_kmeans <- df

#Prepare a list of plots for gridExtra
pl <- list()
for(i in 1:length(kmeans)){
  df_kmeans$cluster <- as.factor(kmeans[[i]]$cluster)
  p <- ggplot(df_kmeans, aes(PC1, PC2, color = cluster)) + 
    geom_point(size = 2) +
    ggtitle(paste("k=", names(kmeans)[i], sep = ""))
    theme_light()
  pl[[i]] <- p
}

#Plots in a grid for K values
CairoPNG(filename = "Fig.K-means.1KGP.7pops.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.K-means.1KGP.7pops.pdf", width = 8, height = 5)
grid.arrange(pl[[1]], pl[[2]], pl[[3]], pl[[4]], pl[[5]], pl[[6]], ncol = 3)
dev.off()


#DBSCAN
dbscan <- dbscan(df[,-1], eps = 4, minPts = 5)
pca_cv$cluster <- as.factor(dbscan$cluster)
umap.comp$cluster <- as.factor(dbscan$cluster)

p<-ggplot(pca_cv, aes(PC1, PC2, color = cluster)) + 
  geom_point() + 
  theme_light() + 
  labs(color = "clusters", x = "PC1", y = "PC2")+
  ggtitle('Epsilon = 4, minPts = 5')
u <- ggplot(umap.comp, aes(X1, X2, color = cluster)) + 
  geom_point() + 
  theme_light() + 
  labs(color = "clusters", x = "UMAP1", y = "UMAP2")+
  ggtitle('Epsilon = 4, minPts = 5')

CairoPNG(filename = "Fig.DBSCAN.1KGP.7pops.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.DBSCAN.1KGP.7pops.pdf", width = 8, height = 5)
grid.arrange(p, u, ncol = 2)
dev.off()

#HDBSCAN
hdbscan <- hdbscan(df[,-1], minPts = 5)
pca_cv$cluster <- as.factor(hdbscan$cluster)
umap.comp$cluster <- as.factor(hdbscan$cluster)

p<-ggplot(pca_cv, aes(PC1, PC2, color = cluster)) + 
  geom_point() + 
  theme_light() + 
  labs(color = "clusters", x = "PC1", y = "PC2")+
  ggtitle('minPts = 5')

u <-ggplot(umap.comp, aes(X1, X2, color = cluster)) + 
  geom_point() + 
  theme_light() + 
  labs(color = "clusters", x = "UMAP1", y = "UMAP2")+
  ggtitle('minPts = 5')

CairoPNG(filename = "Fig.HDBSCAN.1KGP.7pops.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.HDBSCAN.1KGP.7pops.pdf", width = 8, height = 5)
grid.arrange(p, u, ncol=2)
dev.off()

#ConsensusClusterPlus (not used in the final version of the TFM)
#library(ConsensusClusterPlus)
#
#df.data <- merge(df, co.red.dummies, by = "sample")
#df.data <- cbind(genome.matrix, df.data[,-1])
#
#Clustering
#title = tempdir()
# results = ConsensusClusterPlus(df.data, maxK = 26, reps = 50, pItem = 0.8, pFeature = 1, 
#                                 title = title, clusterAlg = "hc", distance = "pearson", plot = "png")
#
#ICL
#icl <- calcICL(results, title = title, plot='png')
#icl[['clustersConsensus']]


#### Supervised - Learning ####
# 
#SVM (not used in the final version of the TFM)
#
#Split data
#index<- 1:nrow(df_hot_encoding)
#N <- trunc(length(index)/3)
#testindex <- sample(index, N)
#testset <- df.data[testindex,]
#trainset <- df.data[-testindex,]
#trainset$super_pop <- as.factor(df_hot_encoding$super_pop[-testindex])
#
#Model
#svm.model <- svm(super_pop~ ., data = trainset[,-1], cost = 100, gamma = 1)
#svm.pred <- predict(svm.model, testset)
#
#conf.matrix <-table(pred = svm.pred, true = df_hot_encoding$super_pop[testindex]) #compute SVM confusion matrix
#conf.matrix
#classAgreement(conf.matrix)
#
#testset$cluster <- svm.pred
#testset$super_pop <- df_hot_encoding$super_pop[testindex]
#
#ggplot(testset, aes(EV2, EV1, color = cluster)) + geom_point()
#ggplot(testset, aes(EV2, EV1, color = super_pop)) + geom_point()
#
#
#GLMNET (not used in the final version of the TFM)
#
#glmnet object
#v.labels <- as.numeric(as.factor(df.labels)) #We transform our labels to a numeric vector
#fit <- glmnet(df.data, v.labels) #As default we use a Gaussian linear model (least squares)
#As default the fit would be a lasso regression (\alpha = 1).
#
#fit shows
#plot(fit)
#print(fit)
#coef(fit)
#
#cv.glmnet object
#df.matrix <- as.matrix(df.data) #We transform our df to a matrix
#cvfit <- cv.glmnet(df.matrix, v.labels) #Cross validation
#
#cv.fit shows
#plot(cvfit)
#cvfit$lambda.min
#coef(cvfit)
#predict(cvfit, newx = df.matrix[1:30,])
#
# Add penalty factors
#p.fac <- rep(1,ncol(df.matrix))
#p.fac[c(11, 17, 10)]<- 0
#pfit <- glmnet(df.matrix, v.labels, penalty.factor = p.fac)
#plot(pfit, label=TRUE)
#
#Confusion matrices
#cfit <- cv.glmnet(df.matrix[-testindex,], v.labels[-testindex], family="multinomial")
#cnf <- confusion.glmnet(cfit, newx = df.matrix[testindex,], newy = v.labels[testindex])
#print(cnf)
#
#ROC *glmnet* (Use only with binomial family)
#bfit <- cv.glmnet(df.matrix, v.labels, family="binomial", type.measure = "auc", keep=TRUE)
#roc <- roc.glmnet(bfit$fit.preval, newy = v.labels)
#
#
#### Deep - Learning with Keras and Tensorflow ####  (not used in the final version of the TFM)
#
# Preparing data
#x.matrix.train <- as.matrix(df.data[-testindex, ])
#x.matrix.test <- as.matrix(df.data[testindex, ])
#
#y.train <- hot_encoding_super_pop[-testindex, -1]
#y.test <- hot_encoding_super_pop[testindex, -1]
#
#write.csv(as.matrix(df.data), "chr22-matrix")
#write.csv(v.labels, "chr22-labels")
#
# Defining model (sequential)
#    
#    model <- keras_model_sequential() 
#    model %>%  #Defining the model for our study
#      layer_dense(units = 256, activation = 'relu', input_shape = c(23,17)) %>% 
#      layer_dropout(rate = 0.4) %>% 
#      layer_dense(units = 128, activation = 'relu') %>%
#      layer_dropout(rate = 0.3) %>%
#      layer_dense(units = 10, activation = 'softmax')
#
#    model %>% compile(
#      loss = 'categorical_crossentropy',
#      optimizer = optimizer_rmsprop(),
#      metrics = c('accuracy')
#    )
#
##End Of Script