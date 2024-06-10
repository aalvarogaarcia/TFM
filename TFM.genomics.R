#Script to process genotype from COVID-19 patients

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
library(ConsensusClusterPlus)
library(MASS)

rm(list = ls())

set.seed(123)

#### DATA INGESTION ####

#Define base path
getwd()
path <- "C:/Users/agarm/Desktop/TFM 2.0/data"
setwd(path)

#File reading
path <- "Data"
fn_phenotypes<- "corticoids_n1487.tsv"
fn_pc <- "corticoids_n1487.geno0.01.mind0.01.n1500.noHighLD.recoded.raw"
fn_ev <- "corticoids_n1487.eigenval"
phenotypes <- read_tsv(paste(path,fn_phenotypes, sep="/"))
ev <- read.csv(paste(path,fn_ev, sep="/"), header = FALSE)
#genotype <- read.ftable(paste(path,fn_pc,sep="/"), row.var.names = FALSE, col.vars = )

#head(phenotypes)  
#head(pc)  

#Colors palette
colors <- c("#00FFFF", "#242278", "#00A1F2", "#00EEFF", "#111691",  #Because of my TFM I've choose a blue gradient
            "#00B7FF", "#00DDFF", "#0300C2", "#00C3FF", "#0042E8")  #From https://coolors.co (Colors palettes IA)


#### DATA AGGREGATION ####

#EXPLORATORY DATA ANALYSIS#

#UNIVARIANT ANALYSES

#Table of data
round(prop.table(table(phenotypes$sex))*100,1)
round(prop.table(table(phenotypes$center))*100,2)
round(prop.table(table(phenotypes$mort))*100,1)

#Search NA's
summary(phenotypes) #Phenotypes  

#Delete NA's
index <- which(phenotypes$mort == -9)
phenotypes <- phenotypes[-index,]


#Group Ages
phenotypes$agegroup <- cut(phenotypes$age, seq(0, 110, 10))


#Plot Covariates
p1 <- ggplot(phenotypes, aes(x = as.factor(sex))) + geom_bar(fill = colors[3:4]) + xlab('Sex') + ylab('Quantity') + theme_light()
p2 <- ggplot(phenotypes, aes(x = as.factor(center))) + geom_bar(fill = colors[1:9]) + xlab('Center') + ylab('Quantity') + theme_light()
p3 <- ggplot(phenotypes, aes(x = as.factor(mort))) + geom_bar(fill = colors[9:10]) + xlab ('Mort') + ylab('Quantity') + theme_light()


CairoPNG(filename = "Fig.genomic.EDA1.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomic.EDA1.pdf", width = 8, height = 5)
grid.arrange(p1, p2, p3, ncol = 2)
dev.off()

p5 <- ggplot(phenotypes, aes(as.integer(age))) + geom_bar(fill = colors[1]) + xlab('Age') + ylab('Quantity') + theme_light()
p6 <- ggplot(phenotypes, aes(x = agegroup)) + geom_bar(fill=colors[2]) + xlab('Age group') + ylab('Quantity') + theme_light()

CairoPNG(filename = "Fig.genomic.EDA2.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomic.EDA2.pdf", width = 8, height = 5)
grid.arrange(p5, p6, ncol=2)
dev.off()


#Principal components Analysis

pl <- list()
i<-1
for (x in phenotypes[,2:11]){
  
  pl[[i]] <- ggplot(phenotypes, aes(x)) +
    geom_density(color = colors[2], fill = colors[1]) +
    labs(caption = paste('Density of PC', i), x = 'PC', y = 'Densidad') +
    theme_light()
  i <- i+1
}

CairoPNG(filename = "Fig.genomic.EDA3.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomic.EDA3.pdf", width = 8, height = 5)
grid.arrange(pl[[1]],pl[[2]],pl[[3]],
             pl[[4]],pl[[5]],pl[[6]],
             pl[[7]],pl[[8]],pl[[9]],pl[[10]], ncol = 5)
dev.off()

#Kurtosis & Skewness (all)
kurtosis <- list('Normal' = c(), 'Leptokurtic' = c(), 'Platykurtic' = c())
skewness <- list('Symmetric' = c(), 'Right' = c(), 'Left' = c())
for (i in 2:11){
  temp <- kurtosis(phenotypes[,i])
  tempt <- skewness(phenotypes[,i])
  if(temp == 3){
    kurtosis[[1]] <- cbind(kurtosis[[1]], colnames(phenotypes)[i])
  }
  else if(temp >3){
    kurtosis[[2]] <- cbind(kurtosis[[2]], colnames(phenotypes)[i])
  }
  else{
    kurtosis[[3]] <- cbind(kurtosis[[3]], colnames(phenotypes)[i])
  }
  if(temp == 0){
    skewness[[1]] <- cbind(skewness[[1]], colnames(phenotypes)[i])
  }
  else if(temp >0){
    skewness[[2]] <- cbind(skewness[[2]], colnames(phenotypes)[i])
  }
  else{
    skewness[[3]] <- cbind(skewness[[3]], colnames(phenotypes)[i])
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
#  10
#Platykurtic:
#  0

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
#  10
#Left:
#  0


#Scree-Plot
var <- ev/sum(ev)
EV <-data.frame(PC = 1:10, V1= var)
CairoPNG(filename = "Fig.genomics.scree-plots.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.scree-plots.pdf", width = 8, height = 5)
ggplot(EV, aes(PC, V1)) + 
  geom_line(col = colors[2]) +
  xlab("PC") +
  ylab("Variance") +
  labs(title ="Scree-plot", caption = " COVID-19") +
  theme_light() +
  ylim(0,0.25)
dev.off()



#BIVARIATE ANALYSIS

CairoPNG(filename = "Fig.genomic.EDA4.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomic.EDA4.pdf", width = 8, height = 5)
ggplot(phenotypes, aes(x = pc1, y = age, col = as.factor(sex))) + 
  geom_point() +
  geom_smooth(method = 'lm', se =FALSE) +
  geom_smooth(aes(pc1, age), method = 'lm',inherit.aes = FALSE, col = colors[10], linetype = "dashed") +
  theme_light()
dev.off()




#MATRIX PREPARATION

#PC
pc <- as.matrix(phenotypes[,2:11])
str(pc)

#Hot-encoding of covariates

covariates <- phenotypes[,12:(ncol(phenotypes)-2)]
covariates$sex <- as.factor(covariates$sex)
covariates$center <- as.factor(covariates$center)

covariates.int <- dummy(covariates, int = TRUE)
covariates.int <- cbind(covariates.int, age = phenotypes$age)

str(covariates.int)

#Phenotype Vector
phenotype.v <-as.factor(phenotypes$mort)


#Merging to data matrix
data <- scale(cbind(pc, covariates.int))

str(data)

#MULTIVARIATE ANALYSIS

#Robust Mean point
mu <- c()
for (i in 2:ncol(pc)){
  temp <- rlm(pc[,i]~1, maxit = 100)$coefficients
  mu <- cbind(mu ,temp)
}

#Robust covariance matrix
sig <- cov.rob(pc)$cov

#Distance to mu (Mahalanobis)
dist.mu <- apply(pc, 1, function(x) mahalanobis(x, mu, sig))
dist.mu <- data.frame(Index = 1:length(dist.mu), distance = dist.mu)

#Plot
CairoPNG(filename = "Fig.genomics.Outlier.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.Outlier.pdf", width = 8, height = 5)
ggplot(dist.mu, aes(Index, distance)) + 
  geom_point(col = colors[2]) + 
  theme_light()+
  ggtitle('Distance Mahalanobis')
dev.off()

#Outlier detection
threshold3 <- 3 * sd(dist.mu$distance)
outlier <- which(dist.mu$distance > threshold3)

fn_outlier <- "outlier.phen-genomic.csv"
write.csv2(phenotypes[outlier,], fn_outlier)

#Outlier factor
outlier.v <- rep(0, nrow(pc))
outlier.v[outlier] <- 1
outlier.v <- factor(outlier.v)
levels(outlier.v) <- c('Normal', 'Outlier')


###DIMENSIONAL REDUCTION

#UMAP
set.seed(123)
data <- scale(data)
pc <- scale(pc)
data.umap <- umap(data)
pc.umap <- umap(pc)

#Plot
data.umap <-data.frame(data.umap$layout, mort = phenotype.v, outlier = outlier.v) 
pc.umap <- data.frame(pc.umap$layout, mort = phenotype.v , outlier = outlier.v)

mort.data.umap <- ggplot(data.umap, aes(X1, X2, col = mort, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'Mort color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) + 
  scale_color_manual(values = colors[1:2])



mort.pc.umap <- ggplot(pc.umap, aes(X1, X2, col = mort, shape = outlier )) + 
  geom_point() + 
  labs(title='PC',
       caption= 'mort color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) +
  scale_color_manual(values = colors[1:2])

CairoPNG(filename = "Fig.genomics.DimRedUMAP.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.DimRedUMAP.pdf", width = 8, height = 5)
grid.arrange(mort.data.umap, mort.pc.umap, ncol = 2)
dev.off()


#t-SNE
set.seed(123)

data.tsne <- Rtsne(data)
pc.tsne <- Rtsne(pc)

#I set up two data frame w/ the t-SNE reduction and the factor mort (if he has diagosis 1 or other)
data.tsne <-data.frame(X1 = data.tsne$Y[,1], X2 = data.tsne$Y[,2], mort = phenotype.v, outlier = outlier.v) 
pc.tsne <- data.frame(X1 = pc.tsne$Y[,1], X2 = pc.tsne$Y[,2], mort = phenotype.v , outlier = outlier.v)


mort.data.tsne <- ggplot(data.tsne, aes(X1, X2, col = mort, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'mort color') +
  xlab('TSNE1') + ylab('TSNE2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8))+
  scale_color_manual(values = colors[1:2])


mort.pc.tsne <- ggplot(pc.tsne, aes(X1, X2, col = mort, shape = outlier )) + 
  geom_point() + 
  labs(title='pc',
       caption= 'mort color') +
  xlab('TSNE1') + ylab('TSNE2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8))+
 scale_color_manual(values = colors[1:2])



CairoPNG(filename = "Fig.genomics.DimRedTSNE.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.DimRedTSNE.pdf", width = 8, height = 5)
grid.arrange(mort.data.tsne, mort.pc.tsne, ncol = 2)
dev.off()

#MDS

set.seed(123)
p <- 1
data.dist <- dist(data, method = 'minkowski', p = p)
pc.dist <- dist(pc, method = 'minkowski', p = p)

data.mds <-data.frame(X1 = cmdscale(data.dist)[,1], X2 = cmdscale(data.dist)[,2], mort = phenotype.v, outlier = outlier.v) 
pc.mds <- data.frame(X1 = cmdscale(pc.dist)[,1], X2 = cmdscale(pc.dist)[,2], mort = phenotype.v , outlier = outlier.v)

mort.data.mds <- ggplot(data.mds, aes(X1, X2, col = mort, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'mort color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) + 
  scale_color_manual(values = colors[1:2])


mort.pc.mds <- ggplot(pc.mds, aes(X1, X2, col = mort, shape = outlier )) + 
  geom_point() + 
  labs(title='pc',
       caption= 'mort color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) + 
  scale_color_manual(values = colors[1:2])

CairoPNG(filename = "Fig.genomics.DimRedMDS.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.DimRedMDS.pdf", width = 8, height = 5)
grid.arrange(mort.data.mds, mort.pc.mds, ncol = 2)
dev.off()



#UNSUPERVISED CLUSTERING

k <- 2:7
data.kmeans <- list()
pc.kmeans <- list()

set.seed(123)
for (i in k){ 
  cluster.data <- kmeans(data, centers = i)
  cluster.pc <- kmeans(pc, centers = i)
  i <- as.character(i)
  data.kmeans[[i]] <- cluster.data 
  pc.kmeans[[i]] <- cluster.pc 
}

v <-list(data.kmeans, pc.kmeans)

#Plot K-mean UMAP
set.seed(123)
data.umap <-data.frame(umap(data)$layout, mort = phenotype.v, outlier = outlier.v) 
pc.umap <- data.frame(umap(pc)$layout, mort = phenotype.v , outlier = outlier.v)

data.umap <- cbind(data.umap, cluster.km = as.factor(v[[1]][[1]]$cluster))
pc.umap <- cbind(pc.umap, cluster.km = as.factor(v[[2]][[1]]$cluster))

km.data.umap <- ggplot(data.umap, aes(X1, X2, col = cluster.km, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'k-means color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) + 
  scale_color_manual(values = colors[1:2])

km.data.umap


km.pc.umap <- ggplot(pc.umap, aes(X1, X2, col = cluster.km, shape = outlier )) + 
  geom_point() + 
  labs(title='PC',
       caption= 'k-means color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) +
  scale_color_manual(values = colors[1:2])


km.pc.umap

CairoPNG(filename = "Fig.genomics.ClKM1.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.ClKM1.pdf", width = 8, height = 5)
grid.arrange(km.data.umap, km.pc.umap, ncol = 2)
dev.off()

#Plot tSNE
data.tsne <- cbind(data.tsne,  cluster.km = as.factor(v[[1]][[1]]$cluster))

km.data.tsne <- ggplot(data.tsne, aes(X1, X2, col = cluster.km, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'k-means color') +
  xlab('TSNE1') + ylab('TSNE2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) + 
  scale_color_manual(values = colors[1:2])

km.data.tsne

CairoPNG(filename = "Fig.genomics.ClKM2.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.ClKM2.pdf", width = 8, height = 5)
km.data.tsne
dev.off()

#Plot MDS
data.mds <-data.frame(X1 = cmdscale(data.dist)[,1], X2 = cmdscale(data.dist)[,2], mort = phenotype.v, outlier = outlier.v) 
pc.mds <- data.frame(X1 = cmdscale(pc.dist)[,1], X2 = cmdscale(pc.dist)[,2], mort = phenotype.v , outlier = outlier.v)

data.mds <- cbind(data.mds, cluster.km = as.factor(v[[1]][[1]]$cluster))
pc.mds <- cbind(pc.mds, cluster.km = as.factor(v[[2]][[1]]$cluster))

km.data.mds <- ggplot(data.mds, aes(X1, X2, col = cluster.km, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'k-means color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) + 
  scale_color_manual(values = colors[1:2])

km.data.mds


km.pc.mds <- ggplot(pc.mds, aes(X1, X2, col = cluster.km, shape = outlier )) + 
  geom_point() + 
  labs(title='pc',
       caption= 'k-means color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) +
  scale_color_manual(values = colors[1:2])

km.pc.mds

CairoPNG(filename = "Fig.genomics.ClKM3.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.ClKM3.pdf", width = 8, height = 5)
grid.arrange(km.data.mds, km.pc.mds, ncol = 2)
dev.off()


#DBSCAN
set.seed(123)
data.dbscan <- as.factor(dbscan(data, eps = 5)$cluster)
pc.dbscan <- as.factor(dbscan(pc, eps = 0.1)$cluster)

levels(data.dbscan) <- c('Noise', 'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5')
levels(pc.dbscan)<-  c('Noise', 'Cluster 1')

set.seed(123)
data.umap <-data.frame(umap(data)$layout, mort = phenotype.v, outlier = outlier.v) 
pc.umap <- data.frame(umap(pc)$layout, mort = phenotype.v , outlier = outlier.v)

data.umap <- cbind(data.umap, cluster.db = data.dbscan)
pc.umap <- cbind(pc.umap, cluster.db = pc.dbscan)

db.data.umap <- ggplot(data.umap, aes(X1, X2, col = cluster.db, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'DBSCAN (eps= 5; minPts= 5) color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) 

db.pc.umap <- ggplot(pc.umap, aes(X1, X2, col = cluster.db, shape = outlier )) + 
  geom_point() + 
  labs(title='PC',
       caption= 'DBSCAN (eps= 0.1; minPts= 5) color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) 

#Plot tSNE
data.tsne <- cbind(data.tsne,  cluster.db = data.dbscan)

db.data.tsne <- ggplot(data.tsne, aes(X1, X2, col = cluster.db, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'DBSCAN (eps= 5; minPts= 5) color') +
  xlab('TSNE1') + ylab('TSNE2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) 

db.data.tsne

CairoPNG(filename = "Fig.genomics.Cldb1.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.Cldb1.pdf", width = 8, height = 5)
grid.arrange(db.data.umap, db.pc.umap, ncol = 2)
dev.off()

CairoPNG(filename = "Fig.genomics.Cldb2.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.Cldb2.pdf", width = 8, height = 5)
db.data.tsne
dev.off()

#Plot MDS
set.seed(123)
data.mds <-data.frame(X1 = cmdscale(data.dist)[,1], X2 = cmdscale(data.dist)[,2], mort = phenotype.v, outlier = outlier.v) 
pc.mds <- data.frame(X1 = cmdscale(pc.dist)[,1], X2 = cmdscale(pc.dist)[,2], mort = phenotype.v , outlier = outlier.v)

data.mds <- cbind(data.mds, cluster.db = data.dbscan)
pc.mds <- cbind(pc.mds, cluster.db = pc.dbscan)

db.data.mds <- ggplot(data.mds, aes(X1, X2, col = cluster.db, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'DBSCAN (eps= 5; minPts= 5) color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) + 
  scale_color_manual(values = colors[1:6])

db.data.mds


db.pc.mds <- ggplot(pc.mds, aes(X1, X2, col = cluster.db, shape = outlier )) + 
  geom_point() + 
  labs(title='PC',
       caption= 'DBSCAN (eps= 0.1; minPts= 5) color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) +
  scale_color_manual(values = colors[1:2])

db.pc.mds

CairoPNG(filename = "Fig.genomics.Cldb3.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.Cldb3.pdf", width = 8, height = 5)
grid.arrange(db.data.mds, db.pc.mds, ncol = 2)
dev.off()


#HDBSCAN
data.hdbscan <- as.factor(hdbscan(data, minPts = 5)$cluster)
pc.hdbscan <- as.factor(hdbscan(pc, minPts = 5)$cluster)

levels(data.hdbscan) <- c('Noise', 'Cluster 1', 'Cluster 2', 'Cluster 3', 
                          'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7',
                          'Cluter 8', 'Cluster 9', 'Cluster 10', 'Cluster 11')
levels(pc.hdbscan)<-  c('Noise', 'Cluster 1', 'Cluster 2', 'Cluster 3')

set.seed(123)
data.umap <-data.frame(umap(data)$layout, mort = phenotype.v, outlier = outlier.v) 
pc.umap <- data.frame(umap(pc)$layout, mort = phenotype.v , outlier = outlier.v)

data.umap <- cbind(data.umap, cluster.hdb = data.hdbscan)
pc.umap <- cbind(pc.umap, cluster.hdb = pc.hdbscan)

hdb.data.umap <- ggplot(data.umap, aes(X1, X2, col = cluster.hdb, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'HDBSCAN (minPts= 5) color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) 

hdb.data.umap


hdb.pc.umap <- ggplot(pc.umap, aes(X1, X2, col = cluster.hdb, shape = outlier )) + 
  geom_point() + 
  labs(title='PC',
       caption= 'HDBSCAN (minPts= 5) color') +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) 


hdb.pc.umap
CairoPNG(filename = "Fig.genomics.Clhdb1.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.Clhdb1.pdf", width = 8, height = 5)
grid.arrange(hdb.data.umap, hdb.pc.umap, ncol = 2)
dev.off()

#Plot MDS
set.seed(123)
data.mds <-data.frame(X1 = cmdscale(data.dist)[,1], X2 = cmdscale(data.dist)[,2], mort = phenotype.v, outlier = outlier.v) 
pc.mds <- data.frame(X1 = cmdscale(pc.dist)[,1], X2 = cmdscale(pc.dist)[,2], mort = phenotype.v , outlier = outlier.v)

data.mds <- cbind(data.mds, cluster.hdb = data.hdbscan)
pc.mds <- cbind(pc.mds, cluster.hdb = pc.hdbscan)

hdb.data.mds <- ggplot(data.mds, aes(X1, X2, col = cluster.hdb, shape = outlier )) + 
  geom_point() + 
  labs(title='PC and Covariates',
       caption= 'HDBSCAN (minPts= 5) color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) 

hdb.data.mds


hdb.pc.mds <- ggplot(pc.mds, aes(X1, X2, col = cluster.hdb, shape = outlier )) + 
  geom_point() + 
  labs(title='PC',
       caption= 'HDBSCAN (minPts= 5) color') +
  xlab('MDS1') + ylab('MDS2') +
  theme_light() +
  scale_shape_manual(values = c('Normal' = 19, 'Outlier' = 8)) 

hdb.pc.mds

CairoPNG(filename = "Fig.genomics.Clhdb2.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.Clhdb2.pdf", width = 8, height = 5)
grid.arrange(hdb.data.mds, hdb.pc.mds, ncol = 2)
dev.off()

CairoPNG(filename = "Fig.genomics.Clhdb3.png", width = 800, height = 500, pointsize = 12)
CairoPDF(file = "Fig.genomics.Clhdb3.pdf", width = 8, height = 5)
set.seed(123)
plot(hdbscan(data, minPts = 5), gradient = colors )
dev.off()





