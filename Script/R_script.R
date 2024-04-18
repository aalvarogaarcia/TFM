#### R - SCRIPT #####

#Este documento recoge ordenadamente y con comentarios el código de R utilizado para el informe dinámico


#### INTRODUCCIÓN ####

#Librerías
library(ggplot2) #Visualización de los datasets
library(tidyverse) #Manejo de los datasets
library(SNPRelate) #Librería básica encargada del análisis principal de los datos genéticos
library(reticulate) #Uso de python
library(RefManageR) #Bibliografía
library(openxlsx) #Gestión de archivos tipo excel
library(stats) #Incluye diferentes estadísticas
library(e1071)#Supervised Learning
library(glmnet) #Supervised Learning
library(umap) #Dimension Reduction
library(Rtsne) #Dimension Reduction
library(tensorflow) #Unsupervised Learning
#library(data.table)

#Prepare the session
rm(list = ls())


#Define paths
getwd()
path <- "write-your-path-here"
setwd(path)

#### DATA MANAGEMENT ####

#Load data
fn_vcf <- "-.vcf"  #Aquí incluiremos la ruta al archivo VCF del cual estraeremos los datos
fn_cv <- "-.xlsx" #Aquí incluiremos las covariables que deseemos estudiar
fn_red <- "-" #This must be a plain text file

fn_GDS <- "chr22" #This will be how we name the GDS file

snpgdsVCF2GDS(fn_vcf, fn_GDS) #This command creates a GDS file
co.var <- read.xlsx(fn_cv) #This command reads an excel file of covariables, be careful in front of whom you open it.
red <- read.csv(fn_red, header= FALSE) #Leemos el archivo de reducción
red <- unlist(unname(red)) 

filter <- c("CEU", "FIN", "GBR", "IBS", "TSI", "CHB", "YRI")

co.var.red<- filter(co.var, co.var$sample %in% red)
co.var.red <- co.var[co.var$sample %in% red & co.var$pop %in% filter, ]

GDS <- snpgdsOpen(fn_GDS)
genome.matrix<- snpgdsGetGeno(GDS)
#Eliminación LD

snpset <- snpgdsLDpruning(GDS, ld.threshold = 0.2, method = "r")

snpset.id <- unlist(unname(snpset)) 

genome.matrix<- snpgdsGetGeno(GDS, snp.id = snpset.id)
#### DIMENSION REDUCTION ####

#PCA

pca <- snpgdsPCA(GDS, snp.id= snpset.id, num.thread = 10)

df <- data.frame(sample = pca$sample.id,
                 EV1=pca$eigenvect[,1],
                 EV2=pca$eigenvect[,2],
                 EV3=pca$eigenvect[,3],
                 EV4=pca$eigenvect[,4],
                 EV5=pca$eigenvect[,5],
                 EV6=pca$eigenvect[,6],
                 EV7=pca$eigenvect[,7],
                 EV8=pca$eigenvect[,8],
                 EV9=pca$eigenvect[,9],
                 EV10=pca$eigenvect[,10],
                 stringsAsFactors=FALSE)
perc.explain <- pca$eigenval/sum(pca$eigenval, na.rm = TRUE)

paste("PCA 1-70 ", round(sum(perc.explain*100, na.rm = TRUE), 2), "%")
paste("PCA1-10", round(sum(perc.explain[1:10]*100, na.rm = TRUE),2), "%")

#CSV

pca_cv <- merge(df, co.var.red, by = "sample")
fn_csv<- "PCA10_Covariables"
write.csv2(pca_cv,fn_csv)

#Ajuste de covariables


hot_encoding_pop <- model.matrix(~pop, data = pca_cv, contrasts = "contr.treatment")
hot_encoding_gender <- model.matrix(~gender, data = pca_cv, contrasts = "contr.treatment")
hot_encoding_super_pop <- model.matrix(~super_pop, data = pca_cv, contrasts = "contr.treatment")

df_hot_encoding <- cbind(df, hot_encoding_pop[,-1], genderMale = hot_encoding_gender[,-1], super_pop =pca_cv$super_pop)

library(dummy)
co.red.dummies <- dummy(co.var.red[,-1], int = TRUE)
co.red.dummies <- cbind(sample = co.var.red[,1], co.red.dummies)
#Plot

#Paleta de colores de 26 clases
colores <- c("#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF",
             "#00FFFF", "#FF8000", "#80FF00", "#8000FF", "#FFFF80", "#FF80FF",
             "#80FFFF", "#008000", "#0080FF", "#800080", "#808080", "#C0C0C0",
             "#A0A0A0", "#606060", "#404040", "#202020", "#AA12AA", "#158752", 
             "#CCCC56", "#949513") 

ggplot(pca_cv, aes(EV2, EV1, color = pop, shape = super_pop)) + geom_point() + scale_color_manual(values=colores) + theme_light()


#UMAP

# Data and label
df.data <- df_hot_encoding[, grep("EV|pop|gender", colnames(df_hot_encoding))]
df.labels <- df.data[, ncol(df.data)]

df.data <- df.data[, -ncol(df.data)]

#Ploting

df.umap <- umap(df.data)
df.umap

coordinates.umap <- data.frame(df.umap$layout, labels = df.labels)

ggplot(coordinates.umap, aes(X1, X2, color = labels)) + geom_point() + theme_light()

#UMAP learn - Python
df.umap.learn <- umap(df.data, method = "umap-learn") #Needs a python module


#T-SNE

df.tsne <- Rtsne(df.data, perplexity = 10, check_duplicates = FALSE)

coordinates.tsne <- data.frame(df.tsne$Y, labels = df.labels)

ggplot(coordinates.tsne, aes(X1, X2, color = labels)) + geom_point() + theme_light()


#### CLUSTERS #####

#### Unsupervised - Learning ####

#K-means

k_valores <- c(3,9)
kmeans_PCA <- list()

set.seed(123)
for (k in k_valores){ #El bucle nos da libertad de hacer varios agrupamientos dependiendo del número de grupos que queramos
  agrupamiento <- kmeans(df[,-1], centers = k)
  k<-as.character(k)
  kmeans_PCA[[k]] <- agrupamiento  
}

for (i in kmeans_PCA){
  print(i$size) #Muestro en pantalla el tamaño de los grupos
}

table(as.factor(pca_cv$super_pop)) #Ha de ajustarse dependiendo del estudio pero comprobamos si coinciden los grupos
table(as.factor(pca_cv$pop))

print(agrupamiento$size)
df_hot_encoding$cluster <- as.factor(agrupamiento$cluster)

#Plot kmeans (es la misma de Plots)
df_kmeans <- df

df_kmeans$pop <- as.factor(kmeans_PCA$'9'$cluster)
df_kmeans$super_pop <- as.factor(kmeans_PCA$'3'$cluster)

par(mfrow=c(1,2))

ggplot(pca_cv, aes(EV2, EV1, color = pop, shape = super_pop)) + geom_point() + scale_color_manual(values=colores) + theme_light()

ggplot(df_kmeans, aes(EV2, EV1, color = pop, shape = super_pop)) + geom_point() + scale_color_manual(values=colores) + theme_light()

ggplot(df_hot_encoding, aes(EV2, EV1, color = cluster)) + geom_point()


#ConsensusClusterPlus
 library(ConsensusClusterPlus)

  df.data <- merge(df, co.red.dummies, by = "sample")
  df.data <- cbind(genome.matrix, df.data[,-1])
  
  
  #Clustering
  title = tempdir()
  results = ConsensusClusterPlus(df.data, maxK = 26, reps = 50, pItem = 0.8, pFeature = 1, 
                                 title = title, clusterAlg = "hc", distance = "pearson", plot = "png")

  
  #ICL
  icl <- calcICL(results, title = title, plot='png')
  icl[['clustersConsensus']]

  #### Supervised - Learning ####


#SVM

#Split data
index<- 1:nrow(df_hot_encoding)
N <- trunc(length(index)/3)
testindex <- sample(index, N)
testset <- df.data[testindex,]
trainset <- df.data[-testindex,]
trainset$super_pop <- as.factor(df_hot_encoding$super_pop[-testindex])

#Model
svm.model <- svm(super_pop~ ., data = trainset[,-1], cost = 100, gamma = 1)
svm.pred <- predict(svm.model, testset)

conf.matrix <-table(pred = svm.pred, true = df_hot_encoding$super_pop[testindex]) #compute SVM confusion matrix
conf.matrix
classAgreement(conf.matrix)

testset$cluster <- svm.pred
testset$super_pop <- df_hot_encoding$super_pop[testindex]

ggplot(testset, aes(EV2, EV1, color = cluster)) + geom_point()
ggplot(testset, aes(EV2, EV1, color = super_pop)) + geom_point()


#GLMNET

#glmnet object
v.labels <- as.numeric(as.factor(df.labels)) #We transform our labels to a numeric vector
fit <- glmnet(df.data, v.labels) #As default we use a Gaussian linear model (least squares)
#As default the fit would be a lasso regression (\alpha = 1).


#fit shows
plot(fit)
print(fit)
coef(fit)

#cv.glmnet object
df.matrix <- as.matrix(df.data) #We transform our df to a matrix
cvfit <- cv.glmnet(df.matrix, v.labels) #Cross validation

#cv.fit shows
plot(cvfit)
cvfit$lambda.min
coef(cvfit)
predict(cvfit, newx = df.matrix[1:30,])

# Add penalty factors
p.fac <- rep(1,ncol(df.matrix))
p.fac[c(11, 17, 10)]<- 0
pfit <- glmnet(df.matrix, v.labels, penalty.factor = p.fac)
plot(pfit, label=TRUE)

#Confusion matrices
cfit <- cv.glmnet(df.matrix[-testindex,], v.labels[-testindex], family="multinomial")
cnf <- confusion.glmnet(cfit, newx = df.matrix[testindex,], newy = v.labels[testindex])
print(cnf)

#ROC *glmnet* (Use only with binomial family)
bfit <- cv.glmnet(df.matrix, v.labels, family="binomial", type.measure = "auc", keep=TRUE)
roc <- roc.glmnet(bfit$fit.preval, newy = v.labels)

#### Deep - Learning ####

# Preparing data
x.matrix.train <- as.matrix(df.data[-testindex, ])
x.matrix.test <- as.matrix(df.data[testindex, ])

y.train <- hot_encoding_super_pop[-testindex, -1]
y.test <- hot_encoding_super_pop[testindex, -1]

write.csv(as.matrix(df.data), "chr22-matrix")
write.csv(v.labels, "chr22-labels")


  # Defining model (sequential)
    
    model <- keras_model_sequential() 
    model %>%  #Defining the model for our study
      layer_dense(units = 256, activation = 'relu', input_shape = c(23,17)) %>% 
      layer_dropout(rate = 0.4) %>% 
      layer_dense(units = 128, activation = 'relu') %>%
      layer_dropout(rate = 0.3) %>%
      layer_dense(units = 10, activation = 'softmax')

    model %>% compile(
      loss = 'categorical_crossentropy',
      optimizer = optimizer_rmsprop(),
      metrics = c('accuracy')
    )

##End Of Script
snpgdsClose(GDS)





