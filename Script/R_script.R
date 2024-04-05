#### R - SCRIPT #####

#Este documento recoge ordenadamente y con comentarios el código de R utilizado para el informe dinámico


#### INTRODUCCIÓN ####

#Librerías
rm(ls())
library(ggplot2) #Visualización de los datasets
library(tidyverse) #Manejo de los datasets
library(SNPRelate) #Librería básica encargada del análisis principal de los datos genéticos
library(reticulate) #Uso de python
library(RefManageR) #Bibliografía
library(openxlsx) #Gestión de archivos tipo excel
library(stats) #Incluye diferentes estadísticas
library(data.table)

#Prepare the session
rm(list = ls())


#Define paths
getwd()
path <- "write-your-path-here"
setwd(path)

fn <- "vcf1.vcf"
filename <- paste(path, fn, sep="/")
data <- fread(filename)

#Load data
fn_vcf <- "-.vcf"  #Aquí incluiremos la ruta al archivo VCF del cual estraeremos los datos
fn_cv <- "-.xlsx" #Aquí incluiremos las covariables que deseemos estudiar
fn_red <- "-" #This must be a plain text file

snpgdsVCF2GDS(fn_vcf, "chr22") #This command creates a GDS file
co.var <- read.xlsx(fn_cv) #This command reads an excel file of covariables, be careful in front of whom you open it.
red <- as.vector(read.csv(fn_red, header= FALSE)) #Leemos el archivo de reducción
red <- unlist(unname(red)) 


co.var.red<- filter(co.var, co.var$sample %in% red)

#Eliminación LD

snpset <- snpgdsLDpruning(GDS, ld.threshold = 0.2, method = "r")

snpset.id <- unlist(unname(snpset)) 

#PCA

pca <- snpgdsPCA(GDS, snp.id= snpset.id, num.thread = 10)

df <- data.frame(sample.id = pca$sample.id,
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

#Excel

pca_cv <- data.frame(df, covariables[,-1])
archivo_xlsx <- "PCA10_Covariables.xlsx"
write.xlsx(pca_cv,archivo_xlsx)

#Plot

#Paleta de colores de 26 clases
colores <- c("#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF",
             "#00FFFF", "#FF8000", "#80FF00", "#8000FF", "#FFFF80", "#FF80FF",
             "#80FFFF", "#008000", "#0080FF", "#800080", "#808080", "#C0C0C0",
             "#A0A0A0", "#606060", "#404040", "#202020", "#AA12AA", "#158752", 
             "#CCCC56", "#949513") 

ggplot(pca_cv, aes(EV2, EV1, color = pop, shape = super_pop)) + geom_point() + scale_color_manual(values=colores) + theme_light()








#### CLUSTERS #####

#K-means

k_valores <- c(5,26)
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



#Plot kmeans (es la misma de Plots)

df_kmeans$pop <- as.factor(kmeans_PCA$'26'$cluster)
df_kmeans$super_pop <- as.factor(kmeans_PCA$'5'$cluster)

par(mfrow=c(1,2))

ggplot(pca_cv, aes(EV2, EV1, color = pop, shape = super_pop)) + geom_point() + scale_color_manual(values=colores) + theme_light()

ggplot(df_kmeans, aes(EV2, EV1, color = pop, shape = super_pop)) + geom_point() + scale_color_manual(values=colores) + theme_light()








##End Of Script






