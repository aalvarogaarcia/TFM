#Rscript to prepare dataset1 (data of proteomics from patients of Interstitial Lung Diseases)
#Please, be aware that this data is confidential and cannot be released without previous consent
#from the owners

#Prepare board
rm(list = ls())
getwd()
path <- getwd()
path <- "/home/jlorsal/Downloads"
setwd(path)

#Load library
library(data.table)

#Define file datasets
fn.data <- "Proteomics.csv"
data <- fread(fn.data, sep=",", header= FALSE)
data[1:10,1:5]
data <- data[-3,]
data[1:10,1:5]
dim(data)
#[1] 2128 2940

#Panels of proteins
panels <- data[1,]
panels
head(panels)
str(panels)

#Proteins
proteins <- data[2,]
proteins
head(proteins)

colnames(data[1:10,1:5])
colnames(data)[1] <- "id"
colnames(data[1:10,1:5])

#Load senstitive data
fn.ids <- "list433"
ids <- fread(fn.ids)
head(ids)
colnames(ids)[1] <- "id"
head(ids)
dim(ids)

#Anomimize sensitive data
fn.ids.anon <- "list433.anonimized"
ids.anon <- fread(fn.ids.anon)
dim(ids.anon)
head(ids.anon)
colnames(ids.anon) <- "id.anon"
colnames(ids.anon)

#Merge data
library(dplyr)
data2 <- dplyr::inner_join(as.data.frame(ids), as.data.frame(data), by = "id")
data2[1:10,1:5]
dim(data2)
#[1]  432 2940
dim(data2)[1]
dim(data2)[2]

#Insert anonimized ids into the first column
data2$id <- ids.anon$id
data2[1:10,1:5]

dim(panels)
colnames(panels)[1] <- "id" 
colnames(panels)
data3 <- rbind(panels, data2)
data3[1:10,1:5]

#Subset data3 to grab only columns with proteins belonging to'Inflammatory' panel
#Transpose data
data4 <- t(data3)
data4[1:10,1:10]
data5 <- data4[data4[,1] == "Inflammation",]
dim(data5)

#Transposte again to return to the original order inds x proteins
data6 <- t(data5)
data6[1:10,1:10]

#Prepare output data
data7 <- cbind(id = ids.anon$id.anon, data6[-1,])
data7[1:10,1:10]
dim(data7)
#[1] 432 370

#Name columns as anonimized proteins
c <- 1:(dim(data7)[2]-1)
colnames(data7)
colnames(data7) <- c("id", c)

#Save data7 as your working dataset and enjoy clustering! Bazinga!
fn <- "Olink.NPX_values.Inflammatory-panels.csv"
write.csv(data7, fn, quote = FALSE, row.names = FALSE)

##End Of Script