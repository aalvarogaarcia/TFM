#SCOURGE individuals with COVID-19 treated with corticoids 
#Genetics defined by 10 PCs
#Covariates: sex, age
#All individuals: treatment with corticoids
#Phenotype: mortality-90 days

#Prepare board
rm(list = ls())
getwd()
path <- "/mnt/756b4178-6cf1-4cab-bdfb-8342cc4c885a/jlorsal/ANALISIS_275GB/Alvaro/03_dataset2"

setwd(path)

#source /etc/profile.d/zx-modules.sh
#module load Python/3.10.8-GCCcore-12.2.0
#module load R/4.2.2-foss-2022b

# Load necessary libraries
library(data.table)
library(ggplot2)
library(umap)
library(dbscan)

# Load and prepare the dataset
fn <- "corticoids_n1487.tsv"
data <- fread(fn, header = TRUE, stringsAsFactors = FALSE, quote = "")
head(data)

#Add a fake mortality variable bootstrapping the values in "mort"
set.seed(123)  # for reproducibility
data$fake.mort <- sample(data$mort)
colnames(data)
table(data$mort)
table(data$fake.mort)

# Remove rows with missing values
data2 <- na.omit(data)
dim(data2)

#Discard data without mortality data (coded as '-9')
data3 <- subset(data2, data$mort != "-9")
dim(data3)
id <- data3[,1]
id
dim(id)

# Convert all columns to numeric
data4 <- data3[,-1]
colnames(data4)
dim(data4)

data5 <- data.frame(lapply(data4, function(x) as.numeric(as.character(x))))
str(data5)
dim(data5)
colnames(data5)
head(data5)

#Stratify age and add to the matrix
data5$agegroup <- cut(data4$age, seq(0, 110, 10))
dim(data5)
colnames(data5)

#Plot PCs
plot(data5$pc1, data5$pc2, col = factor(data5$center))
which(data5$pc2 < -0.4)
#[1] 762 813
#These individuals show be removed!
dim(data5)
data5 <- data5[-c(which(data5$pc2 < -0.4)),]
dim(data5)

#########
### PCA for dimensionality reduction
colnames(data5)
col.to.exclude <- c(1,14,15,16)
pca <- prcomp(data5[,-col.to.exclude], center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca$x)
pca_data

# Visualize PCA
ggplot(data5, aes(x = pca$x[,1], y = pca$x[,2], color = as.factor(center))) +
  geom_point() +
  labs(color = "Center") +
  theme_minimal()
#which(pca$x[,2] < -5)

# This is the MATRIX OF DATA
# Scale the data (avoid the first columns of ids)
data6 <- as.matrix(scale(data5[,c(1:13)]), center = TRUE, scale = TRUE)
colnames(data6)
dim(data6)
data6[1:10,1:dim(data6)[2]]
str(data6)
hist(data6[,1], nclass=100)  #pc1
hist(data6[,10], nclass=100) #pc10
hist(data6[,11], nclass=100) #sex
hist(data6[,12], nclass=100) #age
hist(data6[,13], nclass=100) #center

#Check if missing data are present
which(is.na(data6))
#integer(0)

# Perform K-means clustering with K=3
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(data6, centers = 6)

# Print the cluster assignments
kmeans_result$cluster

# Plot the clusters using ggplot2
ggplot(data6, aes(x = pc1, y = pc2, color = factor(kmeans_result$cluster))) + 
  geom_point() +
  labs(color = "Cluster") +
  theme_minimal()

#Color by sex
ggplot(data6, aes(x = pc1, y = pc2, color = factor(sex))) + 
  geom_point() +
  labs(color = "Sex") +
  theme_minimal()

########
### UMAP
########
# Check if all columns are numeric
colnames(data6)
all(sapply(data6, is.numeric))
str(data6)
data6[1:10,1:dim(data6)[2]]
umap_result <- umap(data6, preserve.seed = TRUE)

# Visualize UMAP
ggplot(data5, aes(x = umap_result$layout[,1], y = umap_result$layout[,2])) +
  geom_point() +
  theme_minimal()

# Visualize UMAP colored by center
ggplot(data5, aes(x = umap_result$layout[,1], y = umap_result$layout[,2],
                  color = as.factor(center))) +
  geom_point() +
  labs(color = "Center") +
  theme_minimal()

# Visualize UMAP colored by mort
ggplot(data5, aes(x = umap_result$layout[,1], y = umap_result$layout[,2],
                  color = as.factor(mort))) +
  geom_point() +
  labs(color = "Mort") +
  theme_minimal()

# Visualize UMAP colored by fake mort
ggplot(data5, aes(x = umap_result$layout[,1], y = umap_result$layout[,2],
                  color = as.factor(data5$fake.mort))) +
  geom_point() +
  labs(color = "Fake ot") +
  theme_minimal()

# Visualize UMAP colored by age group
ggplot(data5, aes(x = umap_result$layout[,1], y = umap_result$layout[,2],
                  color = as.factor(data5$agegroup))) +
  geom_point() +
  labs(color = "Age group") +
  theme_minimal()

write.table(umap_result$layout, "umap.txt", quote = FALSE, sep = "\t")

#Clustering
#See: https://cran.r-project.org/web/packages/dbscan/vignettes/hdbscan.html
cl <- hdbscan(umap_result$layout, minPts = 10)
plot(umap_result$layout, col = cl$cluster, pch=20)

### Details of PCA
##Load eigenvalues (a vector)
eigenvalues <- scan(file = "corticoids_n1487.eigenval")
length(eigenvalues)

# Check eigenvalues are numeric
if (!is.numeric(eigenvalues)) {
  stop("Eigenvalues are not numeric.")
}

# Calculate the total variance
total_variance <- sum(eigenvalues)
total_variance

# Calculate explained variance for each PC
explained_variance <- eigenvalues / total_variance
explained_variance

# Calculate cumulative explained variance
cumulative_explained_variance <- cumsum(explained_variance)
cumulative_explained_variance

# Show explained variance for each PC
for (i in 1:length(explained_variance)) {
  cat(sprintf("Explained variance in PC%d: %.4f\n", i, explained_variance[i]))
}

# Show cumulative explained variance
for (i in 1:length(cumulative_explained_variance)) {
  cat(sprintf("Cumulative explained variance in PC%d: %.4f\n", i, cumulative_explained_variance[i]))
}

# Scatter plot showing variances of each plotted PC
label.x <- paste("PC1 (", round(explained_variance[1] * 100, 2), "%)", sep = "")
label.y <- paste("PC2 (", round(explained_variance[2] * 100, 2), "%)", sep = "")
ggplot(data = data6, aes(x = pc1, y = pc2)) +
  geom_point() +
  xlab(label.x) +
  ylab(label.y)

##End Of Script
