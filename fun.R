#### FUNCTIONS FOR UI #####
#SELECTING MATRIX
get_matrix <- function(x, y, z){
  if(x == "Extended Matrix"){
    var <- merge(y, z, by = colnames(y)[1])
  }
  else if(x == "Genotype"){
    var <- y
  }
  else if(x == 'Phenotype'){
    var <- z
  }
  return(var)
}

#SELECT TRAIN SET
get_train <- function(x, n){
 
}
#GET LAYOUTS DEPENDING ON METHOD
get_layout <- function(x, y){
  if(x == "UMAP"){
    df <- umap(y[,-1])
    coord <- df$layout
  }
  else if(x == "t-SNE"){
    df <- Rtsne(y[,-1], perplexity = 10, check_duplicates = FALSE)
    coord <- df$Y
  }
  else{
    df <- prcomp(y[,-1])
    coord <- cbind(df$x[,1], df$x[,2])
  }
}

