##### PANELS FOR UI ######

#INPUT PANEL
main_panel <-tabsetPanel(
    id = 'main-panel',
    type = 'hidden',
    
    tabPanel('Unsupervised',
             textInput(
               "directory",
               "Directory of Data:",
               "/Users/alvarogarciamunoz/Documents/GitHub/TFM-2.0/Getting started"
               #"C:/Users/agarm/Desktop/TFM 2.0/TFM-2.0/Data"
             ),
             
             #Genotype text input
             textInput(
               "genotype",
               "Genotype-Data file name:",
               "Olink.NPX_values.Inflammatory-panels.csv"
             ),
             
             #Phenotype text input
             textInput(
               "phenotype",
               "Phenotype-Data file name:",
               "phenotypes-int.tsv"
             ),
             
             #Seed
             numericInput("seed", "Set seed:", 123),
             
             #Matrix selected
             selectInput(
               "matrix", 
               "Select a matrix:",
               c("Extended Matrix", "Genotype", 'Phenotype'),
               selected = "Extended Matrix"
             ), 
             
             #Dim Red method selection
             selectInput(
               "dimred",
               "Select dim reduction method:",
               c("UMAP", "t-SNE", "PCA"),
               selected = "UMAP"
             ),
             
             selectInput(
               "cluster",
               "Select clustering method (unsupervised):",
               c("k-means", "ConsensusClusters", "HDBScan", 'None'),
               selected = "k-means"
             )
    ),
    tabPanel('Supervised',
             textInput(
               "directory",
               "Directory of Data:",
               "/Users/alvarogarciamunoz/Documents/GitHub/TFM-2.0/Getting started"
               #"C:/Users/agarm/Desktop/TFM 2.0/TFM-2.0/Data"
             ),
             
             #Genotype text input
             textInput(
               "genotype",
               "Genotype-Data file name:",
               "Olink.NPX_values.Inflammatory-panels.csv"
             ),
             
             #Phenotype text input
             textInput(
               "phenotype",
               "Phenotype-Data file name:",
               "phenotypes-int.tsv"
             ),
             
             numericInput('seed', 'Set seed:', 123),
             
             numericInput('tSize', 'Test group size:', 1/3, max = 1, min = 0),
             
             selectInput(
               "dimred",
               "Select dim reduction method:",
               c("UMAP", "t-SNE", "PCA"),
               selected = "UMAP"
             )
        ),
    tabPanel('None',
             textInput(
               "directory",
               "Directory of Data:",
               "/Users/alvarogarciamunoz/Documents/GitHub/TFM-2.0/Getting started"
               #"C:/Users/agarm/Desktop/TFM 2.0/TFM-2.0/Data"
             ),
             
             #Genotype text input
             textInput(
               "genotype",
               "Genotype-Data file name:",
               "Olink.NPX_values.Inflammatory-panels.csv"
             ),
             
             #Phenotype text input
             textInput(
               "phenotype",
               "Phenotype-Data file name:",
               "phenotypes-int.tsv"
             ),
             
             selectInput(
               "matrix", 
               "Select a matrix:",
               c("Extended Matrix", "Genotype", 'Phenotype'),
               selected = "Extended Matrix"
             ), 
             
             selectInput(
               "dimred",
               "Select dim reduction method:",
               c("UMAP", "t-SNE", "PCA"),
               selected = "UMAP"
             )     
    )
  )
  
  #Text input for directory
  
 




#HIDDEN PANEL FOR PARAMETERS

parameters_panel <- tabsetPanel(
  id = "params",
  type = "hidden",
  
  #K-means Parameters
  tabPanel("k-means",
           numericInput("centers", "How many centers?:", value = 2, min = 2)
  ),
  
  #Consensus Parameters
  tabPanel("ConsensusClusters",
           numericInput("maxk", "Maximum centers?:", value = 20),
           numericInput("reps", "How many reps?:", value = 50),
           numericInput("pitem", "Which p value?:", value = 0.8, min = 0, max = 1),
           selectInput("clusteralg", "What algorithm?:", 
                       c("pam", "km"), "km"),
           selectInput("distance", "What distance?:",
                       c("pearson", "spearman", "euclidean", "binary", "maximum", "canberra", "minkowski"), "euclidean"),
           
  ),
  
  #Parameters for HDBScan
  tabPanel("HDBScan",
           numericInput("minpts", "Minimum points?:", value = 5, min = 1)
  ),
  
  tabPanel("None")
)
