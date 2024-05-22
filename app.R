#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


#Libraries
library(shiny)
library(shinythemes)
library(tidyverse)
library(ggplot2)
library(umap)
library(Rtsne)
library(dbscan)
library(ConsensusClusterPlus)

source('pan.R')
source('fun.R')

  
  

ui <- fluidPage( theme = shinytheme('simplex'),
    # Application title
    titlePanel("Clustering UI."),
    inputPanel(
      #Parameters for clustering
      selectInput(
        'algo',
        'Select type of clustering: ', 
        c('Unsupervised', 'None'),
        selected = 'Unsupervised'),
      
      submitButton('Submit')
    ),

        sidebarPanel(
          main_panel,
          parameters_panel,
          submitButton()
        ),
        

        # Show a plot of the generated distribution
        mainPanel(
          h1("Plotting!!"),
          
          h4("First plot"),
          
          plotOutput("umap")
        )
    )


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  observeEvent(input$algo, {
    updateTabsetPanel(inputId = 'main-panel', selected = input$algo)
  })
  
  observeEvent(input$cluster, {
    updateTabsetPanel(inputId = "params", selected = input$cluster)
  })
  
  observeEvent(input$matrix, {
    updateTabsetPanel(inputId = "clusttype", selected = input$matrix)
  }) 
  
  #Read data
  gen.df <- reactive({
    if(!is.na(input$genotype)){
      fn.gen <- paste(input$directory, input$genotype, sep = "/")
      gen.df <- read.csv(fn.gen)
      
      #No more covariance 
      temp <- cov(gen.df[,-1])
      temp <- eigen(temp)
      temp <- temp$vector
      m <-c() 
      for(i in 1:nrow(gen.df[,-1])){
        pr <- temp * gen.df[i,-1]
        m <- rbind(m, pr)
      }
      gen.df <- cbind("id" = gen.df[,1], m)
    }
    gen.df
  })
  
 phen.df <- reactive({
    if(!is.na(input$phenotype)){
      fn.phen <- paste(input$directory, input$phenotype, sep = "/")
      phen.df <- read_tsv(fn.phen)
    }
  })
  
  output$umap <- renderPlot({
    #Matrix management
      
      var <- get_matrix(input$matrix, gen.df(), phen.df())
      
      var <- scale(var[,-1])
    #Set seed
      set.seed(input$seed)
      
    
    #Dimension reduction method
      coord <- get_layout(input$dimred, var)
      
    #Get clusters
      if(input$algo == 'Unsupervised'){
        if(input$cluster == "k-means"){
          model <- kmeans(var[,-1], centers = input$centers)
          col <- as.factor(model$cluster)
        }
        else if(input$cluster == "ConsensusClusters"){
          model <- ConsensusClusterPlus(as.matrix(var[,-1]), 
                                        maxK = input$maxK, 
                                        reps = input$reps,
                                        pFeature = 1,
                                        title = paste(input$matrix, "ConsensusCluster"),
                                        pItem = input$pitem,
                                        clusterAlg = input$clusteralg,
                                        distance = input$distance,
                                        plot = "png")
          temp <- calcICL(model, title = paste(input$matrix, "ICL cluster"))
          temp <- temp[["itemConsensus"]][which(icl.data[["itemConsensus"]][,ncol(icl.data[["itemConsensus"]])]>input$pitem),]
          col <- as.factor(temp$k)
        }
        else if(input$cluster == "HDBScan"){
          model <- hdbscan(var[,-1], minPts = input$minpts)
          col <- as.factor(model$cluster)
        }
      }
      else{
        col <- NA
      }
      
      
      clusters <- col
      
      
      #Plotting
      
      if(length(clusters)<=1){
        plot(coord)
      }
      else{
        coord <- cbind(coord, clusters)
        plot(coord[,-ncol(coord)], col = coord[, ncol(coord)])
      }
      
    })
    }
    
  
    
    
                                       


# Run the application 
shinyApp(ui = ui, server = server)
