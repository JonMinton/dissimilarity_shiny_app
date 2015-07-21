# Load prerequisites 

require(reshape2)
require(plyr)
require(stringr)
require(ggplot2)
require(maptools)
require(grid)
require(spdep)
require(Rcpp)
require(MASS)
require(CARBayes)
require(shiny)
require(dplyr)

# Load data 

la_to_dz <- read.csv("data/la_to_dz.csv") 

las <- c("all",
         la_to_dz %>%
          group_by(local_authority) %>%
          summarise %>%
          .$local_authority %>%
          as.vector
    )




shinyUI(fluidPage(
  titlePanel("D inference app"),
  sidebarLayout(
    
    sidebarPanel(
      #  This section contains the user inputs
      # Choose a dataset 
      h1("Data selection"),
      cat("ui:selectInput:option\n"), selectInput(
        "option",
        "Choose a dataset",
        choices=c(
          "country of origin, 2001"="coo_2001",
          "religion, 2001"="rel_2001",
          "ethnicity, 2001"="eg_2001",
          "country of origin, 2011"="coo_2011",
          "religion, 2011"="rel_2011",
          "ethnicity, 2011"="eg_2011" 
        ),
        selected=""
      ),
      
      # Select local authority (or all areas)
      cat("ui:selectInput:option_la\n"), 
      selectInput(
        "option_la",
        "Choose a Local Authority",
        choices=las,
        selected="all"
      ),
      
      # select numerator or denominator using uiOutput to make reactive to above selection
      cat("ui:uiOutput:numerator\n"), uiOutput("numerator"),
      cat("ui:uiOutput:denominator\n"), uiOutput("denominator"),
      # confirm selection button
      cat("ui:actionButton:ok_num_demon\n"), actionButton("ok_num_denom", "compile selection"),
      br(),
      
      # load shapefile and generate W matrix
      cat("ui:actionButton:load_shapefile_button\n"), actionButton("load_shapefile_button", "click to load the shapefile"),
      cat("ui:actionButton:make_w_matrix_button\n"), actionButton("make_w_matrix_button", "click to generate the neighbourhood matrix"),

      hr(),
      h1("Model tweaking"),
      br(),
      cat("ui:sliderInput:posterior_sample_size\n"), sliderInput("posterior_sample_size", "choose posterior sample size",
                  min=1000, max=10000, step=1000, value=1000),
      br(),
      cat("ui:actionButton:generate_posterior_button\n"), actionButton("generate_posterior_button", "click to run model"),
      

      cat("ui:sliderInput:seg_k\n"), sliderInput("seg_k", "Choose segregation thresholds",
                   min=0, max=1, value=c(0,1))
      ),
    
    mainPanel(
      h1("Checks"),
      p(
        "This secton shows some outputs that allow you to see whether the ",
        "prerequisites required by the model have been loaded successfully"
      ),
      cat("ui:htmlOutput:all_checks\n"), htmlOutput("all_checks"),

      p("The first few rows of the selected data are shown below"),
      cat("ui:tableOutput:show_combined_input_table\n"), tableOutput("show_combined_input_table"),
      br(),
      hr(),
      h1("Analysis"),
      p("This section will present results from the model run"),
      p("Once the model has been run a figure will be generated"),
      p("Along with other summary statistics"),
      
      cat("ui:textOutput:model_report\n"), htmlOutput("model_report"),
      h2("show plot"),
      cat("ui:plotOutput:show_posterior_distribution\n"), plotOutput("show_posterior_distribution"),
      cat("ui:tableOutput:tabulate_posterior\n"), tableOutput("tabulate_posterior")

      )
    )
    
))
  

