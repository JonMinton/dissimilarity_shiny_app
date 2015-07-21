
# Notes and to-dos


# 1) include checks for if numerators are greater than denominators
# 2) include feedback to tell user to wait when performing large calculations
# 3) change order of US outputs to make more intuitive sense
# 4) 

# server - load prerequisities --------------------------------------------


require(spdep)
require(maptools)
require(CARBayes)
require(Rcpp)
require(MASS)

require(stringr)
require(plyr)
require(tidyr)
require(dplyr)


require(grid)
require(ggplot2)

require(shiny)





# Server - load scripts  --------------------------------------------------
source("scripts/D_compute.r")
sourceCpp("scripts/cppfunctions.cpp")
# 


# server - load data  -----------------------------------------------------


la_to_dz <- read.csv("data/la_to_dz.csv") 

las <- c("all",
         la_to_dz %>%
           group_by(local_authority) %>%
           summarise %>%
           .$local_authority %>%
           as.vector
)

#####################################################################################
#####################################################################################

server <- function(input, output, server){
  
  ###############################################################################################
  #####  Observer FUNCTIONS #####################################################################
  ###############################################################################################
  run_model <- reactive(
    {    
      cat("server:reactive:run_model\n")
      out <- NULL
      w <- generate_w_matrix()
      dta <- combine_input_table()
      n_burnin <- input$select_burnin
      n_sample <- input$select_sample
      
      if (!is.null(w) & !is.null(dta)){
        cat("in inner of run_model\n")
        browser()
        # remove duplicates 
        w_dz <- rownames(w)
        w <- w[!duplicated(w_dz), !duplicated(w_dz)]
        
        # remove areas with no neighbours
        w_rowsum <- apply(w, 1, sum)
        w <- w[w_rowsum > 0, w_rowsum > 0]
        

        w_dz <- rownames(w)
        dta_dz <- dta$datazone %>% as.character

        ss <- intersect(w_dz, dta_dz)
        tmp <- w_dz %in% ss
        w <- w[tmp, tmp]
        dta <- subset(dta, datazone %in% ss)
        
        mdl <- S.CARiar(
          formula= numerator  ~ 1,
          trials = dta$denominator,
          W=w,
          data=dta,
          family="binomial",
          burnin = n_burnin,
          n.sample = n_sample, 
          verbose = TRUE
        )
        
        
        out <- list(
          datazones = dta_dz,
          denominator = dta$denominator,
          beta = mdl$samples$beta[,1],
          phi = mdl$samples$phi
        )
      }
      return(out)
    })
  
  generate_posterior_distribution <- eventReactive(
    input$generate_posterior_button,
    {
      cat("server:eventReactive:generate_posterior_distribution\n")
      model_outputs <- run_model()
      
      K <- input$posterior_sample_size      
      phi <- model_outputs$phi
      beta <- model_outputs$beta
      denominator <- model_outputs$denominator
      
      out <- array(NA, K)
      for(k in 1:K){
        p.current <- exp(phi[k ,] + beta[k])   / (1 + exp(phi[k ,] + beta[k]))
        p.current.overall <- sum(p.current * denominator) / sum(denominator)
        out[k] <- sum(denominator * abs(p.current - p.current.overall)) / 
          ( 2 * sum(denominator) * p.current.overall * (1-p.current.overall))                 
      }
      return(out)
    })
  
  
  summarise_posterior_distributions <- eventReactive (
    input$generate_posterior_button,
    {
      cat("server:eventReactive:summarise_posterior_distribution\n")
      bayes <- generate_posterior_distribution()
      classical <- calc_d_classical()$D.boot
      n_digits <- 4
      
      if (!is.null(bayes) & !is.null(classical)){
        qi_bayes <- bayes %>% quantile(c(0.025, 0.5, 0.975)) %>% round(4)
        qi_classical <- classical %>% quantile(c(0.025, 0.5, 0.975)) %>% round(4)
        
        out <- data.frame(
          method=c("Classical", "Bayesian"),
          lower=c(qi_classical[1], qi_bayes[1]),
          middle=c(qi_classical[2], qi_bayes[2]),
          upper=c(qi_classical[3], qi_bayes[3])
        )
        
      } else {out <- NULL}
      return(out)        
    })
  
  
  load_shapefiles <- eventReactive(
    input$load_shapefile_button, 
    {
      cat("load_shapefiles\n")
      out <- readShapeSpatial(
        "shp/scotland_2001_datazones/scotland_dz_2001.shp"
      )      
      out@data <- rename(out@data, datazone=zonecode)
      
      return(out)
    })
  
  link_shp_with_attributes <- reactive({
    cat("link_shp_with_attributes\n")
    
    shp_data <- load_shapefiles()    
    att_data <- combine_input_table()
    
    if (!is.null(shp_data) & !is.null(att_data)){
      out <- shp_data
      # The data can be loaded
      out@data <- plyr::join(out@data, att_data)
    } else {
      # The data cannot be linked
      out <- NULL
    }
    return(out)
  })
  
  load_data <- reactive({
    cat("load_data\n")    
    file_name <- input$option
    data <- read.csv(paste0("data/", file_name, ".csv"))
    if (input$option_la!="all"){
      
      dzs <-  la_to_dz  %>% 
        filter(local_authority == input$option_la)  %>%
        .$datazone %>%
        as.vector
      
      data <- data %>%
        filter(datazone %in% dzs)
    }
    return(data)
  })
  
  summarise_data <- reactive({
    cat("summarise\n")    
    data <- load_data()
    out <- data %>% group_by(type) %>% summarise(count=sum(count)) %>% data.frame
    return(out)
  })
  
  get_labels <- reactive({
    cat("get_labels\n")
    labels <- levels(load_data()$type) %>% as.character() %>% sort()
    return(labels)
  })
  
  
  
  generate_w_matrix <- eventReactive(
    input$make_w_matrix_button, 
    {
      cat("generate_w_matrix\n")
      out <- NULL
      dta <- link_shp_with_attributes()
      if (!is.null(dta)){
        w_nb <- poly2nb(dta)
        out <- nb2mat(w_nb, style="B", zero.policy=TRUE)
        rownames(out) <- colnames(out) <- dta@data$datazone        
      }  
      return(out)
    })
  
  combine_input_table <- eventReactive(
    input$ok_num_denom,
    {
      cat("combine_input_table\n")
      
      numerators <- input$numerator_selection
      denominators <- input$denominator_selection
      out <- NULL
      if (!is.null(numerators) & !is.null(denominators)){
        cat("numerators: ")
        for (i in 1:length(numerators)) {cat(numerators[i], "\n")}
        cat("\n\n")  
        cat("denominators: ")
        for (i in 1:length(denominators)) {cat(denominators[i], "\n")}
        cat("\n\n")      
        
        data_raw <- load_data()
        data_raw <- data_raw %>% dplyr::select(datazone, type, count)
        
        data_numerator <- data_raw %>% 
          filter(type %in% numerators) %>%
          group_by(datazone) %>% summarise(numerator=sum(count))
        
        data_denominator <- data_raw %>% 
          filter(type %in% denominators) %>%
          group_by(datazone) %>% summarise(denominator=sum(count))
        
        data_out <- inner_join(data_denominator, data_numerator)
        out <- data_out                  
      } 
      return(out)
    })
  
  
  calc_d_classical <- reactive({
    cat("calc_d_classical\n")    
    dta <- combine_input_table()
    out <- Dissimilarity.compute(
      minority=dta$numerator,
      total=dta$denominator
    ) 
    return(out)
  })
  
  #############################################################################################
  ### REACTIVE UIS ############################################################################
  #############################################################################################
  output$numerator <- renderUI({
    cat("output:numerator\n")
    selections <- get_labels()
    selectInput("numerator_selection", "Select numerator", choices=selections, multiple=T)
  })
  
  output$denominator <- renderUI({
    cat("output:denominator\n")
    selections <- get_labels()
    selectInput("denominator_selection", "select denominator", choices=selections, multiple=T)
  })
  
  
  ##############################################################################################
  ### OUTPUTS ##################################################################################
  ##############################################################################################
  
  
  
  output$show_combined_input_table <- renderTable({
    cat("output:show_combined_input_table\n")      
    out <- combine_input_table() %>%
      as.data.frame %>%
      head
    return(out)
  })
  
  output$report_attributes_linked <- renderText({
    cat("output:report_attributes_linked\n")      
    dta <- link_shp_with_attributes() 
    if (is.null(dta)){
      out <- "The data has not been merged yet"
    } else {
      out <- "the data have been merged"
    }
    return(out)
  })
  
  output$all_checks <- renderUI({
    cat("output:all_checks\n")
    
    # This will replace a number of other renderText functions, creating a 
    # single dynamic paragraph that will report on the state of various 
    # prerequisites.
    
    # The structure of the output will read
    
    # (1)"The numerator and denominator [have/have not] been selected"
    # (2)"[The numerator and denominator selection is valid as no numerators are greater
    # than their denominators] / 
    # [The numerator and denominator selection is not valid because 
    # [[XXX]] numerators are greater than their denominators. Please choose again.]"
    # (3) "The shapefile [has/has not] been loaded"
    # (4) "The numerators/denominators [have/have not] been merged to the shapefiles."
    # (5)"[The shapefile is needed to calculate the neighbourhood matrix]/
    # [The shapefile has been loaded but the neighbourhood matrix has not yet been calculated]/
    # [The neighbourhood matrix has been calculated]
    # (6)"[All areal units are included. Warning: This may take some time"] /
    # [Only one local authority has been selected]
    # (7) [The model [has not]/[has] been run"
    
    # (1)"The numerator and denominator [have/have not] been selected"
    tmp <- combine_input_table()
    
    out_01 <- if (is.null(tmp)){
      "The numerator and denominator have not been selected"
    } else {
      paste(
        "The numerator and denominator have been selected and",
        "has", dim(tmp)[1], "rows"
      )
    }
    # (2)"[The numerator and denominator selection is valid as no numerators are greater
    # than their denominators] / 
    tmp2 <- tmp$denominator - tmp$numerator
    
    out_02 <- if (any(tmp2 < 0)){
      paste(
        "The numerator/denominator selection is invalid as",
        length(tmp2[tmp2 < 0]), "numerators are greater than",
        "the corresponding denominators. <strong>Please choose again</strong>"
      )
    } else {
      paste(
        "The numerator/denominator selection is valid as",
        "no denominators are smaller than the corresponding",
        "numerators"
      )
    }
    
    
    # (3) "[The shapefile has not been loaded]/
    # [The shapefile has been loaded and has length [[XXX]]"
    tmp3 <- load_shapefiles()
    
    out_03 <- if (is.null(tmp3)){
      "Shapefiles not yet loaded"
    } else {
      paste(
        "The shapefiles have been loaded and have length", length(tmp3)
      )
    }
    
    # (4) "The numerators/denominators [have/have not] been merged to the shapefiles."
    tmp4 <- link_shp_with_attributes()
    
    out_04 <- if(is.null(tmp4)){
      "The numerators/denominators have not been merged to the shapefiles"
    } else {
      "The numerators/denominators have been merged to the shapefiles"
    }
    
    # (5)"[The shapefile is needed to calculate the neighbourhood matrix]/
    # [The shapefile has been loaded but the neighbourhood matrix has not yet been calculated]/
    # [The neighbourhood matrix has been calculated]
    
    ###
    
    tmp5 <- generate_w_matrix()
    
    out_05 <- if(is.null(tmp4)){
      "The neighbour matrix cannot be created. Check the shapefile has been loaded"
    } else {
      "The neighbour matrix has been calculated"
    }
    
    # (6)"[All areal units are included. Warning: This may take some time"] /
    # [Only one local authority has been selected]
    
    tmp6 <- input$option_la
    
    out_06 <- if(tmp6==""){
      "No areal selection has been made yet"
    } else if (tmp6=="all"){
      "All of Scotland selected. This may take some time"
    } else {
      "Only one local authority selection. This will run fast but may not be representative"
    }
    
    
    output <- HTML(paste(
      "<p><b>Data report:</b></p>",
      "<ul>",
      "<li>", out_01, "</li>",
      "<li>", out_02, "</li>",
      "<li>", out_03, "</li>",
      "<li>", out_04, "</li>",
      "<li>", out_05, "</li>",
      "<li>", out_06, "</li>",
      "</ul>"
    ))
    return(output)
  })
  
  
  output$report_w_matrix_generated <- renderText({
    cat("output:report_w_matrix_generated\n")
    tmp <- generate_w_matrix()
    if (is.null(tmp)){
      out <- "The w matrix has not been generated"
    } else {
      out <- paste("The w matrix has been generated and has dimensions", 
                   dim(tmp)[1], " by ", dim(tmp)[2])
    }
    return(out)
  })
  
  output$show_posterior_distribution <- renderPlot({
    cat("output:show_posterior_distribution\n")      
    samples <- generate_posterior_distribution() 
    
    if (!is.null(samples)){
      samples <- data.frame(value=samples)
      thresholds <- input$seg_k
      samples$filled <- "yes"
      samples$filled[samples$value <= thresholds[1] | samples$value >= thresholds[2]] <- "no"
      prop_in_band <- length(samples$filled[samples$filled=="yes"])/length(samples$filled) %>%
        round(., 2)
      
      out <- samples %>% 
        ggplot(data=., aes(x=value, fill=filled)) + 
        geom_histogram(bindwidth = 0.001) +
        scale_fill_manual(values=c("yes"="black", "no"= "lightgray"), guide=FALSE) +
        geom_vline(
          xintercept=thresholds,
          linetype="dotted",
          linewidth=2.5
        ) + 
        coord_cartesian(xlim=c(0,1)) +
        annotate("text", x=mean(thresholds), y=-1.5, col="red", fontface="bold", label=prop_in_band)
      
    } else {out <- NULL}
    return(out)
  })
  
  output$tabulate_posterior <- renderTable({
    cat("output:tabulate_posterior\n")
    
    out <- summarise_posterior_distributions()
    
    return(out)
    
  })
  
  output$model_report <- renderUI({
    cat("output:model_report\n")
    
    
    # The structure of the output will read
    
    # (1)"Lower and upper threshold"
    # (2) Probability true value within threshold
    # (3) mean using classical and Bayesian
    # (4) CrIs using classical and Bayesian
    thresholds <- input$seg_k
    
    out_01 <- paste("The selected thresholds are", thresholds[1], " to", thresholds[2])
    
    samples <- generate_posterior_distribution()
    
    p_in_t <- length(samples[between(samples, thresholds[1], thresholds[2])]) / length(samples)
    p_in_t <- round(p_in_t, 2)
    
    
    out_02 <- paste("The probability that the true value is within the threshold selected is", p_in_t)
    
    
    output <- HTML(paste(
      "<p><b>Data report:</b></p>",
      "<ul>",
      "<li>", out_01, "</li>",
      "<li>", out_02, "</li>",
      "</ul>"
    ))
    return(output)
  })
  
  
}


ui <- fluidPage(
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
      cat("ui:actionButton:load_shapefile_button\n"), 
      actionButton(
        "load_shapefile_button", 
        "click to load the shapefile"
        ),
      
      cat("ui:actionButton:make_w_matrix_button\n"), 
      actionButton(
        "make_w_matrix_button", 
        "click to generate the neighbourhood matrix"
        ),
      
      ########################################################################
      hr(),
      h1("Model tweaking"),
      
      br(),
      cat("ui:sliderInput:posterior_sample_size\n"), 
      sliderInput("posterior_sample_size", 
                  "choose posterior sample size",
      min=1000, max=10000, step=1000, value=1000
      ),
      

      
      br(),
      cat("ui:sliderInput:select_burnin\n"), 
      sliderInput(
        "select_burnin", "Choose burn-in size",
        min = 1000, max = 100000, value = 20000, step = 2000
      ),

      br(),
      cat("ui:sliderInput:select_sample\n"), 
      sliderInput(
        "select_sample", "Choose sample size (must be larger than burn-in size)",
        min = 20000, max = 1000000, value = 100000, step = 2000
      ),
      
            
      br(),
      cat("ui:actionButton:generate_posterior_button\n"), 
      actionButton("generate_posterior_button", "click to run model"),
      
      br(),
      cat("ui:sliderInput:seg_k\n"), 
      sliderInput(
        "seg_k", 
        "Choose segregation thresholds",
        min=0, max=1, value=c(0,1)
      )
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
      
      #########################################################################
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
)


shinyApp(ui = ui, server = server)