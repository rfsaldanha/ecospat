# Packages
library(shinydashboard)
library(rgdal)
library(leaflet)
library(spdep)

ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "Spatial Econometrics"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "data", icon = icon("th")),
      menuItem("EDA", tabName = "eda", icon = icon("map")),
      menuItem("Models", tabName = "models", icon = icon("signal"))
    )
  ),
  dashboardBody(
    tabItems(
      # Tab Data
      tabItem(tabName = "data",
              h2("Data file"),
              fluidRow(
                box(solidHeader = TRUE,
                    "Load here your spatial data file including variables. Currently, we suggest the GML format.",
                    fileInput("input_data", label = "", multiple = FALSE)
                )
              )
      ),
      
      # Tab EDA
      tabItem(tabName = "eda",
              h2("Exploratory Data Analysis"),
              fluidRow(
                box(solidHeader = TRUE, width = 12,
                    leafletOutput("EDA_map", height = 600)
                )
              )
      ),
      tabItem(tabName = "models",
              h2("Models"),
              fluidRow(
                box(solidHeader = TRUE, width = 4,
                    title = "Options",
                    uiOutput("model_dependent_variable"),
                    uiOutput("model_independent_variable"),
                    sliderInput("models_k", label = "k-Nearest Neighbor Weights", min = 1, max = 50, value = 5, step = 1),
                    actionButton("models_estimate", "Estimate", icon = icon("math"), status = "primary")
                ),
                tabBox(width = 8,
                    tabPanel("OLS",
                      verbatimTextOutput("model_ols_summary"),
                      verbatimTextOutput("model_ols_error")
                    ),
                    tabPanel("SAR",
                             verbatimTextOutput("model_sar_summary"),
                             verbatimTextOutput("model_sar_impacts")
                    ),
                    tabPanel("SEM",
                             verbatimTextOutput("model_sem_summary")
                    ),
                    tabPanel("SAC",
                             verbatimTextOutput("model_sac_summary"),
                             verbatimTextOutput("model_sac_impacts")
                    ),
                    tabPanel("SLX",
                             verbatimTextOutput("model_slx_summary"),
                             verbatimTextOutput("model_slx_impacts")
                    ),
                    tabPanel("SDM",
                             verbatimTextOutput("model_sdm_summary"),
                             verbatimTextOutput("model_sdm_impacts")
                    ),
                    tabPanel("SDEM",
                             verbatimTextOutput("model_sdem_summary")
                    )
                )
              )
      )
    )
  )
)

server <- function(input, output) {
  ### Read spatial data file
  geodata <- reactive({
    req(input$input_data$datapath)
    readOGR(input$input_data$datapath)
  })
  
  output$EDA_map <- renderLeaflet({
    leaflet(geodata()) %>%
      addPolygons(
        stroke = TRUE, 
        fillOpacity = .5, 
        smoothFactor = 0.2
      ) %>%
      addTiles()
  })
  
  
  output$model_dependent_variable <- renderUI({
    variables <- names(geodata()@data)
    selectInput("model_dependent_variable", label = "Dependent variable", choices = variables)
  })
  
  output$model_independent_variable <- renderUI({
    variables <- names(geodata()@data)
    selectInput("model_independent_variable", label = "Independent variables", choices = variables, multiple = TRUE)
  })
  
  # Model specification
  esp <- reactive({
    paste0(as.character(input$model_dependent_variable), " ~ ", paste0(input$model_independent_variable, collapse = " + "))
  })
  
  # Weights
  w_matrix <- eventReactive(input$models_estimate, {
    coords <- coordinates(geodata())
    IDs <- row.names(geodata()@data)
    #nb2listw(knn2nb(knearneigh(coords, k=input$models_k),row.names=IDs),style="B")
    nb2listw(knn2nb(knearneigh(coords, k=input$models_k),row.names=IDs),style="W")
  })
  
  w_matrix_tr <- eventReactive(input$models_estimate, {
    W <- as(w_matrix(), "CsparseMatrix")
    trW(W, type="mult")
  })
  
  # OLS Model
  model_ols <- eventReactive(input$models_estimate, {
    lm(formula = formula(esp()), data = geodata()@data)
  })
  
  output$model_ols_summary <- renderPrint({
    summary(model_ols())
  })
  
  output$model_ols_error <- renderPrint({
    lm.morantest(model_ols(), w_matrix())
  })
  
  # SAR Model
  model_sar <- eventReactive(input$models_estimate, {
    lagsarlm(formula(esp()), data = geodata()@data, listw = w_matrix())
  })
  
  output$model_sar_summary <- renderPrint({
    summary(model_sar())
  })
  
  output$model_sar_impacts <- renderPrint({
    summary(impacts(model_sar(), tr=w_matrix_tr(), R=1000), zstats=TRUE, short=TRUE)
  })
  
  # SEM Model
  model_sem <- eventReactive(input$models_estimate, {
    errorsarlm(formula = formula(esp()), data = geodata()@data, listw = w_matrix())
  })
  
  output$model_sem_summary <- renderPrint({
    summary(model_sem())
  })
 
  # SAC Model
  model_sac <- eventReactive(input$models_estimate, {
    sacsarlm(formula = formula(esp()), data = geodata()@data, listw = w_matrix()) 
  })
  
  output$model_sac_summary <- renderPrint({
    summary(model_sac())
  })
  
  output$model_sac_impacts <- renderPrint({
    summary(impacts(model_sac(), tr=w_matrix_tr(), R=1000), zstats=TRUE, short=TRUE)
  })
  
  # SLX Model
  model_slx <- eventReactive(input$models_estimate, {
    lmSLX(formula = formula(esp()), data = geodata()@data, listw = w_matrix())
  })
  
  output$model_slx_summary <- renderPrint({
    summary(model_slx())
  })
  
  output$model_slx_impacts <- renderPrint({
    summary(impacts(model_slx(), tr=w_matrix_tr(), R=1000), zstats=TRUE, short=TRUE)
  })
  
  # SDM Model
  model_sdm <- eventReactive(input$models_estimate, {
    lagsarlm(formula = formula(esp()), data = geodata()@data, listw = w_matrix(), type = "mixed")
  })
  
  output$model_sdm_summary <- renderPrint({
    summary(model_sdm())
  })
  
  output$model_sdm_impacts <- renderPrint({
    summary(impacts(model_sdm(), tr=w_matrix_tr(), R=1000), zstats=TRUE, short=TRUE)
  })
  
  # SDEM Model
  model_sdem <- eventReactive(input$models_estimate, {
    errorsarlm(formula = formula(esp()), data = geodata()@data, listw = w_matrix(), etype = "emixed")
  })
  
  output$model_sdem_summary <- renderPrint({
    summary(model_sdem())
  })
  
}

shinyApp(ui, server)