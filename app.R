source("Functions_poisson.R")
source("server_function.R")
source("ui_function.R")
source("packages.R")


# UI ----------------------------------------------------------------------
ui = dashboardPage(
  dashboardHeader(title = "SSD for single-arm Poisson trials", titleWidth = 350),
  
  dashboardSidebar(disable = T),
  
  dashboardBody(
    tags$head(tags$style(HTML('
        .skin-blue .main-header .logo {
          background-color: #3c8dbc;
        }
        .skin-blue .main-header .logo:hover {
          background-color: #3c8dbc;
        }
      '))),
    tags$head(tags$style(HTML('
         .skin-blue .left-side, .skin-blue .wrapper {
                        background-color: #ecf0f5;
                        }
         '))),
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
      }
    "))),
    
    fluidRow(
      box(title = "General Setting", width = 4, ui.general.setting(), solidHeader = T),
      box(title = "Analysis Framework", width = 4, ui.test.setting.din(), solidHeader = T),
      box(title = "Power Approach", width = 4, ui.design.din(), solidHeader = T),
      ),
    
    tabBox(
    id = "tabs",
    tab_results(),
    tab_analysis(main = "Analysis prior", id = "analysis"),
    tab_design(main = "Design prior", id = "design"),
    selected = "res_panel",
    width = 12
    ),
    uiOutput("results.ui")
    )
)


# server ------------------------------------------------------------------


server = function(input, output, session){
  
  ## Dynamic ui definition
  dynamic_ui(input, output, session)
  
  ###### RESULTS
  general_results(input, output, session)
  ###### ANALYSIS PRIOR PANEL
  analysis_panel(input, output, session)
  
  # Design prior
  design_panel(input, output, session)
  
  #### Dependencies to deal with conditional ui:
  param = reactiveValues(alphaD_r = NULL, betaD_r = NULL,
                         thetaD_r = NULL,
                         alphaA_r = NULL, betaA_r = NULL, lambda_r = NULL)

  modify = reactiveVal(0)
  observeEvent(list(input$theta0, input$alternative),{
    modify(modify() + 1)
  })

  observeEvent(modify(), {
    updateNumericInput(inputId = "thetaD", value = NA)
    updateNumericInput(inputId = "thetaD_des", value = NA)
    updateNumericInput(inputId = "alphaD", value =  NA)
    updateNumericInput(inputId = "betaD", value = NA)
    updateNumericInput(inputId = "alphaA", value =  NA)
    updateNumericInput(inputId = "betaA", value = NA)
    param$alphaD_r = NULL
    param$betaD_r = NULL
    param$alphaA_r = NULL
    param$betaA_r = NULL
    param$thetaD_r = NULL
  })

  observeEvent(input$alphaD,{
               param$alphaD_r = input$alphaD})
  observeEvent(input$betaD,{
               param$betaD_r = input$betaD})
  observeEvent(input$alphaA,{
               param$alphaA_r = input$alphaA})
  observeEvent(input$betaA,{
               param$betaA_r = input$betaA})
  observeEvent(input$lambda,{
               param$lambda_r = input$lambda})

  observeEvent(input$thetaD,{
    updateNumericInput(inputId = "thetaD_des", value = input$thetaD)
    param$thetaD_r = input$thetaD
  })

  observeEvent(input$frame == "B",{
    updateNumericInput(inputId = "alphaA", value = param$alphaA_r)
    updateNumericInput(inputId = "betaA", value = param$betaA_r)
    updateSliderInput(inputId = "lambda", value = param$lambda_r)
  })

  observeEvent(input$approach == "P",{
    updateNumericInput(inputId = "alphaD", value = param$alphaD_r)
    updateNumericInput(inputId = "betaD", value = param$betaD_r)
  })

  observeEvent(input$approach == "C",{
    updateNumericInput(inputId = "thetaD", value = param$thetaD_r)
  })
  
  ## Modification of thetaD
  
  observeEvent(c(input$alternative == "less", input$theta0),{
    updateNumericInput(session = session, inputId = "thetaD", 
                       max = input$theta0)
    updateNumericInput(session = session, inputId = "thetaD_des", 
                       max = input$theta0, value = NA)
  })
  
  observeEvent(c(input$alternative == "greater", input$theta0),
               {
                 updateNumericInput(session = session, inputId = "thetaD",
                                    min = input$theta0)
                 updateNumericInput(session = session, inputId ="thetaD_des",
                                    min = input$theta0, value = NA)
               })
  
  ### WARNINGS
  
  # theta0 positive real number
  output$check1 = renderText({
    req(input$theta0)
    validate(need(input$theta0 > 0,"Please select a real positive number"))
  })
  
  # thetaD under the alternative hypothesis
  
  output$check2 = renderText({
    req(input$theta0, input$thetaD)
    text = NULL
    if ((input$thetaD >= input$theta0 & input$alternative == "less") | (input$thetaD <= input$theta0 & input$alternative == "greater")) {
      text = "Please select a value under the alternative hypothesis.\n"
    }
    if(input$thetaD <= 0){
      text = paste(text, " Please select a positive real value.")
    }
    validate(text)
  })
  
}


# Run App -----------------------------------------------------------------


# Run the application 
shinyApp(ui, server)

