# UI - Input --------------------------------------------------------------

ui.general.setting <- function() {
  column(
    12,
    h5(
      HTML("Choose &theta;<sub>0</sub>:"),
      span(shiny::icon("info-circle"), id = "theta0_tip")
    ),
    tippy_this(
      elementId = "theta0_tip",
      tooltip = "<span style='font-size:13px;'> <span>True event rate under the null hypothesis or historic control",
      placement = "right", allowHTML = T
    ),
    numericInput(inputId = "theta0", label = NULL, min = 0.001, value = NULL, step = 0.1),
    textOutput("check1"),
    div(
      style = "display: inline-block;", width = "100px",
      radioButtons(
        inputId = "alternative", label = "Alternative hypothesis",
        choiceNames = list(HTML("&theta; &#60;  &theta;<sub>0:"), HTML("&theta; &#62; &theta;<sub>0:")),
        choiceValues = list("less", "greater"), selected = "less", inline = T
      )
    ),
    sliderInput(inputId = "lev.pot", label = div(HTML("Desired power level &gamma;")), 
                min = 0.7, max = 0.999, step = 0.01, value = 0.8),
    numericInput(inputId = "n.max", label = "Maximum sample size", value = 200, min = 1),
    checkboxGroupInput(
      inputId = "criterion", label = "Type of criterion",
      list("Conservative" = "Conservative", "Classic" = "Classic"), selected = "Conservative", inline = T
    ),
    checkboxInput("plot", "Plot power function behavior", value = T),
  )
}

ui.test.setting.din <- function(input, output, session) {
  column(
    12,
    selectInput("frame",
      label = NULL,
      choices = list("Frequentist" = "F", "Bayesian" = "B"),
      selected = 1
    ),
    uiOutput("test")
  )
}

ui.design.din <- function(input, output, session) {
  column(
    12,
    selectInput("approach",
      label = NULL,
      choices = list("Conditional" = "C", "Predictive" = "P"),
      selected = 1
    ),
    uiOutput("approach")
  )
}

# UI - tabs ---------------------------------------------------------------
tab_results <- function() {
  tabPanel(
    title = "Results",
    value = "res_panel",
    fluidRow(
      column(
        7, verbatimTextOutput("print_res"),
        actionButton("save", "Save Results")
      ),
      column(5, uiOutput("plot.ui")),
    )
  )
}


tab_analysis <- function(main, id) {
  tabPanel(
    title = main,
    value = id,
    sidebarLayout(
      sidebarPanel(
        numericInput("thetaA", div(HTML("Prior mode &theta;<sup>A")),
          min = 0, value = NULL, step = 0.1
        ),
        numericInput("nA", div(HTML("Prior sample size n<sup>A")),
          min = 0, value = 1, step = 1
        ),
        verbatimTextOutput("info_analysis"),
        actionButton("update_par_analysis", "Update prior parameters")
      ),
      mainPanel(plotOutput("analysis_plot"))
    )
  )
}

tab_design <- function(main, id) {
  tabPanel(
    title = main,
    value = id,
    sidebarLayout(
      sidebarPanel(
        numericInput("thetaD_des", div(HTML("Prior mode &theta;<sup>D")),
          min = 0, value = NULL, step = 0.1
        ),
        selectInput("des_type", "Assign a fixed probability to",
          choices = c("The alternative hypothesis", "An interval")
        ),
        uiOutput("nD_methods"),
        sliderInput("prob_alt_des",
          label = div(HTML("Probability level")),
          min = 0.95, max = 0.999, step = 0.001, value = 0.999
        ),
        verbatimTextOutput("info_design"),
        actionButton("update_par_design", "Update prior parameters")
      ),
      mainPanel(plotOutput("design_plot"))
    )
  )
}



# Conditional UI ----------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


freq.ui <- function(input, output, session) {
  tagList(
    sliderInput(
      inputId = "alpha", label = div(HTML("Type I error &alpha;")),
      min = 0.01, max = 0.2, step = 0.005, value = input$alpha %||% 0.05
    ),
  )
}

bayes.ui <- function(input, output, session) {
  tagList(
    sliderInput(
      inputId = "lambda", label = div(HTML("Bayesian significance 1 - &epsilon;")),
      min = 0.7, max = 0.999, value = 0.95, step = 0.005
    ),
    h5(HTML("Analysis Prior"), span(shiny::icon("info-circle"), id = "analysis_tip")),
    tippy_this(
      elementId = "analysis_tip", placement = "right",
      tooltip = "<span style='font-size:13px;'> <span>Bayesian prior distribution.
               It can be chosen through the Analysis prior panel."
    ),
    numericInput(
      inputId = "alphaA", label = div(HTML("Shape parameter &alpha;<sup>A")),
      min = 0, value = NULL, step = 0.1
    ),
    numericInput(
      inputId = "betaA", label = div(HTML("Rate parameter &beta;<sup>A")),
      min = 0, value = NULL, step = 0.1
    ),
  )
}


cond.ui <- function(input, output, session) {
  tagList(
    h5(HTML("Design value &theta;<sup>D</sup>"), span(shiny::icon("info-circle"), id = "designv_tip")),
    tippy_this(
      elementId = "designv_tip", placement = "right",
      tooltip = "<span style='font-size:13px;'> <span>
               Clinically relevant value under the alternative hypothesis."
    ),
    numericInput(inputId = "thetaD", label = NULL, min = 0, value = NULL, step = 0.1),
    textOutput("check2"),
  )
}

pred.ui <- function(input, output, session) {
  tagList(
    h5(HTML("Design prior"), span(shiny::icon("info-circle"), id = "design_tip")),
    tippy_this(
      elementId = "design_tip", placement = "right",
      tooltip = " <span style='font-size:13px;'> <span>
               It should assign neglectable probability to the alternative hypothesis.
               It is advisable to choose it through the Design prior panel."
    ),
    numericInput(
      inputId = "alphaD", label = div(HTML("Shape parameter &alpha;<sup>D")),
      min = 0, value = NULL, step = 0.1
    ),
    numericInput(
      inputId = "betaD", label = div(HTML("Rate parameter &beta;<sup>D")),
      min = 0, value = NULL, step = 0.1
    ),
  )
}


nD.int.ui <- function(input, outpit, session) {
  tagList(
    numericInput("delta",
      label = div(HTML("Select the interval width &delta;")),
      min = 0, step = 0.1, value = NULL
    ),
    textOutput("check_delta"),
  )
}

results.ui <- function(input, output, session) {
  box(
    width = 12,
    fluidRow(
      column(
        2, actionButton("clear", "Clear"),
        downloadButton("download", "Download")
      ),
      column(10, tableOutput("results_tab"))
    )
  )
}
