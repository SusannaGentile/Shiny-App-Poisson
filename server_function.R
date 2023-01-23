# App shiny: server function ---------------------------------------------------------------
plot.power.pois.shiny = function(res, theta0, power.type, criterion.type, lev.pot, n.max, hp)
{
  par(mai = c(1.02,0.6,0.3,0.6))
  n.min <- 2
  power <- res[[6]]
  actual.n <- length(power)
  n.plot.co <- res[[3]]
  n.plot.cl <- res[[1]]
  
  ## Setting max n.plot
  n.max.plot <- length(res[[6]])

  if(n.plot.co < n.max){
    if(n.plot.co > 1000)
      n.max.plot <- n.max.plot + 200
    else{
      if(theta0 > 1)
        n.max.plot <- n.plot.co + 20
      else{
        if(theta0 > 0.1)
          n.max.plot <- n.plot.co + 50
        else
          n.max.plot <- n.plot.co + 100
      }
    }
  }
  
  
  step <- 10
  if(n.plot.co >= 300) step <- 50
  if(n.plot.co >= 500) step <- 100
  
  y.label = switch(power.type, 
                   "FC" = "Frequentist conditional power", 
                   "FP" = "Frequentist predictive power",
                   "BC" = "Bayesian conditional power",
                   "BP" = "Bayesian predictive power")
  
  ## What to show if eta is always greater than 0.8
  if(n.plot.co == 0 & min(power) >= lev.pot){
    plot(power[1:(50 - n.min + 1)], x = n.min:50, type = "o", cex =0.5, 
         ylab = y.label, xlab = "n",  xlim = c(n.min-2,50), ylim = c(0,1))
    abline(h = lev.pot, lty = 2)
    return()
  }
  
  plot(y =  power[1:(n.max.plot-n.min+1)], x = n.min:n.max.plot, type = "o", 
       cex =0.5, ylab = y.label, xlab = "n",  xlim = c(n.min-3,n.max.plot), 
       ylim = c(0,1), xaxt = "n", las = 2)

  rm.x = seq(from = round(n.plot.co-10,-1), to = round(n.plot.co+10,-1), by = 10)
  ticks = seq(0,n.max.plot,by = step)[-(rm.x/step + 1)]
  
  axis(1, at = ticks, cex.axis = 0.8)
  abline(h = lev.pot, lty = 2)
  n.plot.cl = res[[1]]
  n.plot.co = res[[3]]
  
  if(actual.n > n.max){
    abline(v = n.max, lty = 1, col = "red")
    text(x = n.max + 8, y = 0.9, 
         labels = bquote(eta(n^{max})~"="~.(round(power[n.max-n.min-1],3))), cex = 0.8)
  }
  
  if("Classic" %in% criterion.type){
    power.cl =  power[n.plot.cl - n.min + 1]
    abline(v = n.plot.cl, lty = 3, lwd = 0.7)
    axis(1, at = n.plot.cl, cex.axis = 0.8, label = bquote(bold(.(as.character(n.plot.cl)))~"  "))
    points(x = n.plot.cl, y = power.cl, pch = 19, type = "o", cex = 0.6)
  }
  
  if("Conservative" %in% criterion.type){
    if(n.plot.cl == n.plot.co & "Classic" %in% criterion.type)
      return()
    power.co = power[n.plot.co - n.min + 1]
    abline(v = n.plot.co, lty = 3, lwd = 0.7)
    points(x = n.plot.co, y = power.co, pch = 19, type = "o", cex = 0.6)
    axis(1, at = n.plot.co,font = 2, cex.axis = 0.8, 
         label = bquote("  "~bold(.(as.character(n.plot.co)))))
  }
}

n.pois.shiny = function(theta0, hp, power.type = "FC",
                        thetaD = NULL,
                        alphaD = NULL, betaD = NULL,
                        alphaA = NULL, betaA = NULL, 
                        alpha, lambda, lev.pot, n.max, n.min = 2)
{
  dim.n = n.min:n.max
  power = eta.v(theta0 = theta0, hp = hp, thetaD = thetaD, alphaD = alphaD, 
                betaD = betaD, alphaA = alphaA, betaA = betaA, n = dim.n, 
                alpha = alpha, lambda = lambda, power.type = power.type)
  while(max(unlist(power[1,]))<lev.pot){
    n.max = 2*n.max
    dim.n = n.min:n.max
    power = eta.v(theta0 = theta0, hp = hp, thetaD = thetaD, 
                  alphaD = alphaD, betaD = betaD, 
                  alphaA = alphaA, betaA = betaA, 
                  n = n.min:n.max, alpha = alpha, lambda = lambda, 
                  power.type = power.type)
  }
  
  vett.prob = unlist(power[1,])
  r.vett = unlist(power[2,]) 
  results = n_min_opt(dim.n = dim.n, vett.prob = vett.prob, r.vett = r.vett, lev.pot = lev.pot)
  return(append(results, list(dim.n, vett.prob, r.vett)))
}



tab.pois = function(res, theta0, hp,  thetaD = NA, alphaA = NA, betaA = NA, 
                    alphaD = NA, betaD = NA, lev.pot = NA, criterion.type, 
                    alpha = NA, lambda = NA, frame = NA, approach = NA)
{
  tab.cl = NULL
  tab.co = NULL
  
  if(frame == "F")
    alphaA = lambda = betaA = NA
  else
    alpha = NA
  
  if(approach == "C")
    alphaD = betaD = NA
  else
    thetaD = NA
  power.type = paste(frame, approach, sep = "")
  if("Classic" %in% criterion.type){
    tab.cl = matrix(c(power.type, theta0, hp, 
                      thetaD, alphaD, betaD, 
                      alphaA, betaA, 
                      alpha, 1 - lambda, lev.pot, "Classic", 
                      res[[1]], res[[2]]), 
                      nrow = 1, dimnames = list(NULL, c("type.power", "theta0", "H1", 
                                                        "thetaD","alphaD","betaD",
                                                        "alphaA","betaA",
                                                        "alpha","epsilon",
                                                        "Power","Criterion","n", "k")))
  }
  
  if("Conservative" %in% criterion.type){
    tab.co = matrix(c(power.type,theta0, hp, thetaD, alphaD, betaD, alphaA, betaA, 
                      alpha,1 - lambda ,lev.pot, "Conservative",res[[3]], res[[4]]), nrow = 1, 
                    dimnames = list(NULL, c("type.power","theta0", "H1", "thetaD",
                                            "alphaD","betaD","alphaA","betaA",
                                            "alpha","epsilon","Power","Criterion","n","k")))
  }
  tab <-  rbind(tab.cl,tab.co)
  return(tab)
}


prior.plot = function(param, theta0, xlim){
  curve(dgamma(x, param[1], param[2]), xlim = xlim, lwd = 1.5, ylab = " ", xlab = expression(theta), n = 1000)
  abline(v = theta0, lty = 2)
}



# Print results -------------------------------------------------

print.result = function(n, theta0, hp, alphaA = NULL, betaA = NULL, 
                        thetaD = NULL, alphaD = NULL, betaD = NULL, 
                        criterion.type, lev.pot, alpha, lambda, frame, appr, n.max){
  
  frame.label = switch(frame, 
                       "F" = "Frequentist",
                       "B" = "Bayesian")
  appr.label = switch(appr,
                      "C" = "Conditional",
                      "P" = "Predictive")
  
  hp.system = switch(hp,
                     "less" = "True rate is less than ",
                     "greater" = "True rate is greater than ")
  
  cat("Alternative hypothesis        : ", hp.system , theta0, "\n", sep = "")
  cat("Power function                : ", frame.label, " ", appr.label, " ", "Power", "\n", sep = "")
  cat("\n")
  switch(frame, "F" = 
           cat("Type I error                  : ", alpha, "\n", sep = ""),
         "B" = {
           cat("Bayesian significance level   : ", lambda, "\n",
               "Analysis prior hyperparameters: ", "\u03B1\u1D2C", " = ", alphaA, "\n",
               "                                ", greeks("beta"), "\u1D2C", " = ", betaA, "\n"
               , sep = "")
         }
  )
  cat("\n")
  switch(appr, 
         "C" = cat("Design value ", greeks("theta"), "\u1D30               : ", thetaD, "\n" , sep = ""),
         "P" = cat("Design prior hyperparameters  : ", greeks("alpha"), "\u1D30", " = ", alphaD, "\n",
                   "                                ", greeks("beta"), "\u1D30", " = ", betaD, "\n", sep = "")
  )
  cat("\n")
  
  if(n[[3]] == 0 & min(n[[7]]) > lev.pot){
    cat("The power is greater than", lev.pot, "for every sample size")
    }
  else{
    if(n.max < n[[3]] | is.null( n[[3]])){
      cat("WARNING: the maximum sample size is too small to reach the desired power level \n")
  }
    else{
      if(is.null(criterion.type))
        cat("No SSD criterion selected")
      if("Conservative" %in% criterion.type){
        cat("The optimal sample size at level ", lev.pot, 
            " according to the conservative criterion is ", n[[3]], "\n", sep = "")
        cat("The corresponding critical value k is ", n[[4]], "\n", sep = "")
        cat("\n")
    }
      if("Classic"  %in% criterion.type){
        cat("The optimal sample size at level ", lev.pot, 
            " according to the classic criterion is ", n[[1]], "\n", sep = "")
        cat("The corresponding critical value k is ", n[[2]], sep="")
    }
  }
  }
}


# Conditional ui ----------------------------------------------------------

dynamic_ui = function(input, output, session){
  ## Test
  output$test <-  renderUI({
    if(input$frame == "B"){
      bayes.ui(input, output, session)
      }
    else{
      freq.ui(input, output, session)
      }
    })
  ## Approach
  output$approach <-  renderUI({
    if(input$approach == "P"){
      pred.ui(input, output, session)
      }
    else{
      cond.ui(input, output, session)
      }
    })
  ## Plot
  output$plot.ui <- renderUI({
    if(input$plot == TRUE)
    {
        plotOutput("plot")
    }

  })
}

# General results ---------------------------------------------------------

general_results <- function(input, output, session)
{
  power.type <- reactive(
    paste0(input$frame, input$approach, sep = "")
  )
  mode_des <- reactive((input$alphaD - 1)/input$betaD)
  
  n <- reactive({
    req(input$theta0, input$tabs == "res_panel", input$n.max)
    if(input$frame == "B"){
      req(input$alphaA, input$betaA)
    }
    if(input$approach == "C"){
      req(input$thetaD,
          if(input$alternative == "less") input$thetaD < input$theta0
          else input$thetaD > input$theta0)
    }
    else{
      req(input$alphaD, input$betaD,
          if(input$alternative == "less") mode_des() < input$theta0
          else mode_des() > input$theta0
      )
          }
    n.pois.shiny(theta0 = input$theta0, thetaD = input$thetaD, power.type = power.type(), 
                 alphaA = input$alphaA, betaA = input$betaA, 
                 alphaD = input$alphaD, betaD = input$betaD, 
                 lambda = input$lambda, alpha = input$alpha, 
                 lev.pot = input$lev.pot, n.max = input$n.max, hp = input$alternative)
  }) %>% debounce(500)
  
  
  output$plot <- renderPlot({
    req(n(), input$n.max)
    if(input$frame == "B"){
      req(input$alphaA, input$betaA)
    }
    if(input$approach == "C"){
      req(input$thetaD,
          if(input$alternative == "less") input$thetaD < input$theta0
          else input$thetaD > input$theta0)
    }
    else{
      req(input$alphaD, input$betaD)
    }
    plot.power.pois.shiny(n(), theta0 =  input$theta0, power.type = power.type(), 
                          criterion.type = input$criterion, lev.pot = input$lev.pot, 
                          n.max = input$n.max, hp = input$alternative)
    })
  
  output$print_res <- renderPrint({
    req(n(), input$n.max)
    if(input$frame == "B"){
      req(input$alphaA, input$betaA)
    }
    if(input$approach == "C"){
      req(input$thetaD,
          if(input$alternative == "less") input$thetaD < input$theta0
          else input$thetaD > input$theta0)
    }
    else{
      req(input$alphaD, input$betaD)
    }
    print.result(n(), theta0 = input$theta0, hp = input$alternative, 
                 alphaA = input$alphaA, betaA = input$betaA, thetaD = input$thetaD, 
                 alphaD = input$alphaD, betaD = input$betaD, criterion.type = input$criterion, 
                 lev.pot = input$lev.pot, alpha = input$alpha, lambda = input$lambda, 
                 frame = input$frame, appr = input$approach, n.max = input$n.max)
    })
  #### General results
  tab <- reactiveVal()
  count <- reactiveVal(0)
  
  observeEvent(input$save,{
    req(input$criterion)
    count(count() + 1)
    tab(rbind(tab(),tab.pois(res = n(), theta0 = input$theta0, hp = input$alternative, 
                             thetaD = input$thetaD, alphaD = input$alphaD, betaD =input$betaD, 
                             alphaA =input$alphaA, betaA = input$betaA, criterion.type = input$criterion, 
                             approach = input$approach, frame = input$frame, alpha = input$alpha, 
                             lev.pot = input$lev.pot, lambda = input$lambda)))
    })
  
  output$results.ui <- renderUI({
    if(count() > 0)
      results.ui(input, output, session)
  })
  observeEvent(input$clear,{
    count(0)
    tab(NULL)})
  
  output$results_tab <-  renderTable({
    req(!is.null(tab()))
    tab()
  })
  
  output$download <- downloadHandler(
    filename = function(){"SSD.csv"}, 
    content = function(fname){
      write.table(tab(), fname, na = "NA", quote = c(1,11))
    }
  )
  
}

# Analysis prior panel ----------------------------------------------------

analysis_panel <- function(input, output, session)
{
  # SHOW/HIDE
observeEvent(input$frame, {
  if(input$frame == "F")
    hideTab(inputId = "tabs", target = "analysis", session = session)
  else
    showTab(inputId = "tabs", target = "analysis", session = session)
  }) 

  #### Hyperparametters computation
  prior_analysis <-  reactive({
      req(input$thetaA, input$nA, input$tabs == "analysis")
      prior.pois(input$thetaA, input$nA)
  })

  output$info_analysis <- renderPrint({
    cat("\u03B1\u1D2C", " = ", prior_analysis()[[1]], "\n",greeks("beta"), "\u1D2C", " = ",
        prior_analysis()[[2]], "\n", sep = "")
    cat(paste("Probability of the alternative hypothesis = ", round(
      if(input$alternative == "less")
        round(pgamma(input$theta0, prior_analysis()[[1]], prior_analysis()[[2]]),4)
      else
        1 -pgamma(input$theta0, prior_analysis()[[1]], prior_analysis()[[2]]),4)))
    })
  
  output$analysis_plot <-  renderPlot({
    req(prior_analysis(), input$theta0)
    prior.plot(prior_analysis(), input$theta0, 
               xlim = c(0, 
                        qgamma(1 - 10^-4, shape = round(prior_analysis()[[1]]),  
                               rate = round(prior_analysis()[[2]]))))
    })
### Parameters update
  observeEvent(input$update_par_analysis,
               {
                 updateNumericInput(inputId = "alphaA", value = prior_analysis()[[1]])
                 updateNumericInput(inputId = "betaA", value = prior_analysis()[[2]])
                 })

}

# Design prior panel ------------------------------------------------------


design_panel <- function(input, output, session)
{
  # SHOW/HIDE
  observeEvent(input$approach, {
    if(input$approach == "C")
      hideTab(inputId = "tabs", target = "design", session = session)
    else
      showTab(inputId = "tabs", target = "design", session = session)
  }) 
  
  ###### Prior sample size selection
  output$nD_methods <- renderUI({
    if(input$des_type == "An interval")
      nD.int.ui(input, output, session)
  })

  ####### Delta selection
  max_delta <- reactive({
    req(input$thetaD_des, input$theta0)
    if(input$alternative == "less"){
      min(input$thetaD_des, input$theta0 - input$thetaD_des)
    }
    else
      input$thetaD - input$theta0
  })
  observeEvent(input$thetaD_des,
               {
                 req(max_delta())
                 updateNumericInput(inputId = "delta", max = max_delta())
               })
  output$check_delta <- renderText({
    req(input$delta, max_delta())
    text = paste(greeks("delta"), "can be at most equal to", max_delta())
    validate(need(input$delta <=  max_delta(), text))
  })

  ### hyperparameters computation
  interval <- reactive({
    req(input$thetaD_des, input$delta)
    list(input$thetaD_des - input$delta, input$thetaD_des + input$delta)
  })
  
  prior_design <-  reactive({
    req(input$thetaD_des, input$theta0, input$tabs == "design")
    if(input$alternative == "less")
      req(input$thetaD_des < input$theta0)
    else
      req(input$thetaD_des > input$theta0)
    
      if(input$des_type == "The alternative hypothesis")
      {
        if(input$alternative == "less"){
          round(hyper.gamma.int(input$thetaD_des, input$prob_alt_des, l.inf = 0, l.sup = input$theta0, upper = 10000),4)
        }
        else
          round(hyper.gamma.int(input$thetaD_des, input$prob_alt_des, l.inf = input$theta0, l.sup = +Inf, upper = 10000),4)
      }
      else
      {
        req(interval())
        round(hyper.gamma.int(input$thetaD_des, input$prob_alt_des, l.inf = interval()[[1]], l.sup = interval()[[2]], upper = 10000),4)
      }
    }) 
  output$info_design <- renderPrint({
    req(prior_design(),input$thetaD_des)
    if(input$des_type == "An interval")
        cat("The interval is [",input$thetaD_des - input$delta , ", ",  input$thetaD_des + input$delta, "]\n", sep = "")
    
    cat(greeks("alpha"), "\u1D30", " = ", prior_design()[[1]], "\n", greeks("beta"), "\u1D30", " = ", prior_design()[[2]], "\n", sep = "")

  })
  
  output$design_plot <-  renderPlot({
    req(prior_design(), input$thetaD_des)
    if(input$alternative == "less")
      prior.plot(prior_design(), input$theta0, xlim = c(max(0, input$theta0/2 - 1),input$theta0 + 0.2))
    else
      prior.plot(prior_design(), input$theta0, xlim = c(input$theta0, max(input$theta0 + 0.5, input$theta0*3)))
    })
  ### Parameters update
  observeEvent(input$update_par_design,
               {
                 updateNumericInput(inputId = "alphaD", value = prior_design()[[1]])
                 updateNumericInput(inputId = "betaD", value = prior_design()[[2]])
               })
  
}




