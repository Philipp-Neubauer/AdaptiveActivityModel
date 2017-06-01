#FA shiny
#shiny::

library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggvis)

JScode <-
  "$(function() {
    setTimeout(function(){
      var vals = [0];
      var powStart = 5;
      var powStop = 1;
      for (i = powStart; i >= powStop; i--) {
        var val = Math.pow(10, i);
        val = parseFloat(val.toFixed(4));
        vals.push(val);
      }
      $('#m').data('ionRangeSlider').update({'values':vals})
    }, 5)})"

theme_default <- function() theme_bw()+theme(panel.grid=element_blank())

# Define UI for slider demo application
  ui = function(request) {
    dashboardPage(
     
    #  Application title
    header = dashboardHeader(title ="Adaptive activity model",titleWidth = 420),
    
    # Sidebar with sliders that demonstrate various available options
    sidebar = dashboardSidebar(width = 420,
                               sliderInput("m", 'Mass', min=1, max=10000, value=100,
                                           animate = animationOptions(interval=200),width = 350),
                               sliderInput("lm", 'dMass', min=100, max=10000, value=1000,
                                           animate = animationOptions(interval=200),width = 350),
                               sliderInput("n_int", 'dTemp', min=10, max=100, value=20,
                                           animate = animationOptions(interval=200),width = 350),
                               fluidRow(column(12,h3("Run Simulations"),offset = 3)),
                               fluidRow(column(1,actionButton("go", "Go"),offset = 5)),
                               fluidRow(column(1,h3(""),offset = 4)),
                               tabsetPanel(id="tabs",
                                           tabPanel(title = 'Trophic',tabName ='Foraging',
                                                    sliderInput("gamma", 'Consumption rate', min=1, max=100, value=40,
                                                                  animate = animationOptions(interval=200)),
                                                    sliderInput("h", 'Maximum ingestion', min=1, max=100, value=20,
                                                                  animate = animationOptions(interval=200)),
                                                    sliderInput("p", 'Consumption rate scaling', min=0.5, max=1, value=0.8,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("q", 'Maximum intake scaling', min=0.5, max=1, value=0.75,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("M", 'M', min=0.1, max=3,step = 0.01,
                                                                value=2,animate = animationOptions(interval=200)),
                                                    sliderInput("v", 'Mortality scaling', min=0, max=3, value=1,step = 0.1,
                                                                animate = animationOptions(interval=200))),
                                           
                                           tabPanel(title = 'Metabolism',
                                                    sliderInput("beta",label = 'SDA', min=0, max=0.5, value=0.15, 
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("phi", 'Excretion loss', min=0, max=0.5, value=0.25,
                                                                  animate = animationOptions(interval=200)),
                                                    sliderInput("delta", "Activity cost", min=0.2, max=10, value=4,
                                                                  animate = animationOptions(interval=200)),
                                                    sliderInput("k", "Std metabolism", min=0.2, max=10, value=1,
                                                                  animate = animationOptions(interval=200)),
                                                    sliderInput("n", 'Scaling', min=0.5, max=1, value=0.88,
                                                                animate = animationOptions(interval=200))),
                                           tabPanel(title = 'O2',tabName ='O2',
                                                    sliderInput("P50", 'P50', min=0, max=10, value=4.2,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("lO", 'Max O2', min=0.05, max=3, value=1.038,step=0.01,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("O2crit", 'Critical O2', min=0, max=10, value=2,step=0.2,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("omega", 'O2 consumption', min=0, max=0.1,step=0.01 ,value=0.03,
                                                                animate = animationOptions(interval=200))
                                           ),
                                           tabPanel(title = 'Temp',tabName ='Temp',
                                                    sliderInput("Ea", 'Activation energy', min=0, max=1,step = 0.02,
                                                                value=0.52,animate = animationOptions(interval=200)),
                                                    sliderInput("temp", 'Temp range (min to lethal)', min=0, max=32, value=c(8,26),
                                                                animate = animationOptions(interval=200),dragRange = T),
                                                    sliderInput("Topt", 'T max scope', min=0, max=32, value=20,step=0.5,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("shape", 'Scope shape', min=0, max=5, value=1,step=0.1,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("Tref", 'W\u221E plot temp', min=0, max=32, value=15,step=0.5,
                                                                animate = animationOptions(interval=200))
                                           )),
                               fluidRow(column(1,bookmarkButton(),offset = 4))),
    # Show a table summarizing the values entered
    dashboardBody(
      #plotOutput("plot"),
      #box(width = 12,height=1500,background = NULL,
               tags$table(
                 tags$td(ggvisOutput("ggvisTempplot")),
                 tags$td(ggvisOutput("ggvisO2plot"))
                 ),
               tags$table(
                 tags$td(ggvisOutput("ggvisMetplot")),
                 tags$td(ggvisOutput("ggvisfeedplot"))
               ),
               tags$table(
                 tags$td(ggvisOutput("ggvisEplot")),
                 tags$td(ggvisOutput("ggvisEffplot"))),
          fluidRow(column(10,tags$table(
                 tags$td(ggvisOutput("ggvisTGvis"))),offset=1))
               
             )
      #)
    
  )}
  
    