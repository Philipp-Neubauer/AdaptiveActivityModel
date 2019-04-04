#FA shiny
#shiny::

library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggvis)

options(shiny.sanitize.errors = F)

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
                               sliderInput("m", 'Mass for physio plots (g; log10)', min=0, max=5, value=2,width = 350),
                               sliderInput("mr", 'Mass range (g; log10)', min=-4, max=6, value=c(-2,3),width = 350,dragRange = T),
                               sliderInput("lm", 'dMass', min=10, max=10000, value=100,width = 350),
                               sliderInput("dt", 'dt', min=100, max=10000, value=100,width = 350),
                               sliderInput("n_int", 'dTemp', min=10, max=100, value=50,width = 350),
                               selectInput("def",'Default values',choices=c('M strategy' = 'M','P strategy' = 'P'),selected='M'),
                               fluidRow(column(12,h3("Run Simulations"),offset = 3)),
                               fluidRow(column(1,actionButton("go", "Go"),offset = 4)),
                               fluidRow(column(1,h3(""),offset = 4)),
                                          
                               tabsetPanel(id="tabs",
                                           tabPanel(title = 'RN',tabName ='RN',
                                                    sliderInput("tmax", 'Max time', min=0, max=100, value=25,width = 350),
                                                    sliderInput("slope", 'Reaction norm slope', min=-2, max=2, value=0,step=0.01,width = 350),
                                                    sliderInput("tr", 'Reaction', min=0, max=2, value=0.5,step=0.05,width = 350),
                                                    sliderInput("c", 'Env change', min=0, max=1, value=0.3,step=0.1,width = 350)
                                           ),
                                           tabPanel(title = 'Trophic',tabName ='Foraging',
                                                    sliderInput("gamma", 'Maximum encountered food', min=1, max=100, 
                                                                value=60,
                                                                  animate = animationOptions(interval=3000)),
                                                    sliderInput("h", 'Maximum ingestion', min=1, max=100, value=30,
                                                                  animate = animationOptions(interval=3000)),
                                                    sliderInput("p", 'Consumption scaling', min=0.5, max=1, value=0.8,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("q", 'Maximum intake scaling', min=0.5, max=1, value=0.8,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("M", 'M', min=0.01, max=1,step = 0.01,
                                                                value=0.1,animate = animationOptions(interval=3000)),
                                                    sliderInput("v", 'Mortality coeff', min=0, max=10, value=6,step = 0.1,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("nu", 'Mortality scaling', min=-1, max=0, value=-0.2,step = 0.05,
                                                                animate = animationOptions(interval=3000))),
                                           
                                           tabPanel(title = 'Metabolism',
                                                    sliderInput("beta",label = 'SDA', min=0, max=0.5, value=0.15, 
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("phi", 'Excretion loss', min=0, max=0.5, value=0.25,
                                                                  animate = animationOptions(interval=3000)),
                                                    sliderInput("delta", "Activity cost", min=0.2, max=10, value=4,
                                                                  animate = animationOptions(interval=3000)),
                                                    sliderInput("k", "Std metabolism", min=0.2, max=10, value=1,
                                                                  animate = animationOptions(interval=3000)),
                                                    sliderInput("n", 'Scaling', min=0.5, max=1, value=0.88,
                                                                animate = animationOptions(interval=3000))),
                                           tabPanel(title = 'O2',tabName ='O2',
                                                    sliderInput("P50", 'P50', min=0, max=10, value=4.2,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("lO", 'Max O2', min=0.05, max=5, value=0.5,step=0.1,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("O2crit", 'Critical O2', min=0, max=10, value=2,step=0.2,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("omega", 'O2 consumption', min=0, max=1,step=0.01 ,value=0.25,
                                                                animate = animationOptions(interval=3000))
                                           ),
                                           tabPanel(title = 'Temp',tabName ='Temp',
                                                    sliderInput("Ea", 'Activation energy', min=0, max=1,step = 0.02,
                                                                value=0.52,animate = animationOptions(interval=3000)),
                                                    sliderInput("temp", 'Temp range (min to lethal)', min=0, max=32, value=c(5,26),
                                                                animate = animationOptions(interval=3000),dragRange = T),
                                                    sliderInput("Topt", 'T max scope', min=0, max=32, value=20,step=0.5,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("shape", 'Scope shape', min=0, max=5, value=3,step=0.1,
                                                                animate = animationOptions(interval=3000)),
                                                    sliderInput("Tref", 'W\u221E plot temp', min=0, max=32, value=15,step=0.5,
                                                                animate = animationOptions(interval=3000))
                                           )),
                               fluidRow(column(1,bookmarkButton(),offset = 4))),
    # Show a table summarizing the values entered
    dashboardBody(
      tabBox(width = "800px",
      tabPanel('Rates',
               plotOutput("Tauplot"),
               plotOutput("O2plot"),
               plotOutput("Eplot")),
      tabPanel('Eco reaction',
               fluidRow(column(6,plotOutput("am")),column(6,plotOutput("alloc"))),
               fluidRow(column(6,plotOutput("ls")),column(6,plotOutput("norm")))),
      tabPanel('Evo reaction',
               fluidRow(column(6,plotOutput("TGvis")),column(6,plotOutput("Gvis"))),
               plotOutput("R0"))
             
    ))
    
    )}
  