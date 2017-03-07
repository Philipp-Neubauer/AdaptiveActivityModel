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
                               
                               selectInput("variable", "Variable:",
                                           c("Mass" = "m",
                                             "Activity cost" = "delta",
                                             "Std metabolism" = "k",
                                             'Consumption rate' = 'gamma',
                                             'Assimilation loss'='phi',
                                             'Maximum ingestion'='h',
                                             'SDA'='beta'),selected = 'm'),
                               
                               fluidRow(column(1,''),column(1,'Min'),column(3,textInput('min',NULL,width = 80,value=1)),column(1,''),column(1,'Max'),column(4,textInput('max',NULL,width = 80,value = 1000))),
                               
                               
                               
                               conditionalPanel(
                                 condition = "input.variable != 'm'",
                                 tags$head(tags$script(HTML(JScode))),
                                 sliderInput("m", "Mass", min=1, max=10000, value=1000,step = 100,
                                             animate = animationOptions(interval=200))
                               ),
                               
                               tabsetPanel(id="tabs",
                                           tabPanel(title = 'Foraging',tabName ='Foraging',
                                                    conditionalPanel(
                                                      condition = "input.variable != 'gamma'",
                                                      sliderInput("gamma", 'Consumption rate', min=1, max=100, value=40,
                                                                  animate = animationOptions(interval=200))
                                                    ),
                                                    conditionalPanel(
                                                      condition = "input.variable != 'h'",
                                                      sliderInput("h", 'Maximum ingestion', min=1, max=100, value=20,
                                                                  animate = animationOptions(interval=200))
                                                    ),
                                                    sliderInput("p", 'Consumption rate scaling', min=0.5, max=1, value=0.8,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("q", 'Maximum intake scaling', min=0.5, max=1, value=0.75,
                                                                animate = animationOptions(interval=200))),
                                           tabPanel(title = 'Metabolism',
                                                    conditionalPanel(
                                                      condition = "input.variable != 'beta'",
                                                      sliderInput("beta",label = 'SDA', min=0, max=0.5, value=0.15, animate = animationOptions(interval=200))
                                                    ),
                                                    conditionalPanel(
                                                      condition = "input.variable != 'phi'",
                                                      sliderInput("phi", 'Excretion loss', min=0, max=0.5, value=0.25,
                                                                  animate = animationOptions(interval=200))
                                                    ),
                                                    conditionalPanel(
                                                      condition = "input.variable != 'delta'",
                                                      sliderInput("delta", "Activity cost", min=0.2, max=10, value=2,
                                                                  animate = animationOptions(interval=200))
                                                    ),
                                                    conditionalPanel(
                                                      condition = "input.variable != 'k'",
                                                      sliderInput("k", "Std metabolism", min=0.2, max=10, value=2,
                                                                  animate = animationOptions(interval=200))
                                                    ),
                                                    sliderInput("n", 'Scaling', min=0.5, max=1, value=0.8,
                                                                animate = animationOptions(interval=200))),
                                           tabPanel(title = 'O2',tabName ='O2',
                                                    uiOutput("o2slide"),
                                                    sliderInput("P50", 'P50', min=0, max=10, value=4.2,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("lO", 'Max O2', min=0.05, max=3, value=1.038,step=0.01,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("O2crit", 'Critical O2', min=0, max=10, value=2,step=0.2,
                                                                animate = animationOptions(interval=200)),
                                                    sliderInput("omega", 'O2 consumption', min=0, max=0.01, value=0.003,
                                                                animate = animationOptions(interval=200))
                                           ),
                                           tabPanel(title = 'Temp',tabName ='Temp',
                                                    sliderInput("Ea", 'Activation energy', min=0, max=1,step = 0.02,
                                                                value=0.52,animate = animationOptions(interval=200)),
                                                    sliderInput("temp", 'Temp range (min to lethal)', min=0, max=32, value=c(1,30),
                                                                animate = animationOptions(interval=200),dragRange = T),
                                                    sliderInput("Topt", 'Optimal Temp', min=0, max=32, value=15,
                                                                animate = animationOptions(interval=200))
                                           ),
                                           tabPanel(title = 'W\u221E',tabName ='Growth',
                                                    sliderInput("r", 'Reproduction allocation', min=0, max=0.6,step = 0.02,
                                                                value=0.2,animate = animationOptions(interval=200)),
                                                    sliderInput("temp_ref", 'Temperature', min=1, max=32,step = 1,
                                                                value=20,animate = animationOptions(interval=200))
                                           )),
                               fluidRow(column(1,bookmarkButton(),offset = 4))),
    # Show a table summarizing the values entered
    dashboardBody(
      #plotOutput("plot"),
      tabBox(width = 12,
             tabPanel('Activity & Scope',tags$table(
               tags$td(ggvisOutput("ggvisplot")),
               tags$td(ggvisOutput("ggvisO2plot"))),
               tags$table(
                 tags$td(ggvisOutput("ggvisTempplot")),
                 tags$td(ggvisOutput("ggvisMetplot")))
             ),
             tabPanel('Outcomes',
                      tags$table(
                        tags$td(ggvisOutput("ggvisfeedplot")),
                        tags$td(ggvisOutput("ggvisEplot"))),
                      tags$table(
                        tags$td(ggvisOutput("ggvisEffplot")),
                        tags$td(ggvisOutput("ggvispredplot"))
                      )
             ),
             tabPanel('W\u221E',
                      tags$td(ggvisOutput("ggvisGvis")),
                      tags$table(
                        tags$td(ggvisOutput("ggvisOGvis")),
                        tags$td(ggvisOutput("ggvisTGvis"))
                      )
             )
      )
    )
  )}
  
    