# This is a Shiny web application for meta-analysis in SyRF
# by Qianying Wang @CAMARADES
# ui.R
# https://camarades.shinyapps.io/meta-analysis-syrf/

library(shiny)
library(metafor)
library(devtools)
# devtools::install_version("htmltools", version = "0.3.6", repos = "http://cran.us.r-project.org")
library(htmltools)

# install.packages("devtools")
# devtools::install_github("guido-s/meta")
library(meta)
library(shinythemes)
library(dplyr)
library(plotly)
library(RCurl)
library(ggplot2)
library(colourpicker)

shinyUI(fluidPage(
  
  theme = shinytheme("spacelab"),
  
  headerPanel(
    list(HTML('<img src="logo.jpg" height=90 width=90 /> '), "Meta-Analysis"),
    windowTitle="SyRFshinyapp"
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "file", label = "Upload CSV File",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
      ),
      
      radioButtons(inputId = "CompareType", label = "Select comparison",
                   choices = list("Analysis of model" = "SC", 
                                  "Analysis of intervention" = "SCT"), 
                   selected = "SC"),
      hr(),  
      radioButtons(inputId = "EffectMeasure", label = "Select effect size measure",
                   choices = list("Normalised mean difference" = "NMD", 
                                  "Standardised mean difference" = "SMD",
                                  "Odds ratio" = "OR"), 
                   selected = "NMD"),
      
      hr(),
      radioButtons(inputId = "HetMethod", label = "Select heterogeneity analysis method",
                   choices = list("Stratified meta-analysis" = "sub", 
                                  "Meta-regression" = "reg"),
                   selected = "sub"),
      
      hr(),
      selectInput(inputId = "HetEstimator", label = "Select heterogeneity estimator",
                  choices = list("DerSimonian-Laird" = "DL", 
                                 "Hedges" = "HE", 
                                 "Hunter-Schmidt" = "HS", 
                                 "Sidik-Jonkman" = "SJ", 
                                 "Maximum-likelihood" = "ML", 
                                 "Restricted maximum-likelihood" = "REML", 
                                 "Empirical Bayes" = "EB"),
                  selected = "REML"),
      
      hr(),
      checkboxInput(inputId = "KnHaTest", label = "Fit with Knapp and Hartung method", value = F),
      # radioButtons(inputId = "KnHaTest", label = "Fit with Knapp and Hartung method",
      #              choices = list("Yes" = T, 
      #                             "No" = F),
      #              selected = F),
      
      hr(),
      uiOutput("CheckVar")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data", 
                 verticalLayout(
                   wellPanel(
                     fluidRow(
                       column(6,
                              radioButtons(inputId = "DataType", label = "Select data you want to use for analysis",
                                           choices = list("Pre-nested data" = "pre", 
                                                          "Data nested by outcome" = "nest"
                                           ), 
                                           selected = "pre")
                       ),
                       column(6,
                              radioButtons(inputId = "FileType", label = "File type",
                                           choices = c("csv", "tsv")),
                              downloadButton('DownTable', 'Download table')
                       )
                     ) # fluidRow
                   ),
                   conditionalPanel(condition = "input.DataType == 'nest'",
                                    wellPanel(
                                      helpText("Select variables for data nesting:"),
                                      uiOutput("NestVar")
                                    )
                   ),
                   HTML("<br><br><br>"),
                   DT::dataTableOutput("DT")
                 )
        ),
        
        
        tabPanel("Meta-analysis",
                 verticalLayout(
                   verbatimTextOutput("GlobalOutput"),
                   wellPanel(
                     fluidRow(
                       column(6,
                              textInput(inputId = "ForXlab", label = "Label of x-axis", value = "Effect & REML"),
                              sliderInput(inputId = "GapLeft", label = "Left gap", min = 0, max = 10, value = 1),
                              sliderInput(inputId = "ForWinWidth", label = "Figure width (px)", min = 500, max = 2000, value = 600),
                              sliderInput(inputId = "ForWidth", label = "Plot width", min = 4, max = 50, value = 8)
                       ),
                       column(6,
                              fluidRow(
                                column(4,
                                       selectInput(inputId = "ForOrder", label = "Order of studies",
                                                   choices = list("original order" = "",
                                                                  "increasing effect size" = "ine",
                                                                  "increasing weight" = "inw",
                                                                  "increasing year" = "iny",
                                                                  "decreasing effect size" = "dee",
                                                                  "decreasing weight" = "dew",
                                                                  "decreasing year" = "dey"),
                                                   selected = "")
                                ),
                                column(4,
                                       colourInput(inputId = "ForSqCol", label = "Square color", value = "#E69F00")
                                ),
                                column(4,
                                       colourInput(inputId = "ForDiaCol", label = "Diamond color", value = "#56B4E9")
                                )
                              ),
                              
                              sliderInput(inputId = "GapRight", label = "Right gap", min = 0, max = 10, value = 1),
                              numericInput(inputId = "ForWinHeight", label = "Figure height (px)", value = 400, min = NA, max = NA, step = NA),
                              # sliderInput(inputId = "ForWinHeight", label = "Figure height (px)", min = 300, max = 6000, value = 400),
                              
                              HTML("<br><br>"),
                              
                              fluidRow(
                                column(6,
                                       checkboxInput("ShowWeight", label = "Show weight", value = T)
                                ),
                                column(6,
                                       downloadButton("DownForest", "Download forest plot")
                                )
                              )
                       )
                     )
                   ),
                   plotOutput("ForestPlot")
                 )
        ),
        
        tabPanel("Heterogeneity",
                 verticalLayout(
                   wellPanel(
                     conditionalPanel(condition = "input.HetMethod == 'sub'",
                                      uiOutput("SubVar")
                     ),
                     conditionalPanel(condition = "input.HetMethod == 'reg'",
                                      fluidRow(
                                        column(6, uiOutput("RegContVar")),
                                        column(6, uiOutput("RegDiscVar"))
                                      )
                     )
                   ),
                   # stratified meta-analysis
                   conditionalPanel(condition = "input.HetMethod == 'sub'",
                                    verbatimTextOutput("SubOutput") ),
                   # meta-regression
                   conditionalPanel(condition = "input.HetMethod == 'reg'",
                                    verbatimTextOutput("RegOutput") )
                 )
        ),
        
        
        
        tabPanel("Heterogeneity plot",
                 verticalLayout(
                   # subgroup forest plot
                   conditionalPanel(condition = "input.HetMethod == 'sub'",
                                    wellPanel(
                                      fluidRow(
                                        column(6,
                                               uiOutput("SubPlotVar"),
                                               sliderInput(inputId = "SubGapLeft", label = "Left gap", min = 0, max = 10, value = 1),
                                               sliderInput(inputId = "SubWinWidth", label = "Width (px)", min = 400, max = 1500, value = 600),
                                               sliderInput(inputId = "SubWidth", label = "Plot width (cm)", min = 4, max = 50, value = 8)
                                        ),
                                        column(6,
                                               fluidRow(
                                                 column(4,
                                                        textInput(inputId = "SubXlab", label = "Label of x-axis", value = "Effect & REML")
                                                 ),
                                                 column(4,
                                                        colourInput(inputId = "SubSqCol", label = "Square color", value = "#E69F00")
                                                 ),
                                                 column(4,
                                                        colourInput(inputId = "SubDiaCol", label = "Diamond color", value = "#56B4E9")
                                                 )
                                               ),
                                               
                                               sliderInput(inputId = "SubGapRight", label = "Right gap", min = 0, max = 10, value = 1),
                                               numericInput(inputId = "SubWinHeight", label = "Height (px)", value = 400, min = NA, max = NA, step = NA),
                                               # sliderInput(inputId = "SubWinHeight", label = "Height (px)", min = 400, max = 6000, value = 400),
                                               HTML("<br><br>"),
                                               fluidRow(
                                                 column(3,
                                                        checkboxInput(inputId = "ShowFixWeight", label = "Show w-fixed", value = T)
                                                 ),
                                                 column(3,
                                                        checkboxInput(inputId = "ShowRandWeight", label = "Show w-random", value = T)
                                                 ),
                                                 column(6,
                                                        downloadButton(outputId = "DownSubForest", label = "Download forest plot")
                                                 )
                                               )
                                               
                                        )
                                      ) # fluidRow
                                    ) # wellPanel
                   ), # conditionPanel
                   
                   
                   HTML("<br>"),
                   # meta-regression plot
                   conditionalPanel(condition = "input.HetMethod == 'reg'",
                                    wellPanel(
                                      fluidRow(
                                        column(6,
                                               uiOutput("RegPlotVar"),
                                               fluidRow(
                                                 column(6,
                                                        textInput(inputId = "RegXlab", label = "Label of x-axis", value = "type the xlabel")
                                                 ),
                                                 column(6,
                                                        textInput(inputId = "RegYlab", label = "Label of y-axis", value = "Effect size")
                                                 )
                                               ),
                                               fluidRow(
                                                 column(6,
                                                        colourInput(inputId = "RegPtCol", label = "Point color", value = "#131414")
                                                 ),
                                                 column(6,
                                                        textInput(inputId = "RegLCol", label = "Line color", value = "#FC130F")
                                                 )
                                               )
                                        ),
                                        column(6,
                                               sliderInput(inputId = "RegWidth", label = "Width (px)", min = 300, max = 2000, value = 700),
                                               HTML("<br><br>"),
                                               sliderInput(inputId = "RegHeight", label = "Height (px)", min = 300, max = 2000, value = 400)
                                        )
                                      ) # fluidRow
                                    ) # wellPanel
                   ), # conditionPanel
                   HTML("<br>"),
                   
                   conditionalPanel(condition = "input.HetMethod == 'sub'",
                                    plotOutput("SubForest") ),
                   conditionalPanel(condition = "input.HetMethod == 'reg'",
                                    plotlyOutput("RegPlot") )
                 )
        ),
        
        tabPanel("Bar plot",
                 verticalLayout(
                   wellPanel(
                     fluidRow(
                       column(6,
                              uiOutput("BarVar"),
                              sliderInput(inputId = "BarHeight", label = "Height (px)", min = 400, max = 2000, value = 400),
                              sliderInput(inputId = "BarTitleSize", label = "Title size", min = 8, max = 30, value = 12),
                              sliderInput(inputId = "BarYlabSize", label = "Y-axis label size ", min = 0, max = 30, value = 12),
                              sliderInput(inputId = "BarLabAngle", label = "Bar label angle", min = 0, max = 90, step = 1, value = 0)
                       ),
                       column(6,
                              fluidRow(
                                column(6,
                                       textInput(inputId = "BarYlab", label = "Label of y-axis", value = "Effect size")
                                ),
                                column(6,
                                       textInput(inputId = "BarTitle", label = "Title", value = "type the title")
                                )
                              ),
                              fluidRow(
                                column(6,
                                       textInput(inputId = "BarYmin", label = "y-min", value = "")
                                ),
                                column(6,
                                       textInput(inputId = "BarYmax", label = "y-max", value = "")
                                )
                              ),
                              sliderInput(inputId = "BarWidth", label = "Width (px)", min = 100, max = 2000, value = 400),
                              sliderInput(inputId = "BarLabSize", label = "Bar label size", min = 8, max = 30, value = 12),
                              sliderInput(inputId = "BarLabPos", label = "Bar label position", min = 0.1, max = 3, value = 0.5),
                              HTML("<br>"),
                              downloadButton("DownBar", "Download bar plot")
                       )
                     )
                   ),
                   tableOutput("BarInfoTable"),
                   plotOutput("BarPlot")
                 )
        ),
        
        
        
        tabPanel("Trim and Fill",
                 verticalLayout(
                   verbatimTextOutput("TafOutput"),
                   wellPanel(
                     fluidRow(
                       column(6,
                              radioButtons(inputId = "TafYaxis", label = "Select y-axis",
                                           choices = list("Inverse of the standard error" = "seinv",
                                                          "Inverse of the square-root sample size" = "sqrtninv"),
                                           selected = "seinv"),
                              
                              radioButtons(inputId = "TafFill", label = "Select funnel plot type",
                                           choices = list("Show imputed studies" = "Yes",
                                                          "Show only published studies" = "No"),
                                           selected = "Yes"),
                              
                              radioButtons(inputId = "TafSide", label = "Side of funnel plot the missing studies imputed",
                                           choices = list("Left" = "left",
                                                          "Right" = "right"),
                                           selected = "left"),
                              downloadButton("DownFunnel", "Download funnel plot")
                       ),
                       column(6,
                              fluidRow(
                                column(6,
                                       textInput(inputId = "TafXlab", label = "Label of x-axis", value = "Effect size")
                                ),
                                column(6,
                                       textInput(inputId = "TafYlab", label = "Label of y-axis", value = "type the y-axis label")
                                )
                              ),
                              sliderInput(inputId = "FunnelWidth", label = "Width (px)", min = 300, max = 1200, value = 600),
                              sliderInput(inputId = "FunnelHeight", label = "Height (px)", min = 300, max = 1200, value = 400)
                              # HTML("<br><br>"),
                              
                       )
                     )
                   ),
                   plotOutput("FunnelPlot")
                 )
        ),
        
        
        
        
        
        tabPanel("Egger's Regression",
                 wellPanel(
                   fluidRow(
                     column(6,
                            fluidRow(
                              column(6,
                                     textInput(inputId = "EggerXlab", label = "Label of x-axis", value = "Precision (1/SE)")
                              ),
                              column(6,
                                     textInput(inputId = "EggerYlab", label = "Label of y-axis", value = "Effect size/SE")
                              )
                            ),
                            HTML("<br>"),
                            fluidRow(
                              column(4,
                                     colourInput(inputId = "EggerPtCol", label = "Point color", value = "#131414")
                              ),
                              column(4,
                                     colourInput(inputId = "EggerLCol", label = "Line color", value = "#3F51B5")
                              ),
                              column(4,
                                     colourInput(inputId = "EggerRibCol", label = "Ribbon color", value = "#3F51B5")
                              )
                            )
                     ),
                     column(6,
                            sliderInput(inputId = "EggerWidth", label = "Width (px)", min = 300, max = 2000, value = 700),
                            sliderInput(inputId = "EggerHeight", label = "Height (px)", min = 300, max = 2000, value = 400)
                     )
                   ) # fluidRow
                 ), # wellPanel
                 plotlyOutput("EggerPlot"),
                 HTML("<br>"),
                 verbatimTextOutput("EggerOutput")
        ) # tabPanel("Egger's Regression"
      ) # tabsetPanel
    ) # mainPanel 
  ) # sidebarLayout
))
