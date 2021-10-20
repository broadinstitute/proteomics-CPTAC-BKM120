#################################################################
## Filename: ui.r
## Created: March 19, 2016
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from the BKM study.
## This file defines the user interface of the Shiny-app.
#################################################################
library(shiny)

#########################################################
## Define UI
#########################################################
shinyUI(fluidPage(

    ## Application title
    titlePanel(HTML(TITLESTRING), windowTitle = 'Mundt et al. Cancer Research 2018'),
    
    HTML('<br>'),
    HTML('<b>Supplemental data:</b>'),

    fluidRow(

        ##################################################
        ## LEFT panel
        ##################################################
        column(4, wellPanel(

                      ## authentification
                      uiOutput("auth.user"),

                      ## input
                      uiOutput("ui.input")
                  )

               ),

        ##################################################
        ## RIGHT panel
        ##################################################
        column(8,
               plotOutput("plot")
               ) ## end column

    ) ## end fluiRow
))
