#################################################################
## Filename: server.r
## Created: March 19, 2016
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data of the BKM120 project.
## This defines the server logic of the Shiny-app.
##
## changelog: 20160613 switched from md5 to bcrypt
#################################################################
library(shiny)

########################################################
## Define server logic
########################################################
shinyServer( function(input, output, session) {

    global <- reactiveValues(auth=T)

    ##############################
    ## text field for passphrase
    output$auth.user <- renderUI({

        if(global$auth) return()
        list(
            passwordInput('passphrase',label='Enter passphrase', width=120),
            actionButton('authbutton', 'GO')
        )
    })

    ##############################
    ## check passphrase
    observeEvent(input$authbutton, {
        global$auth <- authenticateUser(input$passphrase)
    })

    ##############################
    ## ui
    ##############################
    output$ui.input <- renderUI({

        if(!global$auth) return()

        list(
          ## text input
          textInput('genes', label=paste('Paste a list of gene names (max. ', GENEMAX,')', sep=''), value=GENESSTART),
          HTML('<br><br>'),
          ## submit button
          fluidRow(
              ##column(2, submitButton('GO')),
              column(3, radioButtons('zscore', label='Z-score', choices=c('row', 'none'), selected='row')),
              column(3, radioButtons('allsites', label='pSTY sites', choices=c('most variable', 'all'), selected='most variable')),
              column(3, textInput('min.val', label='min', value=-3, width='80%')),
              column(3, textInput('max.val', label='max', value=3, width='80%')),
              HTML('<br><br>'),

              ## download buttons
              fluidRow(
                  column(6, downloadButton('downloadHM', 'Download pdf')),
                  column(6, downloadButton('downloadTab', 'Download Excel'))
              ),
              HTML('<br><br>'),
              HTML('<p><b>Getting started</b></p>'),
              helpText('Enter or paste your gene names of interest (official gene symbols, e.g. PIK3R2) into the text field. The text field accepts lists of up to 20 gene symbols in either comma-, semicolon-, or space-separated format.'),
              HTML('<p>For more details please see our publication <a href="http://cancerres.aacrjournals.org/content/78/10/2732" target="_blank_">Mundt <i>et al.</i> Cancer Research. 2018</a></p>')
          ))
    })

    ##############################
    ## update list of input genes
    observeEvent(input$genes, {

        if(is.null(input$genes)) return()
        if(!global$auth) return()

        global$genes.input <- extractGenes(input$genes)
    })

    ##############################
    ## generate the heatmap
    output$plot <- renderPlot({

        if(!global$auth) return()
        if(is.null(input$genes)) return()

        genes.vec <- extractGenes( input$genes )

        if(length(genes.vec)==0) return()

        hm=makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, zscore=input$zscore, show.sites=input$allsites, min.val=input$min.val, max.val=input$max.val)
        global$expr.select <- hm
    },
    width = function(){ width=950},
    height= function(){ height=ifelse( global$auth, dynamicHeightHM(length( findGenesInDataset(extractGenes( input$genes ), input$allsites) ), length(unique(extractGenes( input$genes ))) ), 0 )}
    )

    #############################
    ## download heatmap
    output$downloadHM <- downloadHandler(

        filename = {
            paste( FILENAMESTRING, "_", length( extractGenes(input$genes) ), "_genes.pdf", sep='')
            ##fn.tmp
            ##if(input$zscore == 'row')
            ##   fn.tmp <- paste(fn.tmp, 'Zscore', sep='_')
            ##fn.tmp
            ## if(input$all.sites == 'all')
           ##     fn.tmp <- paste(fn.tmp, 'all_pSTY', sep='_')
           ## if(input$all.sites == 'most variable')
           ##     fn.tmp <- paste(fn.tmp, 'most_variable_pSTY', sep='_')

            ##paste(fn.tmp, '.pdf', sep='')

            ##paste( FILENAMESTRING, "_", length( extractGenes(input$genes) ), "_genes", paste(ifelse( input$zscore == 'row', '_Zscore_', '_')),
            ##      input$allsites, '_pSTY.pdf', sep='')
        },
        content = function(file){
            genes.vec <- extractGenes( input$genes )
            if(length(genes.vec)==0) return()
            hm=makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, filename=file, main=paste(TITLESTRING, ifelse( input$zscore == 'row', '(Z-scores, ', '('), ifelse( input$allsites == 'all', 'all sites)', 'most variable site)'), sep=' '), height=ifelse(length(genes.vec) < 4, 6.5, NA), zscore=input$zscore, show.sites=input$allsites)
            }
    )
    #############################
    ## download Excel
    output$downloadTab <- downloadHandler(
        filename = function(){paste( FILENAMESTRING, "_", length( extractGenes(input$genes)), "_genes.xlsx", sep='')},
        content = function(file){
            tab=as.data.frame(global$expr.select)
            WriteXLS('tab', ExcelFileName=file, SheetNames=FILENAMESTRING, FreezeCol=6, FreezeRow=3, row.names=T)
            }
    )
})


