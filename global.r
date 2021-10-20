#################################################################
## Filename: global.r
## Created: March 19, 2016
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data of the BKM120 project
## This file imports the underlying data and contains global functions
## and variables used by 'ui.R' and 'server.R'
#################################################################

source('pheatmap.r')
library(scales)
library(gtable)

## import the data
load('BKM_data_2016-07-13.RData')

## global parameters
GENESSTART <<- 'PIK3R2 PIK3R1 AKT1 AKT2 AKT3 MTOR RPS6KB1 RPS6KA1 TSC2 PTEN PDK1 PDPK1 RB1'
GENEMAX <<- 20
TITLESTRING <<- '<font size="5" face="times"><i><b>"Mass spectrometry-based proteomics reveals potential roles of NEK9 and MAP2K4 in resistance to PI3K inhibitors in triple negative breast cancers"</b></i> (<a href="http://cancerres.aacrjournals.org/content/78/10/2732" target="_blank_">Mundt <i>et al.</i> Cancer Research. 2018</a>)</font><br>'
FILENAMESTRING <<- 'BKM120'
GAPSIZEROW <<- 20


##library(pheatmap)
library(RColorBrewer)
library(gplots)
library(WriteXLS)
library(grid)
library(bcrypt)

##################################################################
## 21060613 bcrypt
authenticateUser <- function(passphrase){
    if(nchar(as.character(passphrase)) > 0){
        return(checkpw(as.character(passphrase), "$2a$12$XXsIDMNZMaOgRcaq8JX8muN3gsG93YGltGcqfiqkGGgDkJKV4cwri"))
    } else {
        return(FALSE)
    }
}

##################################################################
## function to extract gene names from a string of
## character
extractGenes <- function(genes.char){

    gene.max=GENEMAX

    ##cat('list: ', genes.char, '\n')
    ##cat('end\n')

    ##if(!global$auth)return()
    if(is.null(genes.char)) return(NULL)
    if( nchar(genes.char) == 0 )return(NULL)


    ## extract genes
    genes.vec= unlist(strsplit(genes.char, ','))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ' '))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ';'))

    ## unique gene names
    genes.vec <- unique(genes.vec)

    ## limit to 'gene.max' genes
    if(length(genes.vec) > gene.max){
        warning(paste('more than', gene.max,'gene ids submitted! Showing results for the first 20 genes in the list.\n'))
        genes.vec <- genes.vec[1:gene.max]
    }
    return(genes.vec)
}


##################################################################
## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n.entries, n.genes){

    height = (n.entries+2)*11 + (n.genes-1)*GAPSIZEROW + 140

    return(height)
}

#######################################################
## find all entries with associated gene name in
## the dataset. returns vector of indices.
findGenesInDataset <- function(gene, show.sites){

    ## remove spaces
    gene <- gsub(' ', '', gene )
    gene <- unique(gene)
    ## remove emtpy strings
    gene.nchar=which(nchar(gene) == 0)
    if(length(gene.nchar) > 0)
        gene <- gene[-gene.nchar]

    if(length(gene) == 0) return()

    ## check whether the genes are present in the dataset
    gene.idx <- grep( paste(paste('(^|,)', gene, '($|,)', sep=''), collapse='|'), gsub(' ', '', row.anno[, 'Gene.id']) )
    if( length(gene.idx) == 0 ){
        stop('None of the gene ids you have entered could be found in the dataset!\n')
    }
    ## use row names
    gene.idx <- rownames(row.anno)[gene.idx]

    ## exract data and remove empty rows
    data.tmp <- tab.expr.all[gene.idx, ]

    ## remove empty rows
    rm.idx <- apply( data.tmp, 1, function(x) ifelse( sum(is.na(x) )/length(x) == 1, 1, 0 ))
    if(sum(rm.idx) > 0)
        gene.idx <- gene.idx[-which(rm.idx == 1)]

    ## update row annotation
    row.anno.tmp <- row.anno[gene.idx, ]

    ## update data
    data.tmp <- data.tmp[gene.idx, ]

    ## most variable site
    if(show.sites=='most variable'){

        ## extract MS phospho
        ms.pSTY.idx <- grep('4_MS_pSTY', row.anno.tmp$Data.type)
        if( length(ms.pSTY.idx) > 0 ){

            ms.pSTY.sd <- apply(data.tmp[ ms.pSTY.idx, ], 1, sd, na.rm=T)

            rm.idx <- tapply(ms.pSTY.sd, row.anno.tmp[ms.pSTY.idx, 'Gene.id'],  function(x) names(x)[which.max(x)])
            rm.idx <- setdiff( names(ms.pSTY.sd), unlist(rm.idx) )

            gene.idx <- setdiff(gene.idx, rm.idx)

            data.tmp <- data.tmp[gene.idx, ]
            row.anno.tmp <- row.anno.tmp[gene.idx, ]

        }
        ## extract RPPA phospho
        rppa.pSTY.idx <- grep('5_RPPA_pSTY', row.anno.tmp$Data.type)

        if( length(rppa.pSTY.idx) > 0 ){

            rppa.pSTY.sd <- apply(data.tmp[ rppa.pSTY.idx, ], 1, sd, na.rm=T)

            rm.idx <- tapply(rppa.pSTY.sd, row.anno.tmp[rppa.pSTY.idx, 'Gene.id'],  function(x) names(x)[which.max(x)])
            rm.idx <- setdiff( names(rppa.pSTY.sd), unlist(rm.idx) )

            gene.idx <- setdiff(gene.idx, rm.idx)

            data.tmp <- data.tmp[gene.idx, ]
            row.anno.tmp <- row.anno.tmp[gene.idx, ]

        }
    }
    return(gene.idx)
}


#################################################################
## draw the actual heatmap
##
#################################################################
makeHM <- function(gene, filename=NA, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, zscore='row', show.sites='all', min.val=-3, max.val=3, ...){

    min.val <- as.numeric(min.val)
    max.val <- as.numeric(max.val)

    n.bins=12

    ## remove any blanks in the gene names
    row.anno$Gene.id <- gsub(' ', '', row.anno$Gene.id )

    #################################
    ## find all entries of submitted gene names
    gene.idx <- findGenesInDataset(gene, show.sites)
    gene.idx=gene.idx[order( row.anno[gene.idx, 'Gene.id'] )]


    #################################
    ## extract genes of interest
    expr.select <- data.matrix(expr[gene.idx, ])
    row.anno.select <- row.anno[gene.idx, ]

    #################################
    ## reorder the columns
    ord.idx.col <- c(

        grep('W12_.*\\.B\\.', colnames(expr.select)),
        grep('W12_.*\\.V\\.', colnames(expr.select)),
        grep('W12_.*\\.WO\\.', colnames(expr.select)),


        grep('W2_.*\\.B\\.', colnames(expr.select)),
        grep('W2_.*\\.V\\.', colnames(expr.select)),
        grep('W2_.*\\.WO\\.', colnames(expr.select)),

        grep('W6_.*\\.B\\.', colnames(expr.select)),
        grep('W6_.*\\.V\\.', colnames(expr.select)),
        grep('W6_.*\\.WO\\.', colnames(expr.select)),

        grep('W21_.*\\.B\\.', colnames(expr.select)),
        grep('W21_.*\\.V\\.', colnames(expr.select)),
        grep('W21_.*\\.WO\\.', colnames(expr.select)),

        grep('W30_.*\\.B\\.', colnames(expr.select)),
        grep('W30_.*\\.V\\.', colnames(expr.select)),
        grep('W30_.*\\.WO\\.', colnames(expr.select)),

        grep('W4_.*\\.B\\.', colnames(expr.select)),
        grep('W4_.*\\.V\\.', colnames(expr.select)),
        grep('W4_.*\\.WO\\.', colnames(expr.select))
    )

    expr.select <- expr.select[ , ord.idx.col]
    column.anno <- column.anno[colnames(expr.select),]


    #################################
    ## extract sample ids
    sampleIDs <- colnames(expr.select)

    #####################################################
    ## labels for the rows in the heatmap
    featureIDs.anno.select <- paste(row.anno.select[ , 'Gene.id'], row.anno.select[ , 'Data.type'])
    ord.idx.row <- order(featureIDs.anno.select)

    ####################################################################
    ## reorder the rows: mRNA, MS prot, RPPA prot, MS pSTY, RPPA pSTY
    featureIDs.anno.select <- featureIDs.anno.select[ord.idx.row]
    row.anno.select <- row.anno.select[ord.idx.row, ]
    expr.select <- expr.select[ord.idx.row, ]

    #####################################################
    ## clean up the names
    featureIDs.anno.select <- sub('[0-9]_','', featureIDs.anno.select)
    featureIDs.anno.select <- gsub('_',' ', featureIDs.anno.select)

    #####################################################
    ## add phosphosite annotation
    #####################################################
    ## RPPA
    rppa.psty.idx <- grep('RPPA pSTY', featureIDs.anno.select)
    if(length(rppa.psty.idx)>0){
        featureIDs.anno.select[ rppa.psty.idx ] <- paste( sub(' pSTY', '', featureIDs.anno.select[ rppa.psty.idx ]),  row.anno.select[ rppa.psty.idx, 'ID'])
    }
    ## MS
    ms.psty.idx <- grep('MS pSTY', featureIDs.anno.select)
    if(length(ms.psty.idx)>0){
        featureIDs.anno.select[ ms.psty.idx ] <- paste( sub(' pSTY', '', featureIDs.anno.select[ ms.psty.idx ]), paste('p', sub('.*_([S|T|Y][0-9]*.*)$', '\\1', row.anno.select[ ms.psty.idx, 'ID']), sep='') )
    }

    ################################
    ## remove any remaining spaces in the site ids
    featureIDs.anno.select <- gsub(' _',' ', featureIDs.anno.select)

    #################################
    ## apply zscore
    rownames(expr.select) <- make.unique(featureIDs.anno.select)

    if(zscore == 'row'){
        expr.select.zscore <-  t(apply(expr.select, 1, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T)))
    } else {
        expr.select.zscore <- expr.select
    }

    ##################################
    ## cap at -3/3
    expr.select.zscore[ which(expr.select.zscore < min.val) ] <- min.val
    expr.select.zscore[ which(expr.select.zscore > max.val) ] <- max.val

    ##############################
    ## column annotation
    column.anno.fig <- data.frame(
        Treatment=as.character( column.anno[,c('Treatment') ]),
        WHIM=as.character(column.anno[,c('WHIM') ])

    )
    rownames(column.anno.fig) <- rownames(column.anno)

    ##############################
    ## colors for column annotation
    column.anno.col <- list(
        Treatment=grey(seq(0,.95,length.out=length(unique(column.anno.fig$Treatment)))),
        WHIM=c(W12='darkred', W2='darkblue', W21='darkgreen', W4='cadetblue2', W6='burlywood3', W30='darkgoldenrod'),
        idx=c(grey='grey', black='black')
    )
    names(column.anno.col[['Treatment']]) <- unique(as.character(column.anno.fig$Treatment))

    #####################################################
    ## row annotation - use alternating colors
    ## to make it easier to distinguish different genes
    row.anno.idx <- vector('character', nrow(row.anno.select))
    for(i in 1:length(unique(row.anno.select$Gene.id))){
        col.tmp=ifelse(i %% 2 == 0, 'black' , 'grey')
        row.anno.idx[which(row.anno.select$Gene.id == unique(row.anno.select$Gene.id)[i])] <- col.tmp
    }
    row.anno.fig <- data.frame( idx=as.factor(row.anno.idx)  )
    ## View(row.anno.fig)
    ##row.anno.col <- list(black='black', grey='grey')
    ################################
    ## gaps
    ################################
    ## only possible because matrix is ordered according to WHIMS
    gaps.column=cumsum(c(  table(column.anno.fig$WHIM) ))
    gaps.row=cumsum(table(sub(' .*', '', featureIDs.anno.select)))

    ################################
    ## colors misc
    color.breaks = seq(min.val, max.val, length.out=n.bins)
    color.hm =  colorRampPalette( c('blue', 'grey', 'red'))(length(color.breaks))
    color.border = 'white'

    ################################
    ## legend breaks
    legend_breaks=seq(min.val, max.val, 1)
    ##legend_labels=c('-3', '-2', '-1', ' 0', '+1', '+2' ,'+3')
    legend_labels=as.character(legend_breaks)

    ###############################
    ## heatmap
    cellwidth=20
    cellheight=10

    pheatmap( expr.select.zscore, cluster_row=F, cluster_col=F,  annotation_col=column.anno.fig, annotation_colors=column.anno.col, scale = "none", labels_row=featureIDs.anno.select, border_color=color.border, gaps_col=gaps.column, gaps_row=gaps.row, color=color.hm, filename=filename, cellwidth=cellwidth, cellheight=cellheight, labels_col=sampleIDs, breaks=color.breaks, legend_breaks=legend_breaks, legend_labels=legend_labels, na_col='white', gapsize_row=GAPSIZEROW, ... )

    ##############################################
    ## assembel table for download
    ## add row annotation
    mat.row <- as.data.frame( cbind( row.anno.select, expr.select, deparse.level=0 ), stringsAsFactors=F )
    rownames(mat.row) <- make.unique(featureIDs.anno.select)

    ## column annotation
    mat.col <- as.data.frame( cbind( matrix('', nrow=ncol(column.anno.fig), ncol=ncol(row.anno.select)), t(column.anno.fig), deparse.level=0), stringsAsFactors=F)
    colnames(mat.col) <- colnames(mat.row)

    ## put everything together
    mat <- rbind( mat.col, mat.row)

    return(mat)
}
