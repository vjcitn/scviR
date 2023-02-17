#' get lmFit for heterogeneity across subclusters
#' @param inlist list of SingleCellExperiments (SCEs) formed by scran::quickSubCluster
#' @param clname character(1) name of cluster SCE to assess
#' @note It is assumed that 'logcounts' is an assay element,
#' and that 'subcluster' is a colData element of each SCE in inlist
#' @examples
#' data(all.sce)
#' lm3 = get_subcl_LM(all.sce, "3")
#' names(lm3)
#' @export
get_subcl_LM = function(inlist, clname) {
  se = inlist[[clname]]
  x = assay(se, "logcounts")
  mm = model.matrix(~subcluster, data=as.data.frame(colData(se)))
  lmFit(x, mm)
}

#' get lmFit F-stat based collection of n genes most varying in mean across subclusters
#' @param inlist list of SingleCellExperiments (SCEs) formed by scran::quickSubCluster
#' @param clname character(1) name of cluster SCE to assess
#' @param n numeric(1) number to preserve
#' @return list with two elements, feat = rowData corresponding to variable genes, stats = topTable result
#' @note Symbol will be taken from feat and placed in stats component if available
#' @examples
#' data(all.sce)
#' scl = get_subclustering_features(all.sce, "3", 10)
#' names(scl)
#' @export
get_subclustering_features = function(inlist, clname, n=20) {
  lm1 = get_subcl_LM( inlist, clname )
  p = seq(2, ncol(lm1$coef)) # to get F stats
  suppressWarnings({
    lm1 = eBayes(lm1) # lots of zeroes
  })
  tt = topTable(lm1, p, n=n)
  en = rownames(tt)
  rd = rowData(inlist[[1]][en,])
  tt$gene = en
  if ("Symbol" %in% colnames(rd)) tt = data.frame(gene=rd$Symbol, tt)
  list(feat=rowData(inlist[[clname]][rownames(tt),]), stats=tt)
}

#' app to explore diversity in RNA-subclusters within ADT clusters
#' @param sce a SingleCellExperiment with altExp with ADT quantification
#' @param inlist list of SingleCellExperiments (SCEs) formed by scran::quickSubCluster
#' @param adtcls vector of ADT cluster assignments
#' @note TSNE should already be available in `altExp(sce)`; follow OSCA book 12.5.2.  If using
#' example, set `ask=FALSE`.
#' @examples
#' if (interactive()) {
#'  data(sce)
#'  data(all.sce)
#'  data(clusters.adt)
#'  explore_subcl( sce, all.sce, clusters.adt )
#' }
#' @export
explore_subcl = function( sce, inlist, adtcls ) {
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("Explore CITE-seq subclusters"),
    selectInput("clpick", "ADT cluster", choices=names(inlist), selected="3"),
    selectInput("baseADT", "Base for smooths", choices=rowData(altExp(sce))$ID, selected="CD127",
       multiple=FALSE),
    uiOutput("feats"), width=2
    ),
   mainPanel(
    tabsetPanel(
     tabPanel("tsne", helpText("Guide to ADT-based clusters"), plotOutput("tsne")),
     tabPanel("heatmap", helpText("Guide to protein abundance profiles"), plotOutput("heatmap")),
     tabPanel("boxplots", plotOutput("boxplots")),
     tabPanel("smooths", plotOutput("smooths")),
     tabPanel("stats", dataTableOutput("stats")),
     tabPanel("about", helpText("This app helps to explore RNA-based subclusters of ADT-based clusters
 formed according to ch 12.6.1 of the OSCA book.  Inputs are a basic SingleCellExperiment with
logcounts for RNA and ADT features, a list of SCE formed using scran::quickSubCluster, and
the vector of assignments from cells to ADT subclusters."),
   helpText(" "),
   helpText("The TSNE map of ADT clusters is
given, along with a heatmap for ADT abundances, as guides."),  
   helpText(" "),
   helpText("Choose an ADT-based cluster using the labeling on the
TSNE map and F tests will be performed (using limma) to identify genes whose mean abundances vary strongly
across RNA-based subclusters.  Boxplot tab presents marginal distributions of expression
of selected genes in RNA-based subclusters.  Smooths tab depicts association between RNA abundance and protein
abundance for selected genes and a given protein in the ADT panel."))
     )
    )
   )
  )
 server = function(input, output, session) {
  output$tsne = renderPlot({
   plotTSNE(altExp(sce), colour_by="label", text_by = "label", text_color="red")
   })
  output$heatmap = renderPlot({
   se.avg = sumCountsAcrossCells(altExp(sce), adtcls, exprs_values = "logcounts", average=TRUE)
   avg = assay(se.avg)
   pheatmap::pheatmap(avg - rowMeans(avg), breaks=seq(-3, 3, length.out=101))
   })
  featdata = reactive({
     get_subclustering_features(inlist, input$clpick, n=10) 
     })
  output$feats = renderUI({
    scl = featdata()$feat
    checkboxGroupInput("genes", "genes for boxplots and smooths", choices=scl$Symbol, selected=scl$Symbol[1:3])
    })
  output$stats = renderDataTable({
    tab = featdata()$stats
    cl = which(sapply(tab, is.numeric))
    for (j in cl) tab[[j]] = round(tab[[j]], 4)
    tab
    })
  output$boxplots = renderPlot({
    featdata()$feat
    plotExpression(inlist[[input$clpick]], x="subcluster", features=input$genes,
         swap_rownames = "Symbol", ncol=length(input$genes))
    })
  output$smooths = renderPlot({
    featdata()$feat
    plotExpression(inlist[[input$clpick]], x=input$baseADT, features=input$genes,
         show_smooth=TRUE, show_se=FALSE,
         swap_rownames = "Symbol", ncol=length(input$genes))
    })
 }
 runApp(list(ui=ui, server=server))
}
   


#zz = lmFit(assay(all.sce[["3"]], "logcounts"), model.matrix(~subcluster, data=as.data.frame(colData(all.sce[["3"]]))#, data=all.sce[[3]])
#)
#library(citeseqApp)
#data(all.sce)
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1"), swap_rownames="Symbol", ncol=3)
#options(bitmapType="cairo")
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1"), swap_rownames="Symbol", ncol=3)
#dev.off()
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1"), swap_rownames="Symbol", ncol=3)
#all.sce[[3]]
#zz = lmFit(logcounts~subcluster, data=all.sce[[3]])
#?lmFit
#zz = lmFit(assay(all.sce[["3"]], "logcounts"), model.matrix(~subcluster, data=as.data.frame(colData(all.sce[["3"]]))#, data=all.sce[[3]])
#)
#)
#ezz = eBayes(zz)
#ezz[[1]][1,]
#topTable(ezz, 2:3)
#rowData(all.sce[[3]])[rownames(.Last.value),]
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1", "LTB", "NKG7"), swap_rownames="Symbol", ncol=5)
