#' helper to get text from python help utility -- may need handling through basilisk
#' @import reticulate
#' @param object a reference to a python module typically with class 'python.builtin.module'
#' @return character vector of lines from python help result
#' @export
pyHelp2 <- function(object) {
  help <- reticulate::py_capture_output(reticulate::import_builtins()$help(object))
  tmp <- tempfile("py_help", fileext = ".txt")
  writeLines(help, con = tmp)
  ans <- readLines(tmp)
  ans = gsub("<", "&lt;", ans)
  ans = gsub(">", "&gt;", ans)
  unlink(tmp)
  ans
}






#' shiny app that helps access documentation on python-accessible components
#' @import shiny
#' @return shinyApp instance
#' @export
scviHelper = function() {
# build 2-level hierarchy in advance; it cannot be done in the reactive
# framework, it seems
scvi = scviR()
oktypes = c("python.builtin.module", "python.builtin.function", "python.builtin.type")

nscvi = names(scvi)
sidebarInds =  which(sapply(nscvi, function(x) class(scvi[[x]])[1] %in% oktypes))
sidebaropts = nscvi[sidebarInds]

nm2 = lapply(sidebaropts, function(x)names(scvi[[x]]))
names(nm2) = sidebaropts
lens = sapply(nm2,length)
if (any(lens==0)) nm2 = nm2[-which(lens==0)]
nmcl = lapply(names(nm2), function(x) {
   kp = which(sapply(names(scvi[[x]]), function(z) class(scvi[[x]][[z]])[1] %in% oktypes))
   names(scvi[[x]])[kp]
})
names(nmcl) = names(nm2)
lens = sapply(nmcl,length)
if (any(lens==0)) nmcl = nmcl[-which(lens==0)]
    
pylist1 = lapply( names(nmcl), function(x)
   scvi[[x]] )
names(pylist1) = names(nmcl)
pylist2 = lapply( seq_len(length(pylist1)), function(x) {
  ans = lapply( nmcl[[x]], function(z) pylist1[[x]][[z]] )
  names(ans) = nmcl[[x]] 
  ans
  })
names(pylist2) = names(pylist1)
   

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("Get an overview of resources available from scvi-tools"),
   radioButtons("toplev", "components", choices=names(nmcl)),
   actionButton("btnSend", "Stop app"),
   width=2),
  mainPanel(
   tabsetPanel(
    tabPanel("Doc",
     uiOutput("seclev"),
     htmlOutput("helptxt")
     ),
    tabPanel("About",
     helpText("This app helps users to survey tha resources available in scvi-tools.
Content is obtained using a variant of reticulate::py_help applied to module
documentation installed in the basilisk cache for scviR."))
    )
   )
  )
)

server = function(input, output) {
  output$seclev = renderUI({
    radioButtons("sec", "components", choices=nmcl[[input$toplev]], inline=TRUE)
    })
  output$helptxt = renderText({
    validate(need(nchar(input$sec)>0,"waiting"))
    m2 = pylist2[[ input$toplev ]][[input$sec]]
    pyHelp2( m2 )
    }, sep="<br>")
    # deal with stop button
  observe({
      if (input$btnSend > 0) {
        isolate({
          stopApp(returnValue = 0)
        })
      }
    })
}

shinyApp(ui=ui, server=server)
}

#' shiny app that helps access documentation on python-accessible components
#' @import shiny
#' @return shinyApp instance
#' @export
scanpyHelper = function() {
# build 2-level hierarchy in advance; it cannot be done in the reactive
# framework, it seems
scanpy = scanpyR()
oktypes = c("python.builtin.module", "python.builtin.function", "python.builtin.type")

nscanpy = names(scanpy)
sidebarInds =  which(sapply(nscanpy, function(x) class(scanpy[[x]])[1] %in% oktypes))
sidebaropts = nscanpy[sidebarInds]

nm2 = lapply(sidebaropts, function(x)names(scanpy[[x]]))
names(nm2) = sidebaropts
lens = sapply(nm2,length)
if (any(lens==0)) nm2 = nm2[-which(lens==0)]
nmcl = lapply(names(nm2), function(x) {
   kp = which(sapply(names(scanpy[[x]]), function(z) class(scanpy[[x]][[z]])[1] %in% oktypes))
   names(scanpy[[x]])[kp]
})
names(nmcl) = names(nm2)
lens = sapply(nmcl,length)
if (any(lens==0)) nmcl = nmcl[-which(lens==0)]
    
pylist1 = lapply( names(nmcl), function(x)
   scanpy[[x]] )
names(pylist1) = names(nmcl)
pylist2 = lapply( seq_len(length(pylist1)), function(x) {
  ans = lapply( nmcl[[x]], function(z) pylist1[[x]][[z]] )
  names(ans) = nmcl[[x]] 
  ans
  })
names(pylist2) = names(pylist1)
   

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("Get an overview of resources available from scanpy"),
   radioButtons("toplev", "components", choices=names(nmcl)),
   actionButton("btnSend", "Stop app"),
   width=2),
  mainPanel(
   tabsetPanel(
    tabPanel("Doc",
     uiOutput("seclev"),
     htmlOutput("helptxt")
     ),
    tabPanel("About",
     helpText("This app helps users to survey tha resources available in scanpy.
Content is obtained using a variant of reticulate::py_help applied to module
documentation installed in the basilisk cache for scviR."))
    )
   )
  )
)

server = function(input, output) {
  output$seclev = renderUI({
    radioButtons("sec", "components", choices=nmcl[[input$toplev]], inline=TRUE)
    })
  output$helptxt = renderText({
    validate(need(nchar(input$sec)>0,"waiting"))
    m2 = pylist2[[ input$toplev ]][[input$sec]]
    pyHelp2( m2 )
    }, sep="<br>")
    # deal with stop button
  observe({
      if (input$btnSend > 0) {
        isolate({
          stopApp(returnValue = 0)
        })
      }
    })
}

shinyApp(ui=ui, server=server)
}
