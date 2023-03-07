#' helper to get text from python help utility -- may need handling through basilisk
#' @import reticulate
#' @param object a reference to a python module typically with class 'python.builtin.module'
#' @export
py_help2 = function (object) 
{
    help <- reticulate::py_capture_output(reticulate::import_builtins()$help(object), 
        type = "stdout")
    tmp <- tempfile("py_help", fileext = ".txt")
    writeLines(help, con = tmp)
    ans = readLines(tmp)
    unlink(tmp)
    ans
}

#' shiny app that helps access documentation on python-accessible components
#' @import shiny
#' @export
scverse_helper = function() {
  scvi = scviR()
  scanpy = scanpyR()
  
  server = function(input, output) {
   current = reactive({
     avail = list(scvi, scanpy) 
     names(avail) = c("scvi", "scanpy") 
     validate(need(input$curmod %in% names(avail), "pick a module"))
     ans = avail[[input$curmod]]
     if (is.list(ans)) {
       ans = lapply(ans, function(x) {names(x) = x; x})
       }
     ans
     })
     
   getnms = function(targmodule) {
    alln = names(targmodule)
    allc = lapply(alln, function(x) class(targmodule[[x]])[1])
    ismod = which(allc == "python.builtin.module")
    isfun = which(allc %in% c("python.builtin.function"))
    ans = list( modules=alln[ismod], functions=alln[isfun], nonmodfun=alln[-c(ismod,isfun)] )
    ans
    }
   output$topmodules = renderUI({
    nml = getnms( current() )
    radioButtons("topmods", "modules", choices=c(nml$modules))# , nml$functions) )
    })
   output$pickedmodule = renderUI({
    validate(need(length(input$topmods)==1L, "make a selection"))
    att = try(current()[[input$topmods]], silent=TRUE)
    validate(need(class(att)[1] != "python.builtin.function", "function"))
    validate(need(!inherits(att, "try-error") && !is.null(att), "cannot extract input$topmods, try another"))
    radioButtons("lev2", "subtop", choices=names(att), inline=TRUE)
    })
   output$toptext = renderText({
    validate(need(length(input$topmods)==1L, "make a selection"))
    py_help2(current()) # [[input$topmods]])
    }, sep="<br>")
   output$subtext = renderText({
    validate(need(length(input$lev2)==1L, "make a level 2 selection"))
    cur = current()[[input$topmods]]
  
    validate(need(input$lev2 %in% names(cur), sprintf("%s not in names(%s)\n", input$lev2, input$topmods)))
    att = try(current()[[input$topmods]][[input$lev2]], silent=TRUE)
 
    validate(need(!inherits(att, "try-error"), "cannot extract second level, try another"))
    py_help2(att)
    }, sep="<br>")
  
  # deal with stop button
    observe({
              if(input$btnSend > 0)
                 isolate({
                   stopApp(returnValue=0)
                        })
             })
  }
  
  
  
  ui = fluidPage(
   sidebarLayout(
    sidebarPanel(
     radioButtons("curmod", "scverse component", choices=c
      ("scanpy", "scvi"), selected="scvi"),
     uiOutput("topmodules"),
     actionButton("btnSend", "Stop app"),
     width=2
    ),
    mainPanel(
     tabsetPanel(
      tabPanel("elements", uiOutput("pickedmodule"), htmlOutput("subtext")),
      tabPanel("aboutComponent", htmlOutput("toptext")),
      tabPanel("aboutApp", helpText("Provided with scviR to help navigate python documentataion
relevant to scviR.  scanpyR() provides access to the scanpy modules.  anndataR() can also
be used but documentation is not available in this app at present."))
      )
     ) 
    )
   )
  
  runApp(list(ui=ui, server=server))
}
