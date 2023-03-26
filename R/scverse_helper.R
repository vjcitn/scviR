#' helper to get text from python help utility -- may need handling through basilisk
#' @import reticulate
#' @param object a reference to a python module typically with class 'python.builtin.module'
#' @return character vector of lines from python help result
#' @export
pyHelp2 <- function(object) {
  help <- reticulate::py_capture_output(reticulate::import_builtins()$help(object),
    type = "stdout"
  )
  tmp <- tempfile("py_help", fileext = ".txt")
  writeLines(help, con = tmp)
  ans <- readLines(tmp)
  unlink(tmp)
  ans
}

scverseHelperBAD <- function() {
  scvi <- scviR()
  scanpy <- scanpyR()

  server <- function(input, output) {
    current <- reactive({
      avail <- list(scvi, scanpy)
      names(avail) <- c("scvi", "scanpy")
      validate(need(input$curmod %in% names(avail), "pick a module"))
      ans <- avail[[input$curmod]]
      if (is.list(ans)) {
        ans <- lapply(ans, function(x) {
          names(x) <- x
          x
        })
      }
      ans
    })

    getnms <- function(targmodule) {
      alln <- names(targmodule)
      allc <- lapply(alln, function(x) class(targmodule[[x]])[1])
      ismod <- which(allc == "python.builtin.module")
      isfun <- which(allc %in% c("python.builtin.function"))
      ans <- list(modules = alln[ismod], functions = alln[isfun], nonmodfun = alln[-c(ismod, isfun)])
      ans
    }
    output$topmodules <- renderUI({
      nml <- getnms(current())
      radioButtons("topmods", "modules", choices = c(nml$modules)) # , nml$functions) )
    })
    output$pickedmodule <- renderUI({
      validate(need(length(input$topmods) == 1L, "make a selection"))
      att <- try(current()[[input$topmods]], silent = TRUE)
      validate(need(!inherits(att[1], "python.builtin.function"), "can't process function"))
      validate(need(!inherits(att, "try-error") && !is.null(att), "cannot extract input$topmods, try another"))
      radioButtons("lev2", "subtop", choices = names(att), inline = TRUE)
    })
    output$toptext <- renderText(
      {
        validate(need(length(input$topmods) == 1L, "make a selection"))
        pyHelp2(current()) # [[input$topmods]])
      },
      sep = "<br>"
    )
    output$subtext <- renderText(
      {
        validate(need(length(input$lev2) == 1L, "make a level 2 selection"))
        cur <- current()[[input$topmods]]

        validate(need(input$lev2 %in% names(cur), sprintf("%s not in names(%s)\n", input$lev2, input$topmods)))
        att <- try(current()[[input$topmods]][[input$lev2]], silent = TRUE)

        validate(need(!inherits(att, "try-error"), "cannot extract second level, try another"))
        pyHelp2(att)
      },
      sep = "<br>"
    )

    # deal with stop button
    observe({
      if (input$btnSend > 0) {
        isolate({
          stopApp(returnValue = 0)
        })
      }
    })
  }



  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        radioButtons("curmod", "scverse component", choices = c
        ("scanpy", "scvi"), selected = "scvi"),
        uiOutput("topmodules"),
        actionButton("btnSend", "Stop app"),
        width = 2
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

  shinyApp(ui = ui, server = server)
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
}

shinyApp(ui=ui, server=server)
}
