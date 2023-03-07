  
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
  
