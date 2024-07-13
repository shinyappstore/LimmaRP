shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("in_file", "Input file:",accept=c("txt/csv", "text/comma-separated-values,text/plain", ".csv")),
      fluidRow(
        checkboxInput(inputId="is_header", label="Column names?", value=TRUE),
        checkboxInput(inputId="row.names", label="Row names?", value=TRUE),
        radioButtons(inputId="delimiter", "Delimiter character", choices=c(",",";","tab"), selected=",",inline=T),
        radioButtons(inputId="digits", "Decimal character", choices=c(".",","), selected=".",inline=T)
      ),
      sliderInput("NumReps",min=2,max=20,value=2,label="Number of replicates",step=1),
      sliderInput("NumCond",min=2,max=20,value=3,label="Number of conditions",step=1),
      sliderInput("refCond",min=1,max=20,value=3,label="Reference condition",step=1),
      numericInput("qval","q-value threshold",value=0.01,min=0,max=1),
      br(),
      actionButton("button","Run analysis"),
      br(),
      downloadLink('downloadData', 'Download results'),br(),
      downloadLink('downloadFigure', 'Download figure'),
      br(),
      br(),
      textOutput("messages")
    ),
    mainPanel(
      h3("Method to determine differentially regulated features"),
      textOutput("description"),
      textOutput("description2"),
      textOutput("description3"),
      textOutput("description4"),
      br(),
      plotOutput("plot1",height="auto")
    ))
))
