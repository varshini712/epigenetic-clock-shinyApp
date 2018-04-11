library(shiny)
library(sqldf)
library(readr)


ui <- fluidPage(theme = "main.css",
  navbarPage("EpGClock",
    tabPanel("Home",fluidRow(column(5,HTML('<h1>Epigenetic Clock</h1>')),
                          column(9,HTML('<img src = "main.jpeg" class=responsive></img>')))
    ),
    tabPanel("Calculate DNAm Age",
            sidebarLayout(
              sidebarPanel(
                  helpText("Upload a .csv file of beta value (beta matrix). The first column has the CpG probe names, the following column(s) has methylation values"),
                  fileInput(
                  inputId = "betaFile",
                  label = "Upload beta matrix .csv file (max 500MB)",
                  accept = c('text/csv',
                            'text/comma-separated-values,text/plain',
                            '.csv'),
                  multiple = FALSE
                          ),
              hr(),
              helpText("Enter a sample id (for example hESCs_1) and press the button"),
              textInput(inputId = "sampleID", label = "Enter a sampleID/sampleName", placeholder = ""),
              actionButton(inputId = "calculate", label = "Predict age"),
              hr(),
              p("Based on",
              a(href="https://doi.org/10.1186/gb-2013-14-10-r115", "\"DNA methylation age of human tissues and cell types\""),
              "by",
              a(href="https://www.biostat.ucla.edu/people/horvath", "Steve Horvath.")
            )
      ),
    
            mainPanel(tableOutput("table")))),

    tabPanel("Correlation with Chronological Age")
  )
)
server <- function(input,output){
  options(shiny.maxRequestSize = 500 * 1024 ^ 2)
  source("horvath.R")
  output$table <- renderTable({
    req(input$betaFile)
    dat0 = read.csv.sql(input$betaFile$datapath,header = TRUE,sep = ",")
    datMiniAnnotation=read.csv("datMiniAnnotation.csv")
    match1=match(datMiniAnnotation[,1], dat0[,1] )
    dat0Reduced=dat0[match1,]
    dat0Reduced[,1]=as.character(dat0Reduced[,1])
    dat0Reduced[is.na(match1),1]= as.character(datMiniAnnotation[is.na(match1),1])
    datout=data.frame(dat0Reduced)
    # make sure you output numeric variables...
    for (i in 2:dim(datout)[[2]] ){datout[,i]=
      as.numeric(as.character(gsub(x=datout[,i],pattern="\"",replacement=""))) }
    #replace "MethylationData" with a filename of your choice
    colnames(datout)[1] <- "ProbeID"
    
    ##Calculate DNAm age tabsetpanel##
    final <- calculateAge(datout)
  })
}
shinyApp(ui = ui, server = server)

#fluidRow(HTML('<img src = "main1.png" class="responsive"></img>')),
#fluidRow(HTML('<p>DNAm Age Calculator</p>'))