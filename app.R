library(shiny)
library(sqldf)
library(readr)
library(DT)

ui <- fluidPage(theme = "main.css",
  navbarPage("EpGClock",
    tabPanel("Home",fluidRow(column(6,style="background-image:url(\"main.jpeg\");background-repeat:no-repeat;",div(style="height:1500px;")),
                             column(6,h1("Epigenetic Clock"),br(),br(),
      p("Epigenetic clock, also called as DNAmAge is used to predict the DNA methylation age of tissues and cell types using data generated from Illumina Infinium platform.
        This application is developed based on the code and instructions provided by Steve Horvath who developed the methylation age calculator."),
      p(style="font-weight:bold","Citation:",a("DNA methylation age of human tissues and cell types.",href="https://www.ncbi.nlm.nih.gov/pubmed/24138928")),
      style="background-color:white;", div(style = "height:870px;"))
                             )
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
              hr(),
              p("Based on",
              a(href="https://doi.org/10.1186/gb-2013-14-10-r115", "\"DNA methylation age of human tissues and cell types\""),
              "by",
              a(href="https://www.biostat.ucla.edu/people/horvath", "Steve Horvath.")
            )
      ),
    
            mainPanel(DT::dataTableOutput("table")))),

    tabPanel("Correlation with Chronological Age",
             sidebarLayout(
               sidebarPanel(
                 helpText("Upload a .csv file of SampleAnnotation that contains the chronological ages for which the DNAm Age was "),
                 fileInput(
                   inputId = "sampleannoFile",
                   label = "Upload sample annotation .csv file",
                   accept = c('text/csv',
                              'text/comma-separated-values,text/plain',
                              '.csv'),
                   multiple = FALSE),
                 hr(),
                 hr(),
                 helpText("Upload the output file generated using DNAmAge calculation."),
                 fileInput(
                   inputId = "dnamageOutput",
                   label = "Upload DNAmAge output .csv file",
                   accept = c('text/csv',
                              'text/comma-separated-values,text/plain',
                              '.csv'),
                   multiple = FALSE),
                 hr(),
                 hr(),
                 p("Based on",
                   a(href="https://doi.org/10.1186/gb-2013-14-10-r115", "\"DNA methylation age of human tissues and cell types\""),
                   "by",
                   a(href="https://www.biostat.ucla.edu/people/horvath", "Steve Horvath.")
                 )
               ),
               mainPanel(plotOutput("plot"))
             ))
  )
)
server <- function(input,output){
  options(shiny.maxRequestSize = 500 * 1024 ^ 2)
  source("horvath.R")
  
  output$table <- renderDataTable({
    req(input$betaFile)
    withProgress(message = 'Calculation in progress', detail = 'This may take a while...',value = 0, {
      for(i in 1:15){
        incProgress(1/20)
        Sys.sleep(0.25)
      }
    dat0 = read.csv.sql(input$betaFile$datapath,header = TRUE,sep = ",")
    datMiniAnnotation=read.csv("datMiniAnnotation.csv",header = TRUE,sep = ",")
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
  })
  
  output$plot <- renderPlot({
    req(input$sampleannoFile)
    req(input$dnamageOutput)
    datSample = read.csv(input$sampleannoFile$datapath,header = TRUE,sep = ",")
    datDNAmAge = read.csv(input$dnamageOutput$datapath,header = TRUE,sep = ",")
    DNAmAge = datDNAmAge$DNAmAge
    medianAbsDev=function(x,y) median(abs(x-y),na.rm=TRUE)
    medianAbsDev1=signif(medianAbsDev(DNAmAge, datSample$Age),2)
    par(mfrow=c(1,1))
    verboseScatterplot(DNAmAge, datSample$Age,xlab="DNAm Age", ylab="Chronological Age",main=paste("All, err=",medianAbsDev1) );abline(0,1)
  },height = 700, width = 700)
}
shinyApp(ui = ui, server = server)