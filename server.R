library(matrixStats)
source("LimmaRankProd.R")
dat <- NULL
shinyServer(function(input, output,clientData,session) {
  heightSize = function() {
    (input$NumCond-1)*400
  }
  output$messages <- renderText("")
  output$thead <- renderText("Determine differentially regulated features")
  output$description <- renderText("Detection of differentially regulated features becomes tricky 
                                   for data with few replicates and high amounts of missing 
                                   values. By combining the commonly used moderated t-test (limma) with
                                   rank products based on well-defined null hypotheses, these
                                   features can still be detected with high confidence.")
  output$description2 <- renderText("We recommend this method for datasets with a minimum of 1000 features (rows of the table)
                                   and at least 3 replicates. The method is based on ratios between
features within the replicate, i.e. the tests are carried out on paired tests.")
  output$description3 <- renderText("For details, see Schwämmle, V.; Leon, I. R. and Jensen, O. N. 
                                  Assessment and improvement of statistical tools for comparative 
                  proteomics analysis of sparse data sets with few experimental replicates J Proteome Res, 2013, 12, 3874-3883.")
  output$description4 <- renderText("Data input: Data table (csv file) with optionally row and column names. The values should be 
                                    log-transformed intensity/abundance values (not ratios) which have been normalized to be comparable. The 
                                   order of the columns is required to be A1, A2, A3, ..., B1, B2, B3, ..., where
                                   1,2,3 ... are the conditions and A,B,... denote replicates. The tests check for differentially regulated features
                                    versus the \"reference\" condition. For each comparison (log-ratio), we plot the histograms of the uncorrected p-values
                                    and volcano plots of the false discovery rates 
                                    (corrected for multiple testing according to Storey JD. All points above the horizontal line have a q-value below 0.05. A direct approach to false discovery rates. Journal of the Royal Statistical Society. 2002;64:479–498).")
  #   output$NumReps <- renderText("Number of replicates")
  #   output$NumCond <- renderText("Number of conditions")
  #   output$refCond <- renderText("Reference condition")
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      validate(
        need(F,"No data")
      )      
      paste("Results", Sys.Date(), ".csv", sep="");
    },
    content = function(file) {
      write.csv(FullReg, file)
    }
    
  )
  
  output$plot1 <- renderPlot({
    input$button
    isolate({
      dev.control(displaylist="enable")
      updateNumericInput(session,"refCond",max=input$NumCond,value=input$refCond)
      output$messages <- renderText("")
      par(mar=c(3,3,3,3))
      dat <- input$in_file
      NumReps <- input$NumReps
      NumCond <- input$NumCond
      refCond <- input$refCond
      qlim <- input$qval
      if (!is.null(dat)) {
        delim <- input$delimiter
        if (delim == "tab")
          delim <- "\t"
        if (input$row.names){
          dat <- read.csv(input$in_file$datapath,row.names=1,header=input$is_header,sep=delim,dec=input$digits)
          # delete row with empty name
          dat <- dat[!rownames(dat)=="",]
        } else {
          dat <- read.csv(input$in_file$datapath,header=input$is_header,sep=delim,dec=input$digits)
        }
        validate(
          need(mode(as.matrix(dat))=="numeric","Wrong file format / setup")
        )
        print(mode(as.matrix(dat)))
        validate(
          need(ncol(dat) == NumReps*NumCond, "Column number doesn't fit with number of replicates and conditions!")
        )
        if (input$button == 0)
          return()
        if (ncol(dat) == NumReps*NumCond) {
          withProgress(message="Calculating ...", min=0,max=1, {
            
            #           output$messages <- renderText("running")
            Against<-seq(refCond,NumCond*NumReps,NumCond)
            RR<-1:(NumReps*NumCond)
            RR<-RR[-Against]
            RR<-rbind(RR,rep(Against,each=NumCond-1))
            MAData<-NULL
            for (i in 1:ncol(RR)) {
              MAData<-cbind(MAData,dat[,RR[1,i]]-dat[,RR[2,i]])
            }
            rownames(MAData)<-rownames(dat)
            MAReps<-rep(1:(NumCond-1),NumReps)
            incProgress(0.1, detail = paste("Running statistical tests"))
            
            qvalues <- LimmaRankProd(MAData, MAReps)
            incProgress(0.7, detail = paste("Preparing data"))
            
            AllHiLo<-NULL
            LogRatios<-matrix(NA,nrow(MAData),NumCond-1,dimnames=list(rows=rownames(MAData),cols=(1:(NumCond-1))))
            WhereReg<-matrix(F,nrow(MAData),NumCond-1,dimnames=list(rownames(MAData),paste("T",RR[1,1:(NumCond-1)]-1,"vsT",RR[2,1:(NumCond-1)]-1,sep="")))
            Pvalue<-matrix(NA,nrow(MAData),NumCond-1,dimnames=list(rownames(MAData),paste("T",RR[1,1:(NumCond-1)]-1,"vsT",RR[2,1:(NumCond-1)]-1,sep="")))
            for (c in 1:(NumCond-1)) {
              HiLoList<-(qvalues[[c]])[!is.na(qvalues[[c]][,2]) & !is.na(qvalues[[c]][,3]) & (qvalues[[c]][,2]<qlim | qvalues[[c]][,3] < qlim),]
              HiLoList<-HiLoList[order(HiLoList[,4]),]
              LogRatios[rownames(qvalues[[c]]),c]<-qvalues[[c]][,4]
              Pvalue[rownames(qvalues[[c]]),c]<-rowMins(cbind(qvalues[[c]][,2],qvalues[[c]][,3]),na.rm=T)
              colnames(HiLoList)<-c("t-test","limma","rank products","log change")
              AllHiLo<-c(AllHiLo,rownames(HiLoList))
              WhereReg[rownames(HiLoList),c]<-T
            }
            FullReg<-merge(LogRatios,Pvalue,by.x=0,by.y=0,all=T)
            FullReg <- merge(FullReg,WhereReg,by.x=1,by.y=0,all=T)
            rownames(FullReg)<-FullReg[,1]
            FullReg<-FullReg[,2:ncol(FullReg)]
            colnames(FullReg)<-c(paste(colnames(WhereReg),"mean log-ratio"),paste(colnames(WhereReg),"q-value (min(limma,rank products)"),paste(colnames(WhereReg),"regulated"))
            pl <- recordPlot()
            output$downloadData <- downloadHandler(
              filename = function() {
                paste("Results", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(FullReg, file)
              })
            output$downloadFigure <- downloadHandler(
              filename = function() {
                paste("Results", Sys.Date(), ".pdf", sep="");
              },
              content = function(file) {
                pdf(file,height=(NumCond-1)*4)
                print(dev.cur())
                replayPlot(pl)
                dev.off()
              })                
            
          })
          
        }
      }
    })
    #       output$messages <- renderText("finished")
  },height=heightSize)
})
