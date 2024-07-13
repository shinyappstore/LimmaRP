# # function to calculate p-values for a data set consisting of at least 3 replicates and missing values
# until now it will only look for relative regulations within ratios
# the algorithm will compare limma, rank products and t-test
# 
# Input values are:
# Data: Data set consisting of log-ratios for abundance values for different replicates and conditions
# Reps: vector with specification of ratios beginning with 1 ; 
# (1,2,3,1,2,3) means replicate 1 ratio 1, replicate 1 ratio 2 ,...
#
# Output values are:
# list for each ratio with 5 columns each, corresponding to feature name, q-values from t-test, 
# q-values from limma, q-values from rank products, mean of log-ratios
# 
# The script also visualizes p-value distributions and vulcano plots for the q-values
# For more information and citations: 
# Veit Schwämmle , Ileana Rodríguez León , and Ole Nørregaard Jensen
# Assessment and Improvement of Statistical Tools for Comparative Proteomics Analysis of Sparse 
# Data Sets with Few Experimental Replicates
# J. Proteome Res., Article ASAP
# DOI: 10.1021/pr400045u

## 29.4.2014: commentted filtering for one replicate as it also deletes all other feature values

LimmaRankProd <- function(Data,Reps,Titles=NA) {
  
#   library(multtest)
  library(limma)
  library(genefilter)
  library(qvalue)
  library(gplots)
  
  ## Generate MA matrix
  MAData<-Data
  MAReps<-Reps
  # deleting entries with only one element
  NumCond<-max(MAReps)
  NumReps<-min(table(MAReps))
#   for (c in 1:NumCond) {
#     MAData<-MAData[rowSums(!is.na(MAData[,MAReps==c]))>1,MAReps==c]
#   }
  
  
  print(paste("NumCond: ",NumCond,"NumReps: ",NumReps))
  
  ##limma with ratios
  design<-NULL
  plvalues<-NULL
  for (c in (1:NumCond)) {
    design<-cbind(design,as.numeric(MAReps==c))
  }
  lm.fittedMA <- lmFit(MAData,design)
  lm.bayesMA<-eBayes(lm.fittedMA)
  topTable(lm.bayesMA)
  if (!is.null(dim(lm.bayesMA$p.value))) {
    for (i in 1:NumCond) 
      plvalues[[i]]<-lm.bayesMA$p.value[,i]
  } else {
    plvalues[[1]]<-lm.bayesMA$p.value
  }
  ###################
  
  ptvalues<-NULL
  pRPvalues<-NULL
  for (vs in 1:NumCond) {
    tMAData<-MAData[,MAReps==vs]
    ptMAvalues<-NULL
    ## MA t-test_pvalues
    for (pep in 1:(dim(tMAData)[1])) {
      if(sum(!is.na(tMAData[pep,]))>1) {
        ptMAvalues<-append(ptMAvalues,t.test(as.vector(tMAData[pep,]))$p.value)
      } else {
        ptMAvalues <- append(ptMAvalues,NA)
      }
    }
    names(ptMAvalues)<-rownames(tMAData)
    ptvalues[[vs]] <- ptMAvalues
    
    ## rank products
    tRPMAData<-MAData[,MAReps==vs]
    NumElements<-rowSums(!is.na(tRPMAData))
    RPMAownUp_pvalues<-RPMAownDown_pvalues<-NULL
    for (d in unique(NumElements)) {
      RPMAData<-tRPMAData[NumElements==d,]
      if(d>1 && length(as.matrix(RPMAData))>ncol(tRPMAData)) {
        RP.own<-0
        Rank<-NULL
        RankNAs<-0
        for (r in 1:NumReps) {
          Rank[[r]]<-rank(RPMAData[,r],na.last="keep")/(sum(!is.na(RPMAData[,r]))+1)
          names(Rank[[r]]) <-rownames(RPMAData)
          Rank[[r]][is.na(Rank[[r]])]<-1
          RP.own<-RP.own+log(Rank[[r]])
          RankNAs<-RankNAs+sum(Rank[[r]]>1)
        }
        RP.own<-exp(RP.own)
        RPownCorr<- -log(RP.own)
        RPMAownUp_pvalues<-c(RPMAownUp_pvalues,pgamma(RPownCorr,d))
        RP.own<-0
        Rank<-NULL
        RankNAs<-0
        for (r in 1:NumReps) {
          Rank[[r]]<-rank(-RPMAData[,r],na.last="keep")/(sum(!is.na(RPMAData[,r]))+1)
          names(Rank[[r]]) <-rownames(RPMAData)
          Rank[[r]][is.na(Rank[[r]])]<-1
          RP.own<-RP.own+log(Rank[[r]])
          RankNAs<-RankNAs+sum(Rank[[r]]>1)
        }
        RP.own<-exp(RP.own)
        RPownCorr<- -log(RP.own)
        RPMAownDown_pvalues<-c(RPMAownDown_pvalues,pgamma(RPownCorr,d))
      }
    }
    RPMAown_pvalues<-2*apply(cbind(RPMAownDown_pvalues,RPMAownUp_pvalues),1,min)
    RPMAown_pvalues<-RPMAown_pvalues[RPMAown_pvalues<1]
    pRPvalues[[vs]]<-RPMAown_pvalues
  }
  
  OutData<-NULL
  
  ## q-values
  qlvalues<-qtvalues<-qRPvalues<-NULL
    par(mfrow=c(2*NumCond,3))
  for (i in 1:NumCond) {
    ttval<-ptvalues[[i]]
    ttval<-ttval[!is.na(ttval)]
    print(paste(length(ttval),"ptvalues"))
    countQts<-qvalue(ttval)
    qtvalues[[i]]<-countQts$qvalues
    
    ttval<-plvalues[[i]]
    ttval<-ttval[!is.na(ttval)]
    print(paste(length(ttval),"plvalues"))
    countQls<-qvalue(ttval)
    qlvalues[[i]]<-countQls$qvalues
    
    ttval<-pRPvalues[[i]]
    ttval<-ttval[!is.na(ttval)]
    print(paste(length(ttval),"pRPvalues"))
    countQRPs<-qvalue(ttval,lambda=seq(0.5,0.05,-0.05))
    qRPvalues[[i]]<-countQRPs$qvalues
    
    hist(ptvalues[[i]],100,main="t-test",xlab="p-values")
    a<-hist(plvalues[[i]],100,main=paste("limma,","Comparison",i),xlab="p-values")
    if (sum(!is.na(Titles))>0 & length(Titles)==NumCond) {
      text(x=0.5,y=0.5*max(a$counts,na.rm=T),Titles[i],cex=1.2)
    }
    hist(pRPvalues[[i]],100,main="rank products",xlab="p-values")
    
    qvulcdat<-merge(qtvalues[[i]],qlvalues[[i]],all=T,by=0)
    qvulcdat<-merge(qvulcdat,qRPvalues[[i]],all=T,by.x=1,by.y=0)
    qvulcdat<-merge(qvulcdat,rowMeans(MAData[,Reps==i],na.rm=T),all=T,by.x=1,by.y=0)
    rownames(qvulcdat)<-qvulcdat[,1]
    qvulcdat<-qvulcdat[,2:ncol(qvulcdat)]
    maxy<-max(-log(qvulcdat[,1:3]),na.rm=T)
    plot(qvulcdat[,ncol(qvulcdat)],-log(qvulcdat[,1]),main="t-test",xlab="log(ratio)",ylab="-log(q-value)",cex=0.2,ylim=c(1,maxy))
    abline(-log(0.05),0)
    plot(qvulcdat[,ncol(qvulcdat)],-log(qvulcdat[,2]),main="limma",ylab="-log(q-value)",cex=0.2)
    abline(-log(0.05),0)
    plot(qvulcdat[,ncol(qvulcdat)],-log(qvulcdat[,3]),main="rank products",xlab="log(ratio)",ylab="-log(q-value)",cex=0.2,ylim=c(1,maxy))
    abline(-log(0.05),0)
    # dev.off()
    colnames(qvulcdat) <- c("t-test","limma","rank products","mean log-ratio")
    OutData[[i]]<-qvulcdat 
    
  }
    par(mfrow=c(1,1))
  
  return(OutData)
}

