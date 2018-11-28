# # function to calculate p-values for a data set consisting of at least 3 replicates and missing values
# until now it will only look for relative regulations within ratios
# the algorithm will compare limma, rank products and t-test
# see ...
# Input values are:
# Data: Data set consisting of ratios for abundance values for different replicates and conditions
# Reps: vector with specification of ratios beginning with 1 ; (1,2,3,1,2,3) means replicate 1 ratio 1, replicate 1 ratio 2 ,...
#

LimmaRankProd <- function(Data,Reps) {

library(multtest)
library(limma)
library(genefilter)
library(qvalue)
library(gplots)

## Generate MA matrix
MAData<-Data
MAReps<-Reps

NumCond<-max(Reps)
NumReps<-min(table(Reps))

print(paste("NumCond: ",NumCond,"NumReps: ",NumReps))

##limma with ratios
design<-NULL
plvalues<-NULL
for (c in (1:NumCond)) {
  design<-cbind(design,as.numeric(MAReps==c))
}
#design<-cbind(i114vs113=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),i115vs113=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0),i116vs113=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1))
lm.fittedMA <- lmFit(MAData,design)
lm.bayesMA<-eBayes(lm.fittedMA)
topTable(lm.bayesMA)
if (!is.null(dim(lm.bayesMA$p.value))) {
  for (i in 1:NumCond) 
    plvalues[[i]]<-lm.bayesMA$p.value[,i]
} else {
#   hist(lm.bayesMA[["p.value"]],breaks=100)
  plvalues[[1]]<-lm.bayesMA$p.value
#   qqt(lm.bayesMA$t,df=lm.bayesMA$df.residual+lm.bayesMA$df.prior)
}
###################

## MA t-test_pvalues
ptvalues<-NULL
pRPvalues<-NULL
for (vs in 1:NumCond) {
  tMAData<-MAData[,MAReps==vs]
  #filter for at least two points
  tMAData<-tMAData[rowSums(is.finite(tMAData))>1 & rowSds(tMAData,na.rm=T)>0,]
  ptMAvalues<-NULL
  for (pep in 1:(dim(tMAData)[1])) {
#     print(tMAData[pep,])
    ptMAvalues<-append(ptMAvalues,t.test(as.vector(tMAData[pep,]))$p.value)
  }
  names(ptMAvalues)<-rownames(tMAData)
  ptvalues[[vs]] <- ptMAvalues
#   hist(ptMAvalues,100)

  ## rank products
  # substitute missing values by mean of other replicates
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
  #  RPownCorr<- -log(RP.own)+d*log(sum(Rank[[r]]>1)+2)
    RPownCorr<- -log(RP.own)#+d*log(RankNAs/NumReps+2)
    hist(RPownCorr,100)
    qqplot(qgamma((1:100)/100,d),RPownCorr)
    abline(0,1)
    RPMAownUp_pvalues<-c(RPMAownUp_pvalues,pgamma(RPownCorr,d))
    hist(RPMAownUp_pvalues,100)
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
  #  RPownCorr<- -log(RP.own)+d*log(sum(Rank[[r]]>1)+2)
    RPownCorr<- -log(RP.own)#+d*log(RankNAs/NumReps+2)
    hist(RPownCorr,100)
    qqplot(qgamma((1:100)/100,d),RPownCorr)
    abline(0,1)
    RPMAownDown_pvalues<-c(RPMAownDown_pvalues,pgamma(RPownCorr,d))
    hist(RPMAownDown_pvalues,100)
  }
  }
  RPMAown_pvalues<-2*apply(cbind(RPMAownDown_pvalues,RPMAownUp_pvalues),1,min)
  RPMAown_pvalues<-RPMAown_pvalues[RPMAown_pvalues<1]
  hist(RPMAown_pvalues,100)
  pRPvalues[[vs]]<-RPMAown_pvalues
}

OutData<-NULL

## q-values
qlvalues<-qtvalues<-qRPvalues<-NULL
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

  par(mfrow=c(2,3))
  hist(ptvalues[[i]],100)
  hist(plvalues[[i]],100)
  hist(pRPvalues[[i]],100)

  qvulcdat<-merge(qtvalues[[i]],qlvalues[[i]],all=T,by=0)
  qvulcdat<-merge(qvulcdat,qRPvalues[[i]],all=T,by.x=1,by.y=0)
  qvulcdat<-merge(qvulcdat,rowMeans(MAData[,Reps==i],na.rm=T),all=T,by.x=1,by.y=0)
  rownames(qvulcdat)<-qvulcdat[,1]
  qvulcdat<-qvulcdat[,2:ncol(qvulcdat)]
  maxy<-max(-log(qvulcdat[,1:3]),na.rm=T)
  plot(qvulcdat[,ncol(qvulcdat)],-log(qvulcdat[,1]),main="t-test",xlab="log(ratio)",ylab="-log(p)",cex=0.2,ylim=c(1,maxy))
  abline(-log(0.05),0)
  plot(qvulcdat[,ncol(qvulcdat)],-log(qvulcdat[,2]),main="limma",xlab="log(ratio)",ylab="-log(p)",cex=0.2)
  abline(-log(0.05),0)
  plot(qvulcdat[,ncol(qvulcdat)],-log(qvulcdat[,3]),main="rank products",xlab="log(ratio)",ylab="-log(p)",cex=0.2,ylim=c(1,maxy))
  abline(-log(0.05),0)
  par(mfrow=c(1,1))
  # dev.off()
  OutData[[i]]<-qvulcdat 
  
}

return(OutData)
}
  
