#############################################################
# Peptide rollup functions (reference based) used in DAnTE  #
#                                                           #
# By Ashoka D. Polpitiya     
# further changes by V. Schwammle                               #
#############################################################
# minPresence: counts for the most replicated peptide, not all of the individually!!!
# minOverlap: just not include peptide for scaling but still use for median/mean calculation -> problematic!
## Changed: 
## new parameter minPres: now removes all peptides with a presence below the threshold!!!!!
## new parameter minPep: removes all proteins with less than minPep proteins

                     
RRollup.dialog <- list(
  label = "\nCrosstab Rollup. Select a Data Crosstab with associated Row Metadata and the field to roll up.\n",
  show.progress=TRUE,
  Data.dataframeItem="T_Data", label="Data Source (must be log transformed)", 
    tooltip = "This method assumes that the data is in log scale.",
  signal = c("default", "get.dataset.row.metadata.fields", "rollup.field"),
  rollup.field.choiceItem = NULL, label = "Choose Field Containing Proteins",
  minPresence.numericItem=50, label = "Minimum presence % of at least one peptide",
    tooltip="Peptides with less values present than this percentage are dropped.",  
  Mode.radiobuttonItem=c(value="median", "mean"), labels=c("Median", "Mean"), label="Mode", BREAK=T,
  minOverlap.integerItem=3, label="Minimum dataset presence",                       
    tooltip = "Within a group of peptides for a specific protein, the ones that do not overlap well (controlled by this value) are not scaled but they are kept to calculate the final protein abundance.",    
  oneHitWonders.trueFalseItem=FALSE, label="Include One Hit Wonders",
   tooltip = "Protein with only one observed peptide will be included in the final list of proteins. The rationale behind this is that if a particular protein may have only one peptide but it may be quite abundant and present throughout giving some strong confidence on the presence of the protein.",
  reportCount.trueFalseItem=TRUE, label="Report Count Table",  
    tooltip = "Include a table containing peptide counts per protein in the output",
  gpvalue.numericItem=0.05, label="p-value Cutoff for Grubbs' Test",
  gminPCount.integerItem=5, label="Minimum Number of Peptides required for Grubbs' Test",
    tooltip = "Grubb's test for outliers is performed for peptide groups in each dataset. This value controls the minimum number required to perform the test." ,
  center.trueFalseItem=FALSE, label="Mean Center Proteins to Zero Mean"
)

##############################################################################
CountRollup <- function(MSdat1, ProtInfo){
  myTfrm <- function(x) {
    rownames(x) <- x[,1]
    x[,-1,drop=F]
  }
  library(reshape)
  print("Reshaping data..."); flush.console()
  MSdat <- cbind(ProtInfo, MSdat1)
  print("Melting data..."); flush.console()
  cc.melt <- melt(MSdat, id.vars = c("Reference", "Mass Tag ID"), variable_name = "Job")
  print("Data melt done. Casting..."); flush.console()
  MSpivot <- cast(cc.melt, Reference~Job, function(x) sum(!is.na(x)))
  print("Data cast done."); flush.console()  
  MSpivot <- myTfrm(MSpivot)
  MSpivot <- cbind(Total = MSpivot[,1], MSpivot[,colnames(MSdat1),drop=F])
  return(MSpivot)
}

##############################################################################
RRollup <- function(Data, Accs, rollup.field=NULL, minPresence=50, Mode="median",minPep=3,minPres=50,
                        minOverlap=3, oneHitWonders=TRUE, reportCount=FALSE, outfolder="C:/",
                        plotflag=FALSE, gpvalue=0.05, gminPCount=5, center=TRUE,
                        progressbar=NULL, progresslabel=NULL
                        )
#
# Ashoka Polpitiya - July 2007
#
{

    ## Veit: filtering for number of peptides and presence
   Accs<-Accs[rowSums(!is.na(Data))>(minPres/100*ncol(Data))]
    Data<-Data[rowSums(!is.na(Data))>(minPres/100*ncol(Data)),]
  

    # from metadata
    #ProtInfo <- get.ProtInfo(Data)
    #if(!rollup.field%in%colnames(ProtInfo)) stop("Rollup field not found in peptide metadata table")
    data_attributes <- attributes(Data)

    #ProtInfo <- unique(ProtInfo)
    MassTags <- rownames(Data)
    protIPI <- Accs
    

  ## Veit: include parameter for min. number of peptides
    if (oneHitWonders){ 
      minCountPerP <- 2
    }else{
      minCountPerP <- minPep
    }
    minPresence <- minPresence/100

    uProtCounts <- table(protIPI)
    sigIPI <- names(uProtCounts[uProtCounts >= minCountPerP])
    if(!length(sigIPI)) stop(paste("No proteins were seen with more than", minCountPerP, "peptides"))
    restIPI <- names(uProtCounts[(uProtCounts < minCountPerP) & (uProtCounts > 0)])
    #print(restIPI)
    threshold <- round(dim(Data)[2] * minPresence)

    scaledData <- rep(numeric(0),dim(Data)[2])
    orData <- rep(numeric(0),dim(Data)[2])
    protData <- rep(numeric(0),dim(Data)[2]+1)
    sdsData<-NULL
    singlePepProtData <- rep(numeric(0),dim(Data)[2]+1)
    
    # Check that all mass tags are matched
    if(!any(row.names(Data)%in%MassTags)) stop("ProtInfo first column doesn't match any row names of data")
    if(!all(row.names(Data)%in%MassTags)) stop("ProtInfo first column doesn't match all row names of data")

    proteinNames <- "empty"
    oneHitProtNames <- "empty"    

    k = 1
    Nprots <- length(sigIPI)
    for (prot in 1:length(sigIPI)){
    if(prot%%10==0) {
      if(!missing(progressbar)) progressbar$setFraction(prot/Nprots)
      if(!missing(progresslabel)) progresslabel$setText(paste("Status:", prot, "/", Nprots))      
    }
      #print(prot); flush.console()
        pidx <- which(sigIPI[prot]==protIPI)
        data_idx <- is.element(row.names(Data),MassTags[pidx])
        currProtData <- Data[data_idx,,drop=F]
        
        if (is.matrix(currProtData) || is.data.frame(currProtData))
            if (dim(currProtData)[1] > 1)
            {
                xPresenceCount <- rowSums(!is.na(currProtData), na.rm=TRUE)
                if (max(xPresenceCount, na.rm=TRUE) >= threshold)
                {

                    pData <- protein.rollup1(currProtData,minOverlap=minOverlap,
                             Mode=Mode, minPs=gminPCount, pvalue=gpvalue, center=center)
                    scaledData <- rbind(scaledData,pData$sData)
                    orData <- rbind(orData,pData$orData)
                    protData <- rbind(protData,pData$pData)
                    sdsData <- rbind(sdsData,pData$sds)
                    proteinNames[k] <- sigIPI[prot]

                    if (plotflag)
                    {
                        outfile = paste(outfolder,k,".png",sep="")
                        plotCurrProt.RefRup(currProtData, pData, file=outfile,IPI=sigIPI[prot])
                    }
                    k = k + 1
                }
            }
    }
    rownames(protData) <- proteinNames
    scaledData <- remove.duplicates(scaledData)
    orData <- remove.duplicates(orData)
    
    if (oneHitWonders)
    {
        k = 1
        for (prot in 1:length(restIPI))
        {
            pidx <- which(restIPI[prot]==protIPI)
            data_idx <- is.element(row.names(Data),MassTags[pidx[1]])
            currProtData <- Data[data_idx,,drop=F]
            
            xPresenceCount <- sum(!is.na(currProtData))
            if (xPresenceCount >= threshold)
            {
                #singlePepProt <- c(PepCount=1,currProtData)
                singlePepProtData <- rbind(singlePepProtData,currProtData)
                oneHitProtNames[k] <- restIPI[prot]
                k = k + 1
            }
        }
        rownames(singlePepProtData) <- oneHitProtNames
        if (center)
        {
            singlePepProtData1 <- singlePepProtData[,-1]
            scaledSinglePepProtData <- t(scale(t(singlePepProtData1),center=T,scale=F))
            singlePepProtData <- cbind(singlePepProtData[,1],scaledSinglePepProtData)
        }
        outProtData <- rbind(protData,singlePepProtData)
    }else{
        outProtData <- protData
    }
    
    attr(outProtData, "Column_Metadata") <- attr(scaledData, "Column_Metadata") <- attr(Data, "Column_Metadata")
    attr(scaledData, "Row_Metadata") <- attr(Data, "Row_Metadata")
    attr(outProtData, "Row_Metadata") <- c(table=attr(Data, "Row_Metadata")[["table"]], key = rollup.field)
#    out <- list(ScaledData=scaledData, OriginalData=orData, ProteinData=outProtData)
    out <- list(Scaled=scaledData, RolledUp=outProtData, Sds=sdsData)
    if(reportCount) {         
#    uProtCounts <<- uProtCounts
#   outProtData <<- outProtData                                
      countArr <- cbind(Count=uProtCounts[names(uProtCounts)%in%rownames(outProtData)])
#      print(dim(countArr))
      attr(countArr, "Row_Metadata") <- attr(outProtData, "Row_Metadata")
      out <- list(Scaled=scaledData, RolledUp=outProtData, Count = countArr, Sds=sdsData)
    }
#    out <<- out
    return(out)
}                                

##############################################################################
protein.rollup2 <- function(xx, Nmax=5, ...){
  if(nrow(xx) > Nmax)
    xx <- xx[order(apply(xx, 1, median, na.rm=T), decreasing=T)[1:Nmax],,drop=F]

  f_exp <- factor(rep(1:ncol(xx), nrow(xx)))
  f_pep <- factor(rep(1:nrow(xx), each=ncol(xx)))
  fit2 <- rlm(Y ~ f_pep + f_exp - 1, contrasts = list(f_pep = "contr.sum", f_exp = "contr.sum"))
  ff1 <- fit2$coef[nrow(xx)+1:(ncol(xx)-1)]

  median(fit2$coef[1:nrow(xx)]) + c(ff1, -sum(ff1))
}

##############################################################################
protein.rollup1 <- function(currSel, minOverlap=3, Mode="median",
                        minPs=5, pvalue=0.05, center=TRUE)
{
    pepCounts <- rowSums(!is.na(currSel)) # presence of a peptide sample-wise
    pepMaxCounts <- which(pepCounts == max(pepCounts)) # maximally present ones
    if (length(pepMaxCounts) > 1){
        possibleRefs <- currSel[pepMaxCounts,] # multiple refrence peptides?
        totalAbundances <- rowSums(possibleRefs, na.rm=TRUE) # pick the one with highest overall abundance
        refs <- which(totalAbundances == max(totalAbundances,na.rm=TRUE))
        reference <- pepMaxCounts[refs[1]]
    }else{
        reference <- pepMaxCounts[1]
    }

    currSelAdj <- t(t(currSel) - as.vector(t(currSel[reference,]))) #get ratios
    overlapCount <- rowSums(!is.na(currSelAdj)) # count the number of non-missing
    overlapMedians <- apply(currSelAdj,1,median,na.rm=T) # median ratios for each peptide
    overlapMedians[which(overlapCount < minOverlap)] <- 0 # get rid of medians for sparse peptides
    currSel <- currSel - overlapMedians # adjust the originals with median ratios

    orCurrSel <- remove.outliers(currSel, minPs=minPs, pvalue=pvalue)

    if (center)
    {
        scaled <- t(scale(t(currSel),center=T,scale=F))
        scaled.or <- t(scale(t(orCurrSel),center=T,scale=F))
    }
    else
    {
        scaled <- currSel
        scaled.or <- orCurrSel
    }

    if (Mode=="median")
        proteinVal <- apply(scaled.or, 2, median, na.rm=T) # rollup to proteins as medians
    else
        proteinVal <- apply(scaled.or, 2, mean, na.rm=T)

    sdVal <- apply(scaled.or, 2, sd, na.rm=T)

    #proteinVal <- c(PepCount=dim(scaled.or)[1],proteinVal) # append the column with peptide counts
    out <- list(sData=scaled, orData=scaled.or, pData=proteinVal,sds=sdVal)
    return(out)
}

##############################################################################
remove.outliers <- function(Data, minPs=5, pvalue=0.05)
# internal function used by normalize.proteins()
# Calculates the protein abundances after removing outliers
# Depends on package "outliers"
{
    library(outliers)
    ColNames = colnames(Data)
    xPeptideCount <- colSums(!is.na(Data))
    for (i in 1:dim(Data)[2])
    {
        if (xPeptideCount[i] >= minPs)
        {
            repeat
            {
                grubbs <- grubbs.test(Data[,i]) # Grubb's test
                if ( (grubbs$p.value < pvalue) && (!is.nan(grubbs$statistic[2])) &&
                        (grubbs$statistic[2] != 0) ) # pass the p-value cutoff
                {
                    Data[,i] <- rm.outlier.1(Data[,i],fill=TRUE,median=TRUE)
                    # fill the outlier with the median value
                }
                else { break }
            }
        }
    }
    return(Data)
}


##############################################################################
plotCurrProt.RefRup <- function(currData,pdata,file="deleteme.png",
                            bkground="transparent",
                            IPI="IPI:IPI00009793.1")
# pdata : output from protein rollup (scale.proteins)
{
    require(Cairo)
    CairoPNG(filename=file,width=IMGwidth,height=IMGheight,pointsize=FNTsize,bg=bkground,res=600)
    #png(filename=file,width=1024,height=768,pointsize=12,bg=bkground,
    #        res=600)
    tryCatch(
    {
    scaledData <- pdata$sData
    protData <- pdata$pData
    orData <- pdata$orData

    par(mfrow=c(3,1))
    matplot(t(currData),type="b",main=IPI,ylab="Raw Data")
    matplot(t(scaledData),type="b",ylab="Scaled Data")
    matplot(t(orData),type="b",ylab="Scaled and Outlier removed",
            xlab=paste(dim(currData)[1],"Peptides",sep=" "))
    lines(protData[-1],type="l",lwd=2)
    },
    interrupt = function(ex)
    {
      cat("An interrupt was detected.\n");
      print(ex);
    },
    error = function(ex)
    {
      cat("An error was detected.\n");
      print(ex);
    },
    finally =
    {
      cat("Releasing tempfile...");
      dev.off()
      cat("done\n");
    }) # tryCatch()
}


###############################################################################
# modified R outier functions from "outlier" package to handle missing
rm.outlier.1 <- function (x, fill = FALSE, median = FALSE,
                         opposite = FALSE,
                         na.rm = TRUE)
{
    if (is.matrix(x))
        apply(x, 2, rm.outlier.1, fill = fill, median = median,
            opposite = opposite, na.rm = na.rm)
    else if (is.data.frame(x))
        as.data.frame(sapply(x, rm.outlier.1, fill = fill, median = median,
            opposite = opposite, na.rm = na.rm))
    else {
        res <- x
        if (!fill)
            res[-which(x == outlier.1(x, opposite))]
        else {
            if (median)
                res[which(x == outlier.1(x, opposite))] <- median(x[-which(x ==
                  outlier.1(x, opposite))], na.rm = na.rm)
            else res[which(x == outlier.1(x, opposite))] <- mean(x[-which(x ==
                outlier.1(x, opposite))], na.rm = na.rm)
            res
        }
    }
}

##############################################################################
# modified R outier functions to handle missing
outlier.1 <- function (x, opposite = FALSE, logical = FALSE, na.rm = TRUE)
{
    if (is.matrix(x))
        apply(x, 2, outlier.1, opposite = opposite, logical = logical)
    else if (is.data.frame(x))
        sapply(x, outlier.1, opposite = opposite, logical = logical)
    else {
        if (xor(((max(x, na.rm = na.rm) - mean(x, na.rm = na.rm)) <
            (mean(x, na.rm = na.rm) - min(x, na.rm = na.rm))), opposite)) {
            if (!logical)
                min(x, na.rm = na.rm)
            else x == min(x, na.rm = na.rm)
        }
        else {
            if (!logical)
                max(x, na.rm = na.rm)
            else x == max(x, na.rm = na.rm)
        }
    }
}
###############################################################################

remove.duplicates <- function(Data)
{
     rowIDs <- rownames(Data)
     Data <- Data[!duplicated(rowIDs),]
     return(Data)
}


#------------------------------------------------------------------

