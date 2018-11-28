#------------------------------------------------------------------------------
# Function  : prep.file                                                        
# Aim       : convert input file to DGE list object in prep for DE analysis    
# Input     :                                                                  
#   file    : name of file that contains raw read counts and is in folder "Input" eg: raw_count.tsv
#   sep     : separator used in input file eg: "\t","\\," etc...Default = "\t"
#   exclude : character string of samples to exclude from analysis. Give unique names as uses grep. Default = ""
#   meta    : columns that contain metadata such as group, lane, library, date, timepoint. Default = 1:4
#   strsp   : String to split column names by. This is then used to make the sample information object eg : "\\." 
# Output    : DGEList object with counts, sample information and gene information
#------------------------------------------------------------------------------

prep.file <- function(file=filename,sepr="\t",exclude=NA,meta=1:4,strsp="\\."){
  
  # Read input file
  temp = read.delim(paste("Input",file,sep="/"),sep=sepr,header=T)
  
  # Create genes vector
  genes = temp[,meta]
  
  # Create the data matrix without any metadata or excluded samples
  # Rounding counts with decimal numbers (why are they fractions ?)
  
  temp2 = temp[,-(meta)]
  data = temp2[,!grepl(paste(exclude,collapse="|"),colnames(temp2))]
  
  # Rename young samples....
  colnames(data) = gsub("^control","young.control",colnames(data))
  colnames(data) = gsub("^cooled","young.cooled",colnames(data))
  colnames(data) = gsub("^rewarmed","young.rewarmed",colnames(data))
  
  rownames(data) = temp$gene_id
  
  # Extracting sample information from column names
  # If a given column has more than one catergory, change it into a variable column
  
  s = strsplit(colnames(data),strsp)
  cnums = NULL
  for(y in 1:max(lengths(s))){
    if(length(unique(sapply(s,"[[",y)))>1){
      #print(y)
      cnums = c(cnums,y)
    }
  }

  # The unique columns numbers are stored in 'cnums'
  sample.list = NULL
  for(k in cnums){
    sample.list = cbind(sample.list,sapply(s,"[[",as.numeric(k)))
  }
  sample.list = as.data.frame(sample.list)
  rownames(sample.list) = apply(sample.list,1,paste,collapse=".")
  #sample.list
  
  # Create a DGEList object using data, genes and samples objects made with the code above
  x <- DGEList(counts = round(data,0),samples=sample.list,genes=genes)
  x
  
  return(x)
}


#------------------------------------------------------------------------------
# Function  : clean.data.qc
# Aim       : clean input data and generate QC plots
# Input     : 
#   x       : DGEList object, usually returned from 
#   feat    : features of interest in x$samples - names or columns. Used to draw boxplots against each factor
#   qc.dir  : name of the directory to which you want to print your output. Default is 'outdir'
# Output    : DGEList object with counts, normalisation and only those genes that are expressed
#------------------------------------------------------------------------------

clean.data.qc <- function(x,feat,qc.dir=outdir,filt){
  
  # Create the QC output directory to make it easier to look through results
  #-------------------------------------------------------------------------
  if (exists(qc.dir)){
    print("QC directory exists")
  }else{
    dir.create(qc.dir)
  }
  
  # If 'feat' is numeric, extract the column names, else leave them as is
  #-------------------------------------------------------------------------
  foi = feat
  if(is.numeric(feat)){
    foi = colnames(x$samples[,feat])
  }
  
  # Plot summary statistics for library size (looking to see if there are group specific library effects)
  #-------------------------------------------------------------------------
  for(h in foi){
    print(paste("Parameter tested is: ",h,sep=""),quote = F,row.names=F)
    print(aggregate(x$samples$lib.size,by=list(x$samples[,h]),FUN=summary))
  }
  
  # Generate boxplots of library size for various paramenters
  #-------------------------------------------------------------------------
  pdf(paste(qc.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Boxplot-of-library-size-by-group.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  par(mfrow=c(1,1))
  for(i in foi){
    boxplot(x$samples$lib.size~x$samples[,i],ylab = "Library size",las=2,xlab = i,main=paste("Comparison of RNAseq library sizes across",i,sep=" "),col=rainbow(length(unique(x$samples[,i]))))
  }
  dev.off()
  
  # Look for rows that have all zero values. These are genes that are not expressed so don't add value to analysis
  # We want to look for rows where >= 70% of samples are not expressed
  table(rowSums(x$counts==0)==(round(0.7*ncol(x$counts),0)))
  
  # Check if some samples have more non-expressing genes than others
  x$samples$Tot.genes = nrow(x$counts)
  x$samples[,"Non.exp.genes"] = colSums(x$counts==0)
  x$samples[,"Exp.genes"] = x$samples$Tot.genes - x$samples$Non.exp.genes
  x$samples[,"Median.non.exp.genes"] = median(colSums(x$counts==0))
  x$samples[,"Non.exp.diff.from.median"] = colSums(x$counts==0)-median(colSums(x$counts==0))
  x$samples[,"Perc.exp.genes"] = round(100*(x$samples$Exp.genes)/x$samples$Tot.genes,2)
  
  # Write the output to file
  write.table(x$samples,paste(qc.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"sample-metadata.txt",sep="_"),sep="/"),quote=F,sep="\t")
  
  
  # Filtering the cpm table for only those genes that are expressed across the group
  # Keep those genes whose count per million > 1 for atleast 3 samples out of 12 which is the size of the smallest group
  # cpm > 1 means that a gene is expressed if there are at least 15 reads in the sample with the smallest library size (A019,14.6M) and at least 43 reads in samples with largest library size (A005,42.8M)
  # A cpm > 1 means log(cpm) > 0
  #-------------------------------------------------------------------------
  
  # Read counts need to be normalised and logged to proceed using limma. Normalisation converts counts to counts-per-million
  prefilt.cpm <- cpm(x)
  prefilt.lcpm <- cpm(x, log=TRUE)
  
  # Save a copy of original data as 'xi' in case you need to go back to it.
  # In the modular version of this script, dge.1 is xi so it is OK not to make a back up....
  xi = x
  
  # Plotting the results of removing non-expressing genes
  #-------------------------------------------------------------------------
  nsamples <- ncol(x)
  coul <- brewer.pal(12, "Paired")
  col = colorRampPalette(coul)(nsamples)
  
  # Using values in the "filt" option
  # If no value is provided for filter, run for a range of values and return plots
  #-------------------------------------------------------------------------
  if(is.null(filt)){
    
    pdf(paste(qc.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Pre-post-filtering-density-plots.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
    for(k in 1:ncol(x)){
      x = xi
      keep.exprs <- rowSums(prefilt.cpm>1)>=k
      x <- x[keep.exprs,, keep.lib.sizes=FALSE]
      
      # Recalculate cpm for filtered data
      postfilt.lcpm <- cpm(x, log=TRUE)
      
      par(mfrow=c(1,2))
      plot(density(prefilt.lcpm[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2, main="", xlab="")
      title(main=paste("A. Raw data (n = ",k," ; ",nrow(xi)," genes)", sep=""), xlab="Log-cpm")
      abline(v=0, lty=3)
      for (i in 2:nsamples){
        den <- density(prefilt.lcpm[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
      }
      legend("topright", rownames(x$samples), text.col=col, bty="n",cex=0.4)
      
      plot(density(postfilt.lcpm[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,main="", xlab="")
      title(main=paste("B. Filtered data (n = ",k," ; ",nrow(x)," genes)", sep=""), xlab="Log-cpm")
      abline(v=0, lty=3)
      for (i in 2:nsamples){
        den <- density(postfilt.lcpm[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
      }
      legend("topright", rownames(x$samples), text.col=col, bty="n",cex=0.4)
    }
    dev.off()
    final.filt = as.numeric(readline(prompt="Review your plots and enter cut-off:"))
    x = xi
    keep.exprs <- rowSums(prefilt.cpm>1)>=final.filt
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  }
  else{
    final.filt = filt
    pdf(paste(qc.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),filt,"Pre-post-filtering-density-plots.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
    x = xi
    keep.exprs <- rowSums(prefilt.cpm>1)>=filt
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    
    # Recalculate cpm for filtered data
    postfilt.lcpm <- cpm(x, log=TRUE)
    
    par(mfrow=c(1,2))
    plot(density(prefilt.lcpm[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2, main="", xlab="")
    title(main=paste("A. Raw data (filter = ",0," ; ",nrow(xi)," genes)", sep=""), xlab="Log-cpm")
    abline(v=0, lty=3)
    for (i in 2:nsamples){
      den <- density(prefilt.lcpm[,i])
      lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", rownames(x$samples), text.col=col, bty="n",cex=0.4)
    
    plot(density(postfilt.lcpm[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,main="", xlab="")
    title(main=paste("B. Filtered data (filter = ",filt," ; ",nrow(x)," genes)", sep=""), xlab="Log-cpm")
    abline(v=0, lty=3)
    for (i in 2:nsamples){
      den <- density(postfilt.lcpm[,i])
      lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", rownames(x$samples), text.col=col, bty="n",cex=0.4)
    dev.off()
  }
  #-------------------------------------------------------------------------
  # Normalisation by trimmed means of M-values (TMM)
  # Used to make the range of expression values similar across all samples
  # Hopefully will remove any unwanted variation by using this method
  # In microarray lingo, M = log of the ratio of expression between red and green channels (tmt and control)
  # If reference column is unspecified, then the sample whose upper quartile is closest to 
  # the mean upper quartile is used as a reference and all others are scaled relative to it
  #-------------------------------------------------------------------------
  x.norm <- calcNormFactors(x, method = "TMM")
  x.norm$samples[order(x.norm$samples$norm.factors),]
  
  # Plotting before and after normalisation
  #-------------------------------------------------------------------------
  plot.lab.x = apply(x$samples[,c(4:6,8)],1,paste,collapse=".")
  
  pdf(paste(qc.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Pre-post-normalisation-box-plots.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  par(mfrow=c(1,1))
  lcpm.prenorm <- cpm(x, log=TRUE)
  boxplot(lcpm.prenorm, las=2, col=col, main="",names=plot.lab.x,cex=0.6,cex.axis=0.8)
  title(main="A. Example: Unnormalised data",ylab="Log-cpm")
  
  lcpm.postnorm <- cpm(x.norm, log=TRUE)
  boxplot(lcpm.postnorm, las=2, col=col, main="",names=plot.lab.x,cex=0.6,cex.axis=0.8)
  title(main="B. Example: Normalised data",ylab="Log-cpm")
  dev.off()
  
  # In an MDS plot, distances on the plot correspond to the leading fold-change, 
  # which is the average (root-mean-square) log2-fold-change for the 500 genes most divergent between each pair of samples by default.
  #-------------------------------------------------------------------------
  pdf(paste(qc.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Pre-post-filtering-MDS-PCA-plots.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  
  mains = c("Pre-filtering","Post-filtering-for-non-expressing-genes")
  c = 1
  # Loop through pre and post filtered data
  for(data in list(xi,x.norm)){
    m = plotMDS(data, pch=c(16,17)[as.numeric(as.factor(data$samples$Rep))], gene.selection = "pairwise",cex=2,col=c("green3","blue3","red3")[as.numeric(as.factor(data$samples$Tmt))],main=mains[c])
    text(x=m$x, y=m$y, labels = names(m$x),cex=0.6)
    legend("topright",legend=unique(data$samples$Tmt),fill = c("green3","blue3","red3"))
    legend("bottomleft",legend=unique(data$samples$Rep),pch = c(16,17))
    c=c+1
  } 
  
  dev.off()
  
  # Find outlier samples - 4 samples with largest variance
  post.mds = plotMDS(lcpm.postnorm,plot=F)
  pp = post.mds$cmdscale.out[order(post.mds$cmdscale.out[,1]),]
  outl = rownames(tail(pp,n=4))
  
  # MDS plot by feature....using post-normalisation data
  pdf(paste(qc.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Post-normalisation-MDS-plots-by-feature.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  par(mfrow=c(1,1))
  for(i in feat){
    getPalette = colorRampPalette(brewer.pal(9, "Set3"))
    col.group <- as.factor(x.norm$samples[,i])
    levels(col.group) <- getPalette(nlevels(col.group))
    col.group <- as.character(col.group)
    par(mfrow=c(1,1))
    plotMDS(lcpm.postnorm, labels=x.norm$samples[,i], col=col.group,main = paste(colnames(x.norm$samples)[i],": Dimensions 1 & 2"))
  }
  dev.off()
  return(list(x.norm,outl))
}

#----------------------------------------------------------------------------
# Function   : runEdgeRsingleF
# Aim        : Run a differential expression analysis on the data using edgeR for a single factor design
# Input      : 
#     x.norm : DGEList object, after TMM normalisation and removal of non-expressing genes
#     des    : design matrix that contains the factor you want to compare
#     suf    : Suffix to be added to output files eg: "Induced.vs.Uninduced"
#     pair   : The two groups that are to be compared. Provide the control group first. eg : c("Uninduced","Induced")
#     edge.dir  : Directory into which you want your edgeR results to go
#     pval   : p-value cut off for adjusted p-values to determine significant genes
#     logfc  : logFC cut off to determine significant genes
# Output     : Text files with most DE genes, GO terms for the comparison.
#              Returns significant genes and top25 GO terms
#----------------------------------------------------------------------------

runEdgeRsingleF <- function(x.norm,des,pair,suf,edge.dir,pval=0.05,logfc=1){
  
  if (exists(edge.dir)){
    print("Outdir exists")
  }else{
    dir.create(edge.dir)
  }
  
  # Using edgeR rather than limma for DE - single factor being Aged or Young
  # By default, it always uses the first column as control which in our case is "aged" mice
  # Therefore, the comparison is going to be "Young" - "Aged". 
  # Single factor, two group comparison so can use exactTest and no need for glm
  
  disp = estimateDisp(x.norm,design=des)
  pdf(paste(edge.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),suf,"edgeR-plots-showing-biological-COV.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  plotBCV(disp)
  dev.off()
  
  # Run and exact test (based on Fisher's for DE)
  et.y.all = exactTest(disp, pair=pair,dispersion = disp$common.dispersion)
  dt.y = decideTestsDGE(et.y.all)
  colnames(dt.y) = suf
  rownames(dt.y) = rownames(x.norm)
  top = as.data.frame(topTags(et.y.all,n=Inf,p.value=pval))
  all = as.data.frame(topTags(et.y.all,n=Inf,p.value=Inf))
  write.table(all,paste(edge.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste(suf,"edgeR-based-DE-test-results.txt",sep="_"),sep="_"),sep="/"),sep="\t",quote=F,row.names=F)
  
  return(list(et.y.all,top))
}


#------------------------------------------------------------------------------
# Function  : plotDE
# Aim       : Draw plots of DE genes 
# Input     : 
#   results : A results object from exactTest (DGEExact object), glmLRT(DGEGLM object) or lmFit (MArrayLM object)
#   x.norm  : A DGEList object ocntaining normalised count data and genes/samples metadata
#   logfc.col : Column name that contains log fold change values
#   cpm.col : column that contains counts per million
#   fdr.col : column that contains adjusted p.values - doesn't have to be FDR corrected
#   logfc   : The logFC value to be used to determine genes of interest
#   pval    : The pvalue cut-off to be used to determine genes of interest
#   suf     : suffix for filename
#   out.dir : Directory into which you want your results to go
# cont.col.names : Samples that are involved in the current comparison. Eg : "y.cooled","y.control"
# Output    : Text files with most DE genes for each comparison as well as volcano plots, heatmaps
#------------------------------------------------------------------------------

plotDE <- function(results,x.norm,logfc.col,cpm.col= NULL,pval.col,fdr.col,pval=0.05,logfc=1,suf,out.dir=out.dir,cont.col.names){
  
  pdf(paste(out.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste(suf,pval,logfc,"edgeR-based-cloud-volcano-plots.pdf",sep="_"),sep="_"),sep="/"),paper="a4r",width=12,height=8)
  
  # Generating a test results object for cloud plots. Needs to be decideTests for limma and decideTestsDGE for edgeR 
  if(grepl("limma",suf)){
    
    # Run 'decideTests' for all contrasts
    dt.z = decideTests(results)
    
    # Get column number for contrast of interest
    p = grep(strsplit(suf,"_")[[1]][2],colnames(results))
    
    # Pick top genes with p-values <= pval
    top = as.data.frame(topTable(results,coef=p,n=Inf,p.value=pval))
    
    # Collect all genes for a given contrast
    all = as.data.frame(topTable(results,coef=p,n=Inf,p.value=Inf))
  }else{
    dt.z = decideTestsDGE(results)
    rownames(dt.z) = rownames(results)
    colnames(dt.z) = suf
    
    # Drawing cloud plots
    plotWithHighlights(results$table[,cpm.col],results$table[,logfc.col],status = dt.z, hl.col=c("red3","blue3"),xlab="Average-log-CPM",ylab="logFC",main=paste("Cloud plot for",suf,sep=" "))
    abline(h=c(-1, 1), col="green3",lty=2,lwd=2)
    
    # Pick top genes with p-values <= pval
    top = as.data.frame(topTags(results,n=Inf,p.value=pval))
    
    # Highlight top genes on cloud plots
    with(subset(top,abs(get(logfc.col))>=logfc), textxy(get(cpm.col),get(logfc.col), labs=sym, cex=.6,offset=0.6))
    
    # Collect all genes for a given contrast
    all = as.data.frame(topTags(results,n=Inf,p.value=Inf))
  }
  
  #---------------------
  # Draw volcano plots
  #---------------------
  
  # Log10(P-value) equivalent of the highest adj.P.Value <=0.05
  log.edge = max(-log10(all[which(all[,fdr.col] == max(all[which(all[,fdr.col]<=pval),fdr.col])),pval.col]))
  print(paste("-log10(P-value) corresponding to adj.p.val of 0.05 is",round(log.edge,2),"for this sample.",sep=" ")) 
  
  # Drawing volcano plot
  # If FDR p-value = 0, then -log10(FDR) = Infinity and cannot be plotted, so need to correct for that
  if(min(all[,fdr.col]) == 0){
    all[which(all[,fdr.col]==0),fdr.col] = min(all[-which(all[,fdr.col] == 0),fdr.col])/10
    all[which(all[,pval.col]==0),pval.col] = min(all[-which(all[,pval.col] == 0),pval.col])/10
  }
  
  with(all, plot(get(logfc.col), -log10(get(pval.col)),ylab="-log10(P.Value)",xlab="logFC", pch=20, main=paste("Volcano plot : DE genes for ",suf,"comparison",sep=" ")))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  abline(h=log.edge,col="red3",lty=3)
  abline(v=logfc,col="purple",lty=1)
  abline(v=-logfc,col="purple",lty=1)
  with(subset(all, get(fdr.col)<=pval ), points(get(logfc.col), -log10(get(pval.col)), pch=20, col="blue3")) # Significant p-value
  with(subset(all, abs(get(logfc.col))>=logfc), points(get(logfc.col), -log10(get(pval.col)), pch=20, col="orange")) # Large fold change
  with(subset(all, get(fdr.col)<=pval & abs(get(logfc.col))>logfc), points(get(logfc.col), -log10(get(pval.col)), pch=20, col="green3")) # Both
  
  legend("topleft",c("Adj.P.Value <= 0.05",paste("abs(logFC) >",logfc),paste("Adj.P.Value<=0.05 & abs(logFC) >",logfc),"Not significant"),fill=c("blue3","orange","green3","black"))
  
  # Label a subset of points with the textxy function from the calibrate plot
  with(subset(all, get(fdr.col)<=pval & abs(get(logfc.col))>logfc), textxy(get(logfc.col), -log10(get(pval.col)), labs=sym, cex=.6,offset=0.6))
  
  # Add text for categorisation
  up.label = strsplit(suf,"vs")[[1]][1]
  text((max(all[,logfc.col])-0)/2,max(-log10(all[,pval.col])),paste("Up in",up.label,"mice",sep=" "),col = "blue3")
  text((0+min(all[,logfc.col])/2),max(-log10(all[,pval.col])),paste("Down in",up.label,"mice",sep=" "),col = "orange2")
  
  # Extracting the top (adj.p.val < 0.05) most DE genes for all comparisons
  topgenes <- all$sym[which(all[,fdr.col]<=pval & abs(all[,logfc.col]) >= logfc)]
  k <- which(results$genes$sym %in% topgenes)
  
  # If number of genes 'k' is < 2, do not proceed any further.
  if(length(k) > 2){
  
    # Set up colour vector for celltype variable
    mypalette <- brewer.pal(11,"PRGn")
    morecols <- colorRampPalette(rev(mypalette))
    
    # Subset data for heatmaps from counts matrix
    hm = x.norm[which(x.norm$genes$sym %in% top$sym),]
    hm.cpm = cpm(hm, normalized.lib.sizes=TRUE, log=T)
    hmap.dat = as.matrix(hm.cpm)
    #rownames(hmap.dat) = rownames(x.norm$counts)[k]
    
    # Include only those samples that are part of this comparison.
    #colnums = grep(cont.col.names,colnames(hmap.dat))
    #hmap.dat.sub = hmap.dat[,colnums]
    hmap.dat.sub = hmap.dat
    
    # Add a column colour option to distinguish the two comparison groups
    col.cols = c("blue","darkorange")[as.numeric(as.factor(sapply(strsplit(as.character(colnames(hmap.dat.sub)),"\\.[0-9]+"),"[[",1)))] # Column side colours
    
    # Finally, draw the heatmap
    heatmap.2(hmap.dat.sub,scale='row',Rowv=T,Colv=F, ColSideColors = col.cols, col= morecols(50),labCol=rownames(x.norm$samples),cexCol=1.5,cexRow=0.5,margin=c(8,6), lhei=c(2,10), dendrogram="none",trace="none",key=T, density.info="none",main=paste("Top DE genes across",suf,sep=" "))
    }else{
      print("Not enough data for heatmap")
    }
  dev.off()
}
  
#------------------------------------------------------------------------------
# Function  : runDE
# Aim       : Run a differential expression analysis on the data using edgeR or limma
# Notes     : edgeR uses count data directly rather than conevrting to logCPM. 
#             Hopefully, gives more power to analysis i.e more DE genes
# Input     : 
#   x.norm  : DGEList object, after TMM normalisation and removal of non-expressing genes
#   des     : design matrix for DE analysis containing factor of interest
#   contr.matrix: Matrix of contrasts that you are interested in
#   res.dir : Directory into which you want your results to go
#   logfc   : The logFC value to be used to determine genes of interest
#   pval    : The pvalue cut-off to be used to determine genes of interest
# test.name : which analysis you want to run on your data - edgeR or limma ? 
# Output    : Text files with most DE genes for each comparison as well as volcano plots, heatmaps
#------------------------------------------------------------------------------

runDE <- function(x.norm,des,contr.matrix,res.dir=out.dir,logfc=1,pval=0.05, test.name="edgeR"){
  
  # This function helps switch between edgeR and limma tests for DE
  # Whatever your last action is will be returned to "fit"
  
  switch(test.name,
         edgeR = {
           print("Running edgeR analysis")
           disp = estimateDisp(x.norm, design = des)
           pdf(paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"edgeR-plots-showing-biological-COV.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
           plotBCV(disp)
           dev.off()
           gfit <- glmFit(disp,design = des)
           
           edgefit = NULL
           sum.fit = NULL
           
           for(c in 1:ncol(contr.matrix)){
             # Which contrast do you want the DE run for ?
             cont.name = colnames(contr.matrix)[c]
             up.label = strsplit(cont.name,"vs")[[1]][1]
             
             # Perform a likelihood ratio test
             # Correct for family-wise error
             lrt <- glmLRT(gfit,contrast = contr.matrix[,c])
             rownames(lrt$genes) = lrt$genes$gene_id
             lrt$table$adjPVal = p.adjust(lrt$table$PValue,method = "BH")
             
             cont.col.names = paste(names(which(contr.matrix[,c] != 0)),collapse="|")
             temp = plotDE(lrt,x.norm,"logFC","logCPM","PValue","adjPVal",pval=0.05,logfc=1,paste("edgeR",cont.name,sep="_"),res.dir, cont.col.names)
             
             # Identify significant genes
             dt = decideTestsDGE(lrt,adjust.method = "BH")
             colnames(dt) = cont.name
             rownames(dt) = rownames(lrt$table)
             edgefit = cbind(edgefit,dt)
             
             # Summarise DE
             sum.edge = summary(dt)
             sum.fit = cbind(sum.fit,sum.edge)
             
             # Write out topDE genes
             edge.top = data.frame(topTags(lrt,n=Inf,p.value = Inf))
             write.table(edge.top,paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste(cont.name,pval,"edgeR-based-DE-test-results.txt",sep="_"),sep="_"),sep="/"),sep="\t",quote=F,row.names=F)
           }
         },
         limma = {
           print("Running limma analysis")
           pdf(paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Heteroscedascity-VOOM-expression-vs-variance-plots.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
           par(mfrow=c(1,2))
           
           # Convert count data to logCPM
           v.norm <- voom(x.norm, des, plot=T)
           
           # Fit a linear regression model to the data with Bayesian correction
           vfit <- lmFit(v.norm, des)
           vfit.contr <- contrasts.fit(vfit,contrasts=contr.matrix)
           efit <- eBayes(vfit.contr)
           plotSA(efit)
           title("Final model : Mean-variance trend")
           dev.off()
           
           # Look at the summary of differentially expressed genes based on adj.p.val cut-off
           edgefit = decideTests(efit)
           sum.fit = summary(edgefit)
           
           for(c in 1:ncol(contr.matrix)){
             
             # Name of the comparison
             cont.name = colnames(contr.matrix)[c]
             
             # Primary group of interest
             up.label = strsplit(cont.name,"vs")[[1]][1]
             
             # Samples in the comparison group
             cont.col.names = paste(names(which(contr.aged.tmt[,c] != 0)),collapse="|")
             
             # Write out topDE genes
             lim.top = data.frame(topTable(efit,coef=c,n=Inf,p.value = Inf))
             write.table(lim.top,paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste(cont.name,pval,"limma-based-DE-test-results.txt",sep="_"),sep="_"),sep="/"),sep="\t",quote=F,row.names=F)
             
             # Plot DE genes
             temp = plotDE(efit,x.norm,"logFC",NULL,"P.Value","adj.P.Val",pval=0.05,logfc=1,paste("limma",cont.name,sep="_"),res.dir, cont.col.names)
           }
         },
         stop("Enter something that switches me!")
  )
  sum.edge.file = paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste("Summary-of-",test.name,"-results.txt",sep=""),sep="_"),sep="/")
  write.table(sum.fit,sum.edge.file,row.names=F,quote=F,sep="\t") 
  return(edgefit)
}

#------------------------------------------------------------------------------
# Function  : prep.goseq
# Aim       : Prepare data for goseq on a list of genes obtained from DE analysis
# Input     : 
#   v       : An object returned from running eBayes or edgeR (TestResults matrix). 
#           : Alternatively, this can be a list resulting from running 'getVenn'
#   x.norm  : DGEList object that contains normalised data counts and metadata
#   org.kegg: Mapping of KEGG ids to their pathway descriptions obtained usign the function 'mapPathwayToName' 
#   suf     : Suffix to add to the output file name
# Output    : File containing gene lists for various overlaps
#------------------------------------------------------------------------------

prep.goseq <- function(v,x.norm,org.kegg,pval){
  res.go = NULL
  res.kk = NULL
  sum.up = 0
  sum.down = 0
  
  if(is.list(v)){
    for(i in 1:length(v)){
      g = v[i][[1]]
      suf = names(v)[i]
      list.enrich = enrich.go.kegg(x.norm,g,dirx="na",suf,gen.build = "mm10", gene.id.type = "ensGene",pval=pval,mmu.kegg)
      res.go = rbind(res.go,list.enrich[[1]])
      res.kk = rbind(res.kk,list.enrich[[2]])
    }
  }else{
    for(i in 1:ncol(v)){
      suf = colnames(v)[i]
      g.up = rownames(v)[which(v[,i]==1)]
      g.down = rownames(v)[which(v[,i]==-1)]
      print(paste(length(g.up),length(g.down)))
      sum.up = sum.up+length(g.up)
      sum.down = sum.down+length(g.down)
      if(length(g.up) > 0)
        res.up = enrich.go.kegg(x.norm,g.up,dirx="up",suf,gen.build = "mm10", gene.id.type = "ensGene",pval=pval,mmu.kegg)
      if(length(g.down) > 0)
        res.down = enrich.go.kegg(x.norm,g.down,dirx="down",suf,gen.build = "mm10", gene.id.type = "ensGene",pval=pval,mmu.kegg)
      
      res.go = rbind(res.go,res.up[[1]],res.down[[1]])
      res.kk = rbind(res.kk,res.up[[2]],res.down[[2]])
    }
  }
  print(paste(sum.up,sum.down))
  return(list(res.go,res.kk))
}

#------------------------------------------------------------------------------
# Function  : enrich.go.kegg
# Aim       : Run goseq on a list of genes obtained from DE analysis using GO and KEGG
#           : Additionally, run 'enrichKEGG' from clusterProfiler
# Input     : 
#   x.norm  : A DGEList object that contains normalisation factors for the data
#   g       : character vector of Ensembl ids for genes that are up/downregulated
#   dir     : Are the genes "up" or "down" regulated
#   suf     : Name of the comparison for which enrichment is being performed eg: young.cool.vs.ctrl
#   gen.build : Which genome build do you want to use for mapping GO terms ? 
#   gen.id.type : What is the type of ID you are using to map GO terms (values: knownGene (UCSC), vegaGene, ensGene, geneSymbol)
#   pval    : What is the p-value to determine significant GO terms/categories ? 
#   org.kegg: Mapping of KEGG ids to their pathway descriptions obtained usign the function 'mapPathwayToName' 
# Output    : list with 2 elements - GO enrichment and KEGG enrichment terms, significant only
#------------------------------------------------------------------------------


enrich.go.kegg<-function(x.norm,g,dirx,suf,gen.build = "mm10", gene.id.type = "ensGene",pval=0.05,org.kegg){
  
  # Enrichment performed using GOseq and clusterProfiler separately
  enrich = NULL 
  kk = NULL 
  
  # Form the vector for testing DE genes against the gene universe
  all.genes = rep(0,nrow(x.norm))
  names(all.genes) = rownames(x.norm)
  all.genes[which(names(all.genes) %in% g)] = 1
  print(table(all.genes))
  
  # Calculate probability weighting function for set of genes of interest
  pwf = nullp(all.genes,gen.build,gene.id.type,plot.fit = F)
  #plotPWF(pwf,main=paste(suf,dirx,sep="_"))
  
  # Apply goseq whose basis is a Wallenieus distribution which is a hypergeometric distribution where items are samples with a bias
  # The bias in RNAseq data comes from the fact that genes of larger length seem to contribute more to diff expression than those that are shorter
  # Apply multiple testing corrections as 10s to 100s of test are being performed
  # Can subset significant outcomes but for most of our lists, they won't be.
  
  #------------------------------
  # goseq using GO annotations
  #------------------------------
  GO.wall = goseq(pwf,gen.build,gene.id.type)
  GO.wall$over_represented_adj.p.val = p.adjust(GO.wall$over_represented_pvalue,method = "BH")
  GO.wall = GO.wall[,c(1,6:7,4:5,2,8)]
  GO.wall$comparison = suf
  GO.wall$direction = dirx
  
  # Filtering for DE terms. If no DE terms, don't add
  GO.enriched = GO.wall[which(GO.wall$over_represented_adj.p.val<= pval),]
  if(nrow(GO.enriched) > 0){
    enrich = rbind(enrich,GO.enriched)
  }
  
  # Do we need to correct for length bias ?
  #GO.nobias = goseq(pwf, gen.build, gene.id.type, method = "Hypergeometric")
  #GO.nobias$over_represented_adj.p.val = p.adjust(GO.nobias$over_represented_pvalue,method = "BH")
  
  # Looking at contribution of length bias
  # Look at the plot to assess if the y-axis value is consistently greater than the x-axis value, then there is length bias
  #plot(log10(GO.wall[, 2]), log10(GO.nobias[match(GO.nobias[, 1],GO.wall[, 1]), 2]), ,main = paste(suf,dirx,"bias-vs-no.bias",sep="_"), xlab = "log10(Wallenius p-values)", ylab = "log10(Hypergeometric p-values)",xlim = c(-3, 0), ylim = c(-3, 0))
  #abline(0, 1, col = 3, lty = 2, lwd=3)
  
  #------------------------------
  # goseq using KEGG annotations
  #------------------------------
  KEGG = goseq(pwf, gen.build, gene.id.type, test.cats = "KEGG")
  KEGG$term = org.kegg[KEGG$category,1]
  KEGG$ontology = "KEGG"
  KEGG$category = paste("mmu",KEGG$category,sep="")
  
  # Filtering for DE terms. If no DE terms, don't add
  KEGG$over_represented_adj.p.val = p.adjust(KEGG$over_represented_pvalue,method = "BH")
  KEGG = KEGG[,c(1,6:7,4:5,2,8)]
  KEGG$comparison = suf
  KEGG$direction = dirx
  
  # Filtering for DE terms. If no DE terms, don't add
  KEGG.enriched = KEGG[which(KEGG$over_represented_adj.p.val<=pval),]
  if(nrow(KEGG.enriched) > 0){
    enrich = rbind(enrich,KEGG.enriched)
  }
  
  #-------------------------------------------
  # enrichKEGG analysis using clusterProfiler
  #-------------------------------------------
  
  # Setting up the gene universe - needs to include DE and non-DE genes
  univ = names(all.genes)
  univ.ids.names = bitr(univ, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL","GENENAME"), OrgDb="org.Mm.eg.db")
  
  # Subsetting the DE genes from the universe
  ids.names <- univ.ids.names[which(univ.ids.names$ENSEMBL %in% g),"ENTREZID"]
  
  # Performing functional enrichment using the DE list(g) and universe list as background
  kk.temp <- data.frame(enrichKEGG(gene = ids.names,universe = univ.ids.names$ENTREZID,organism = 'mmu',pvalueCutoff = pval))
  
  # Add extra annotation
  if(nrow(kk.temp)>0){
    kk.temp$comparison = suf
    kk.temp$direction = dirx
    kk.temp$symbol = unlist(lapply(kk.temp$geneID,function(x) paste(univ.ids.names$SYMBOL[which(univ.ids.names$ENTREZID %in% unlist(strsplit(x,"/")))],collapse=";")))
    kk.temp$description = unlist(lapply(kk.temp$geneID,function(x) paste(univ.ids.names$GENENAME[which(univ.ids.names$ENTREZID %in% unlist(strsplit(x,"/")))],collapse=";")))
  }
  
  # Add to kk object
  kk = rbind(kk,kk.temp)
  return(list(enrich,kk))
}

#----------------------------------------------------------------------------------------------
# Aim : Given two experiments of interest, find out how many genes overlap in up, down and neutral category and annotate them functionally
# Input :
#   x.norm  : DGEList object post normalisation
#   dt      : A decideTests object that contains experiments of interest eg: edgeR.age.tmt
#   comp    : Name of comparison
#   comp.cols : Columns that contain the comparisons of interest. eg. young.cool.vs.ctrl = 1 and aged.cool.vs.ctrl = 4 so comp.cols = c(1,4)
#----------------------------------------------------------------------------------------------
runComps <- function(x.norm,dt,comp,comp.cols){
  
  # Will contain list of genes
  res = NULL
  
  # Will contain name of comparison -1 = Down; 0 = Neutral, 1 = Up
  names = NULL
  
  for(g in c(-1,0,1)){
    for(h in c(-1,0,1)){
      print(paste(g,h))
      temp = names(which(edgeR.age.tmt[,comp.cols[1]] == g & edgeR.age.tmt[,comp.cols[2]]==h))
      name = paste(g,h,sep="_")
      names = c(names,name)
      res = c(res,list(temp))
    }
  }
  
  # Name list members
  names(res) = names
  
  # Print length of each member list
  print(lapply(res,length))
  
  # Perform enrichment for each level of overlap
  enrich = prep.goseq(res,dge.aged.norm,mmu.kegg)
  
  # Write the output to file
  write.table(enrich[[1]],paste(res.dir.aged,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste(comp,"young.vs.aged-edgeR-based",pval,logfc,"goseq-results.txt",sep="_"),sep="_"),sep="/"),sep="\t",quote=F,row.names=F)
  write.table(enrich[[2]],paste(res.dir.aged,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste(comp,"young.vs.aged-edgeR-based",pval,logfc,"kegg-results.txt",sep="_"),sep="_"),sep="/"),sep="\t",quote=F,row.names=F)
  
  return(res)
}