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

runDE <- function(x.norm,des,contr.matrix,res.dir=out.dir,logfc=1,pval=0.05){
  
  print("Running limma analysis")
  
  # Fit a linear regression model to the data with Bayesian correction
  vfit <- lmFit(exprs(x.norm), des)
  vfit.contr <- contrasts.fit(vfit,contrasts=contr.matrix)
  efit <- eBayes(vfit.contr)
  
  # Look at the summary of differentially expressed genes based on adj.p.val cut-off
  edgefit = decideTests(efit)
  sum.fit = summary(edgefit)
  
  for(c in 1:ncol(contr.matrix)){
    
    # Name of the comparison
    cont.name = colnames(contr.matrix)[c]
    
    # Primary group of interest
    up.label = strsplit(cont.name,"vs")[[1]][1]
    
    # Samples in the comparison group
    cont.col.names = paste(names(which(contr.matrix[,c] != 0)),collapse="|")
    
    # Write out topDE genes
    lim.top = data.frame(topTable(efit,coef=c,n=Inf,p.value = Inf))
    write.table(lim.top,paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste(cont.name,pval,"limma-based-DE-test-results.txt",sep="_"),sep="_"),sep="/"),sep="\t",quote=F,row.names=T)
    
    # Plot DE genes
    temp = plotDE(efit,x.norm,"logFC",NULL,"P.Value","adj.P.Val",pval=0.05,logfc=1,paste("limma",cont.name,sep="_"),res.dir, cont.col.names)
  }
  sum.edge.file = paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste("Summary-of-results.txt",sep=""),sep="_"),sep="/")
  write.table(sum.fit,sum.edge.file,row.names=T,quote=F,sep="\t") 
  return(edgefit)
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
  
  # Run 'decideTests' for all contrasts
  dt.z = decideTests(results)
  
  # Get column number for contrast of interest
  p = grep(strsplit(suf,"_")[[1]][2],colnames(results))
  
  # Pick top genes with p-values <= pval
  top = as.data.frame(topTable(results,coef=p,n=Inf,p.value=pval))
  
  # Collect all genes for a given contrast
  all = as.data.frame(topTable(results,coef=p,n=Inf,p.value=Inf))
  
  #---------------------
  # Draw volcano plots
  #---------------------
  
  # Log10(P-value) equivalent of the highest adj.P.Value <=0.05
  log.edge = max(-log10(all[which(all[,fdr.col] == max(all[which(all[,fdr.col]<=pval),fdr.col])),pval.col]))
  print(paste("-log10(P-value) corresponding to adj.p.val of 0.05 is",round(log.edge,2),"for this sample.",sep=" ")) 
  
  # Drawing volcano plot
  with(all, plot(get(logfc.col), -log10(get(pval.col)) ,xlab="logFC",ylab="-log10(p.value)",pch=20, main=paste("Volcano plot : DE genes for ",suf,"comparison",sep=" ")))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  abline(h=log.edge,col="red3",lty=3)
  abline(v=logfc,col="purple",lty=1)
  abline(v=-logfc,col="purple",lty=1)
  with(subset(all, get(fdr.col)<=pval ), points(get(logfc.col), -log10(get(pval.col)), pch=20, col="blue3")) # Significant
  with(subset(all, abs(get(logfc.col))>=logfc), points(get(logfc.col), -log10(get(pval.col)), pch=20, col="orange")) # Large fold change
  with(subset(all, get(fdr.col)<=pval & abs(get(logfc.col))>logfc), points(get(logfc.col), -log10(get(pval.col)), pch=20, col="green3")) # Both

  legend("topright",c("Adj.P.Value <= 0.05",paste("abs(logFC) >",logfc),paste("Adj.P.Value<=0.05 & abs(logFC) >",logfc),"Not significant"),fill=c("blue3","orange","green3","black"))
  
  # Label a subset of points with the textxy function from the calibrate plot
  with(subset(top,abs(get(logfc.col))>=logfc), textxy(get(logfc.col), -log10(get(pval.col)), labs=rownames(top), cex=.8,offset=0.6))
  
  # Add text for categorisation
  up.label = strsplit(suf,"vs")[[1]][1]
  text(-5,-2,paste("Up in",up.label,"mice",sep=" "),col = "blue3")
  text(5,-2,paste("Down in",up.label,"mice",sep=" "),col = "orange2")
  #text((max(all[,logfc.col])-0)/2,max(-log10(all[,pval.col])),paste("Up in",up.label,"mice",sep=" "),col = "blue3")
  #text((0+min(all[,logfc.col])/2),max(-log10(all[,pval.col])),paste("Down in",up.label,"mice",sep=" "),col = "orange2")
  
  # Extracting the top (adj.p.val < 0.05) most DE genes for all comparisons
  topgenes <- rownames(all)[which(all[,fdr.col]<=pval & abs(all[,logfc.col]) >= logfc)]
  k <- which(rownames(results) %in% topgenes)
  
  # If number of genes 'k' is < 2, do not proceed any further.
  if(length(k) > 2){
    
    # Set up colour vector for celltype variable
    mypalette <- brewer.pal(11,"RdBu")
    morecols <- colorRampPalette(rev(mypalette))
    
    # Subset data for heatmaps from counts matrix
    hmap.dat = as.matrix(x.norm$counts[k,])
    rownames(hmap.dat) = rownames(x.norm$counts)[k]
    
    # Include only those samples that are part of this comparison.
    colnums = grep(cont.col.names,colnames(hmap.dat))
    hmap.dat.sub = hmap.dat[,colnums]
    
    # Add a column colour option to distinguish the two comparison groups
    col.cols = palette()[as.numeric(as.factor(sapply(strsplit(as.character(colnames(hmap.dat.sub)),"\\.[0-9]+"),"[[",1)))] # Column side colours
    
    # Finally, draw the heatmap
    heatmap.2(hmap.dat.sub,scale='row',Rowv=F,Colv=F, ColSideColors = col.cols, col=rev(morecols(50)),labRow=rownames(results)[k],labCol=rownames(x.norm$target),cexCol=0.7,margin=c(8,6), lhei=c(2,10), dendrogram="none",trace="none",key=T, density.info="none",main=paste("Top DE genes across",suf,sep=" "))
  }else{
    print("Not enough data for heatmap")
  }
  dev.off()
}
