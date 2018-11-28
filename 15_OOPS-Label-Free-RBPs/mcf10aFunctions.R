#--------------------------------------------------------------------------------------------
# Function  : vargenes
# Aim       : Looking for highly/lowly variable genes in imputed, normalised, batch corected data
# Input     :
#           :  qntp = Normalised,Batch corrected MsnSet object
#           :  suf = Prefix for the output file name to say what data it is
#           :  reps = a vector of replicates in your experiment
#           :  idtype = id type that denotes each protein. Eg: "uniprot.Swiss.Prot","symbol"
#           :  labely = gene names to be used for the plot. Default = rownames(qntp)
#--------------------------------------------------------------------------------------------

vargenes <- function(qntp, suf, reps, idtype){
  
  pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),suf,"Finding-most-and-least-variable-genes-in-dataset.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  
  # Calculate means/variance
  exp = NULL
  if(class(qntp) == "MSnSet")
    exp = exprs(qntp)
  else
    exp = qntp
  
  # Calculating means and variance
  vars <- apply(exp,1,var,na.rm=T)
  means <- rowMeans(exp,na.rm=T)
  cv2 <- vars/means^2

  # Find genes with highest variance and lowest variance
  high.var <- names(sort(vars, decreasing=TRUE))[1:50]
  head(high.var)
  low.var <- names(sort(vars))[1:50]
  head(low.var)
  
  
  # Plot smoothed density plot
  labely = gsub("_HUMAN","",sapply(strsplit(high.var,"\\."),"[[",1))
  labelx = gsub("_HUMAN","",sapply(strsplit(low.var,"\\."),"[[",1))
  
  smoothScatter(log2(means),log2(vars),main="Colour density representation of mean and variance",xlab = "log2(Mean relative expression)",ylab = "log2(variance)")
  points(log2(means[high.var]),log2(vars[high.var]),col="red")
  text(log2(means[high.var[1:10]]),log2(vars[high.var[1:10]]),labels=labely[1:10])
  points(log2(means[low.var]),log2(vars[low.var]),col="purple")
  text(log2(means[low.var[1:10]]),log2(vars[low.var[1:10]]),labels=labelx[1:10])
  legend("topright",legend = c("High variability","Low variability"),fill = c("red","purple"))

  #------------------------------
  # Plots for most variable genes
  #------------------------------
  cols = c("#edf8b1","#7fcdbb","#2c7fb8","#c994c7")
  
  high.dat = as.matrix(exp[high.var,]/apply(exp[high.var,],1,max))
  high.dat[is.na(high.dat)] = 0
  heatmap(high.dat,zlim=c(0,1),col=grey.colors(20),scale="none",ColSideColors=cols[reps],labRow = labely[as.numeric(high.var)],main = "Heatmap of 50 most variable proteins")
  
  # Melt m and draw line graphs by gene
  mmod = as.data.frame(cbind(Gene=labely,exp[high.var,]))
  mmelt = melt(mmod,id = "Gene")
  colnames(mmelt) = c("Gene","Condition","Exp")
  mmelt$Exp = as.numeric(mmelt$Exp)
  mg = ggplot(data=mmelt,aes(x=Condition, y=Exp, group = Gene,color = Gene)) + geom_line() + theme(axis.text.x=element_text(angle=90, hjust=1),legend.position = "none",plot.title = element_text(hjust = 0.5))+facet_wrap(~Gene,scales="free_y")+labs(title="Expression of 50 MOST variable proteins across conditions")
  print(mg)
  
  # Annotation of most variable genes
  m2 = data.frame(cSplit(mmod,"Gene","; ", direction="long"))
  tmp = data.frame(queryMany(m2$Gene,scopes=idtype,fields=c("uniprot","go.BP.term","go.CC.term","go.MF.term","interpro.short_desc"),species=9606))
  m3 = merge(m2, tmp, by.x = "Gene", by.y = "query", all.x=T)
  
  for(k in c("go.BP","go.CC","go.MF","interpro")){
    name = paste(k,"mod",sep=".")
    m3[,name] = sapply(m3[,k],function(x) paste(unlist(x),collapse="; "))
  }

  #------------------------------
  # Plots for most least variable genes
  #------------------------------
  
  low.dat <- as.matrix(exp[low.var,]/apply(exp[low.var,],1,max))
  heatmap(low.dat,zlim=c(0,1),col=grey.colors(20),labRow=labelx,scale="none",ColSideColors=cols[reps],main = "Heatmap of 50 least variable proteins")
  
  # Melt m and draw line graphs by gene
  lowmod = as.data.frame(cbind(Gene=labelx,exp[low.var,]))
  lowmelt = melt(lowmod,id = "Gene")
  colnames(lowmelt) = c("Gene","Condition","Exp")
  lowmelt$Exp = as.numeric(lowmelt$Exp)
  ng = ggplot(data=lowmelt,aes(x=Condition, y=Exp, group = Gene,color = Gene)) + geom_line() + theme(axis.text.x=element_text(angle=90, hjust=1),legend.position="none",plot.title = element_text(hjust = 0.5))+facet_wrap(~Gene,scales="free_y")+labs(title="Expression of 50 LEAST variable proteins across conditions")
  print(ng)
  
  # Annotation of most variable genes
  n2 = data.frame(cSplit(lowmod,"Gene","; ", direction="long"))
  tnp = queryMany(n2$Gene,scopes=idtype,fields=c("symbol","go.BP.term","go.CC.term","go.MF.term","interpro.short_desc"),species=9606,returnall=T)
  n3 = merge(n2, tnp, by.x = "Gene", by.y = "query", all.x=T)
  for(k in c("go.BP","go.CC","go.MF","interpro")){
    name = paste(k,"mod",sep=".")
    n3[,name] = sapply(n3[,k],function(x) paste(unlist(x),collapse="; "))
  }

  dev.off()
  
  return(list(high=m3,low=n3))  
  
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
    #temp = plotDE(efit,x.norm,"logFC",NULL,"P.Value","adj.P.Val",pval=0.05,logfc=1,paste("limma",cont.name,sep="_"),res.dir, cont.col.names)
  }
  sum.edge.file = paste(res.dir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),paste("Summary-of-results.txt",sep=""),sep="_"),sep="/")
  write.table(sum.fit,sum.edge.file,row.names=T,quote=F,sep="\t") 
  return(list(edgefit,efit,lim.top))
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
  print(suf)
  p = grep(paste("^",suf,"$",sep=""),colnames(results))
  
  # Pick top genes with p-values <= pval
  top = as.data.frame(topTable(results,coef=p,n=Inf,p.value=pval))
  
  # Collect all genes for a given contrast
  all = as.data.frame(topTable(results,coef=p,n=Inf,p.value=Inf))
  
  #---------------------
  # Draw volcano plots
  #---------------------
  
  # If number of genes 'k' is < 2, do not proceed any further.
  if(nrow(top) > 2){
    
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
    
    legend("topleft",c("Adj.P.Value <= 0.05",paste("abs(logFC) >",logfc),paste("Adj.P.Value<=0.05 & abs(logFC) >",logfc),"Not significant"),fill=c("blue3","orange","green3","black"))
    
    # Label a subset of points with the textxy function from the calibrate plot
    with(subset(top,abs(get(logfc.col))>=logfc), textxy(get(logfc.col), -log10(get(pval.col)), labs=rownames(top), cex=.8,offset=0.6))
    
    # Add text for categorisation
    up.label = strsplit(suf,"vs")[[1]][1]
    text(10,5,"Up",col = "blue3")
    text(-10,5,"Down",col = "orange2")
    #text((max(all[,logfc.col])-0)/2,max(-log10(all[,pval.col])),paste("Up in",up.label,"mice",sep=" "),col = "blue3")
    #text((0+min(all[,logfc.col])/2),max(-log10(all[,pval.col])),paste("Down in",up.label,"mice",sep=" "),col = "orange2")
    
    # Extracting the top (adj.p.val < 0.05) most DE genes for all comparisons
    topgenes <- rownames(all)[which(all[,fdr.col]<=pval & abs(all[,logfc.col]) >= logfc)]
    k <- which(rownames(results) %in% topgenes)
    
    # Set up colour vector for celltype variable
    mypalette <- brewer.pal(11,"RdBu")
    morecols <- colorRampPalette(rev(mypalette))
    
    # Subset data for heatmaps from counts matrix
    hmap.dat = as.matrix(exprs(x.norm)[k,])
    rownames(hmap.dat) = rownames(exprs(x.norm))[k]
    
    # Include only those samples that are part of this comparison.
    colnums = grep(paste(toupper(cont.col.names),sep=""),colnames(hmap.dat))
    hmap.dat.sub = hmap.dat[,colnums]
    
    # Add a column colour option to distinguish the two comparison groups
    col.cols = palette()[as.numeric(as.factor(sapply(strsplit(as.character(colnames(hmap.dat.sub)),"\\_[A-Z]+"),"[[",1)))] # Column side colours
    
    # Finally, draw the heatmap
    heatmap.2(hmap.dat.sub,scale='row',Rowv=F,Colv=F, ColSideColors = col.cols, col=rev(morecols(50)),labRow=rownames(results)[k],labCol=rownames(x.norm$target),cexCol=1,srtCol=45,margin=c(8,6), lhei=c(2,10), dendrogram="none",trace="none",key=T, density.info="none",main=paste("Top DE genes across",suf,sep=" "))
  }else{
    print("Not enough data for heatmap")
  }
  dev.off()
}

#----------------------------------------------------------------------------------------------
# Function: enricherPlot
# Aim : Modify DOSE::plot to use colours I like for plotting results of compareCluster
# Default : Will only plot the dot size to show GeneRatio and colour to show adjusted.p.val in grey and gold
# Input : 
#   data : object from compareCluster function
#   N    : Number of top terms to show on the plot Default = 5
#   colorBy : What numeric value do you want the colour scale based on ? Default = p.adjust
#   low.col : What colour would you like your low 'colorBy' values to be ? Default = grey
#   high.col : What colour would you like your high 'colorBy' values to be ? Default = gold
#   trunc.len : At what length do you want your GO/Interpro/KEGG terms truncated ? Default = 40
#   suf : Suffix for output file
#   all.size : What is the size that you want your legend and label text to be ? Default = 10
#   y.size : What is the size that you want for your y-axis labels ?
#   x.size : What is the size that you want for your x-axis labels ?
#----------------------------------------------------------------------------------------------

enricherPlot<-function(data,suf,N=5,colorBy = "p.adjust",sizeBy="obsRatio",low.col="#E69F00", high.col="#999999",trunc.len=40,all.size=10,y.size=12,x.size=14){
  
  #--------------------------------------------------------------------------
  # Function : topN 
  # Aim : Picks the top "N" terms in an enrichment analysis for each cluster
  #--------------------------------------------------------------------------
  
  topN <- function(res, showCategory){
    ddply(.data = res,
          .variables = .(Cluster),
          .fun = function(df, N) {
            if (length(df$Count) > N) {
              if (any(colnames(df) == "pvalue")) {
                idx <- order(df$pvalue, decreasing=FALSE)[1:N]
              } else {
                ## for groupGO
                idx <- order(df$Count, decreasing=T)[1:N]
              }
              return(df[idx,])
            } else {
              return(df)
            }
          },
          N=showCategory)
  }
  
  # Convert 'compareCluster' result to a data.frame
  df = data.frame(data)
  
  # 'gcsize' is the number of proteins in each dataset that could be mapped to GO/Interpro/KEGG. It is the denominator in 'GeneRatio'
  # 'size' = GeneRatio is a text field - split its elements and calculate the actual GeneRatio or proportion of genes contributing to term enrichment
  # 'tot.size' = Modified x-axis labels to contain count for each cluster
  # 'mod.desc' = Modify the length of the description of terms to be 40 characters long. Anything longer will be truncated and followed by "..."
  
  df$size =  sapply(df$obsRatio,function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))
  df$tot.size <- paste(as.character(df$Cluster),"\n", "(", df$Foreground, ")", sep="")
  df$mod.desc = as.character(df$Description)
  df$mod.desc[which(nchar(df$mod.desc)>trunc.len)] = paste(substring(df$mod.desc[which(nchar(df$mod.desc)>trunc.len)],1,trunc.len),"...",sep="")
  
  # Once you've modified the main data frame, subset it to only include the top 'N' terms for each cluster
  # Order this data.frame such that the most enriched terms are at the top of the figure
  df.sub.org = topN(df,N)
  df.sub = df[which(df$mod.desc %in% unique(df.sub.org$mod.desc)),]
  
  idx <- order(df.sub$size, decreasing = F)
  df.sub$mod.desc <- factor(df.sub$mod.desc, levels=unique(df.sub$mod.desc[idx]))
  
  # Draw the plot
  #pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),suf,N,"enricher-dotplot.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  gp = ggplot(df.sub, aes_string(x="tot.size", y="mod.desc", size=df.sub$size, color=colorBy)) + geom_point() + scale_size(range=c(3,8))+ scale_color_gradient(low=low.col, high=high.col) +xlab("")+ylab("")+guides(size=guide_legend("GeneRatio",order=1),color=guide_colorbar(label.theme = element_text(angle = -45)))+theme_bw()+theme(text = element_text(size=all.size),axis.text.x=element_text(size=x.size),axis.text.y=element_text(size=y.size),legend.direction = "horizontal", legend.position = "top",legend.box = "vertical")
  #print(gp)
  #dev.off()
  
  return(gp)
}

# ------------------------------------------------------------------------------------------------------------
# Function	: myProtMapper 
# Aim		: To use the function 'queryMany' from Bioconductor package mygene as fast and most up-to-date
# Input 
# 	: ids = a character list of ids which can be uniprot, ensembl gene, gene symbol,etc
#       : id.type = what type of ids have you provided in the 'ids' list. Default = "uniprot"
#       : outlist = list of ids you want as an output. Default = c("interpro","ensembl.gene","go")
#       : modify = Logical, Default = T; Would you like to modify fields such as interpro, enseml, go to make them more human readable.
# Output  	: A dataframe with required ids and input ids 
# ------------------------------------------------------------------------------------------------------------------------

myProtMapper <- function(ids,id.type="uniprot",out.fields=c("interpro.short_desc","ensembl.gene","go.MF.id","go.CC.id","go.BP.id","kegg"),species=9606,modify=T){
  
  # Get the mapping
  qm = queryMany(ids,scopes=id.type,fields=out.fields,species=9606)
  
  # Returning variable
  ret.qm = NULL
  
  # Resolve the mappings to make them human readable
  if(modify == T){
    qm$go.all = NULL
    
    # Interpro mappings
    if(is.element("interpro",colnames(qm))){
      qm$domains = sapply(qm$interpro,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No Interpro domains")
    }
    
    # GO mappings
    if(!is.na(grep("go",colnames(qm)))){
      
      # Grep all the go columns 'go.CC','go.MF','go.BP'
      f = grep("go",colnames(qm), value=T)
      
      qm$go.all = apply(qm[,f], MARGIN=1, FUN = function(x) paste0(as.character(unique(unlist(x))), collapse=";"))
      qm$go.all = gsub("^;","",gsub(";;",";",qm$go.all))
      qm$go.count = lengths(strsplit(qm$go.all,";"))
    }
    else{
      print("No GO terms")
    }
    
    # Ensembl.gene mappings
    if(is.element("ensembl",colnames(qm))){
      qm$ens = sapply(qm$ensembl,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No Ensembl Ids")
    }
    
    # Return mapped structure with tidy columns
    ret.qm = qm
  }
  else{
    ret.qm = qm
  }
  
  return(data.frame(ret.qm))
}

# ------------------------------------------------------------------------------------------------------------
# Function  : 'makeGene2Cat' to produce a 1:1 mapping of uniprot/ensembl/symbols to GO/Interpro terms. 
#              Will be used as input into the 'goseq' function in the gene2cat slot
# Input 
#           : dat = dataframe with ids and go/interpro terms (obtained from myProtMapper)
#           : from.id =  ids you want to map 'from'. Default = "uniprot"
#           : to.id =  ids you want to map to c("interpro","ensembl.gene","go")
#           : splt = symbol you want to split by if there are multiple ids
# Output  : A two column dataframe with Uniprot ids in the first and Go/Interpro in the second
# ------------------------------------------------------------------------------------------------------------------------

makeGene2Cat <- function(dat,from.id,to.id,splt){
  
  cat.frame = dat[,c(from.id,to.id)]
  d.dt = data.table(cat.frame,key=colnames(cat.frame[,from.id]))
  cat.out = data.frame(d.dt[, list(to.id = unlist(strsplit(get(to.id), splt))), by=from.id])
  cat.out = unique(cat.out)
  
  return(cat.out)
}

#--------------------------------------------------------------------------------------------
# Function: rungoseq
# Aim:  goseq analysis
# Input : genes = genelist of interest; g2c = gene to category mapping i.e univ.go or univ.pro
# b = bias data ; bh = Bonferroni p-value cutoff. 
#--------------------------------------------------------------------------------------------

rungoseq<-function(genes,g2c, univ, b, bh){
  
  all.genes = rep(0,length(univ))
  names(all.genes) = univ
  all.genes[which(names(all.genes) %in% genes)] = 1
  t = table(all.genes)
  
  if(is.null(b)){
    pwf = nullp(all.genes,bias.data = NULL,genome="hg19",id="ensGene",plot.fit = T)
  }else{
    pwf = nullp(all.genes,bias.data = b,plot.fit = T)
  }
  
  GO.wall = goseq(pwf,gene2cat = g2c)
  GO.wall$BH = p.adjust(GO.wall$over_represented_pvalue,method = "BH")
  GO.wall$Foreground = length(intersect(genes,g2c$query))
  GO.wall$Background = length(unique(g2c$query))
  GO.wall$obsRatio = paste(GO.wall$numDEInCat,GO.wall$Foreground, sep="/")
  GO.wall$bgRatio = paste(GO.wall$numInCat,GO.wall$Background, sep="/")
  GO.wall$expectDE = ceiling(GO.wall$Foreground*(GO.wall$numInCat/GO.wall$Background))
  GO.wall$expRatio = paste(GO.wall$expectDE,GO.wall$Foreground, sep = "/")
  GO.wall$foldEnrich = round(GO.wall$numDEInCat/GO.wall$expectDE,2)
  
  GO.enriched = GO.wall[which(GO.wall$BH <= bh),]
  GO.enriched$geneID = sapply(GO.enriched$category,function(x) paste(intersect(genes,g2c$query[which(g2c$to.id == x)]),collapse="/"))
  GO.enriched$Count = sapply(GO.enriched$category,function(x) length(intersect(genes,g2c$query[grep(x,g2c$to.id)])))
  
  GO.enriched$neg.log10.BH = -log10(GO.enriched$BH)
  GO.enriched$neg.log10.BH[which(GO.enriched$BH == 0)] = max(GO.enriched$neg.log10.BH[is.finite(GO.enriched$neg.log10.BH)])+1
  
  return(list(GO.wall,GO.enriched))
}

#--------------------
# reverseMapping
#--------------------
reversemapping=function(map){
  tmp=unlist(map,use.names=FALSE)
  names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
  return(split(names(tmp),as.vector(tmp)))
}

#----------------------------------------------------------------------------------------------
# Function: enricherPlot
# Aim : Modify DOSE::plot to use colours I like for plotting results of compareCluster
# Default : Will only plot the dot size to show GeneRatio and colour to show adjusted.p.val in grey and gold
# Input : 
#   data    : object from compareCluster function
#   N       : Number of top terms to show on the plot Default = 5
#   colorBy : What numeric value do you want the colour scale based on ? 
#             Default = BH or Benjamini-Hochberg adjusted p-value
#   sizeBy  : What numeric value do you want the size of the dots based on ? Default = obsRatio
#   low.col : What colour would you like your low 'colorBy' values to be ? Default = grey
#   high.col: What colour would you like your high 'colorBy' values to be ? Default = gold
#   trunc.len: At what length do you want your GO/Interpro/KEGG terms truncated ? Default = 40
#   suf     : Suffix for output file
#   all.size: What is the size that you want your legend and label text to be ? Default = 10
#   y.size  : What is the size that you want for your y-axis labels ?
#   x.size  : What is the size that you want for your x-axis labels ?
#----------------------------------------------------------------------------------------------

enricherPlot<-function(data,suf,N=5,colorBy = "BH",sizeBy = "obsRatio",low.col="#E69F00", high.col="#999999",trunc.len=40,all.size=10,y.size=12,x.size=14){
  
  #--------------------------------------------------------------------------
  # Function : topN 
  # Aim : Picks the top "N" terms in an enrichment analysis for each cluster
  #--------------------------------------------------------------------------
  
  topN <- function(res, showCategory){
    ddply(.data = res,
          .variables = .(Cluster),
          .fun = function(df, N) {
            if (length(df$obsRatio) > N) {
              if (any(colnames(df) == "neg.log10.BH")) {
                idx <- order(df$neg.log10.BH, decreasing=FALSE)[1:N]
              } else {
                ## for groupGO
                idx <- order(df$numDEInCat, decreasing=T)[1:N]
              }
              return(df[idx,])
            } else {
              return(df)
            }
          },
          N=showCategory)
  }
  
  # Convert 'compareCluster' result to a data.frame
  df = data.frame(data)
  df$Cluster = factor(df$Cluster,levels=c(paste("F",1:20,sep="")))
  
  # 'gcsize' is the number of proteins in each dataset that could be mapped to GO/Interpro/KEGG. It is the denominator in 'GeneRatio'
  # 'size' = GeneRatio is a text field - split its elements and calculate the actual GeneRatio or proportion of genes contributing to term enrichment
  # 'tot.size' = Modified x-axis labels to contain count for each cluster
  # 'mod.desc' = Modify the length of the description of terms to be 40 characters long. Anything longer will be truncated and followed by "..."
  
  gcsize = sapply(df$obsRatio,function(x) strsplit(x,"/")[[1]][2])
  #df$size =  sapply(df$obsRatio,function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))
  df$tot.size <- paste(as.character(df$Cluster),"\n", "(", df$Direction,gcsize, ")", sep="")
  df$mod.desc = as.character(df$Description)
  df$mod.desc[which(nchar(df$mod.desc)>trunc.len)] = paste(substring(df$mod.desc[which(nchar(df$mod.desc)>trunc.len)],1,trunc.len),"...",sep="")
  
  # Once you've modified the main data frame, subset it to only include the top 'N' terms for each cluster
  # Order this data.frame such that the most enriched terms are at the top of the figure
  df.sub.org = topN(df,N)
  df.sub = df[which(df$mod.desc %in% unique(df.sub.org$mod.desc)),]
  
  idx <- order(df.sub[,colorBy], decreasing = F)
  df.sub$mod.desc <- factor(df.sub$mod.desc, levels=unique(df.sub$mod.desc[idx]))
  df.sub$tot.size = factor(df.sub$tot.size,levels = unique(df.sub$tot.size))
  
  # Draw the plot
  pdf(paste(outdir,paste(suf,N,"enricher-dotplot.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  gp = ggplot(df.sub, aes_string(x="tot.size", y="mod.desc", size=sizeBy, color=colorBy)) + geom_point() + scale_size(breaks = c(0,2,4,8,16), limits=c(0,20))+ scale_color_gradient2(low=low.col,mid=high.col, high = "#56B4E9",midpoint = quantile(df.sub[,colorBy],0.98))+xlab("")+ylab("")+guides(size=guide_legend("Fold enrichment",order=1),color=guide_colorbar(title = "-log10(adj.p.value)", title.vjust = 1.0, label.position = "bottom", reverse = F,label.theme = element_text(angle = -90)))+theme_bw()+theme(text = element_text(size=all.size),axis.text.x=element_text(size=x.size),axis.text.y=element_text(size=y.size),legend.direction = "horizontal", legend.position = "top",legend.box = "vertical")
  print(gp)
  dev.off()
  
  return(gp)
}

