# Functions needed for OOPS

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

myProtMapper <- function(ids,id.type="uniprot",out.fields=c("interpro.short_desc","ensembl.gene","go.MF.id","go.CC.id","go.BP.id"),species=9606,modify=T){
  
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

enricherPlot<-function(data,suf,N=5,colorBy = "p.adjust",low.col="#E69F00", high.col="#999999",trunc.len=40,all.size=10,y.size=12,x.size=14){
  
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
  
  gcsize = sapply(df$GeneRatio,function(x) strsplit(x,"/")[[1]][2])
  df$size =  sapply(df$GeneRatio,function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))
  df$tot.size <- paste(as.character(df$Cluster),"\n", "(", gcsize, ")", sep="")
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

