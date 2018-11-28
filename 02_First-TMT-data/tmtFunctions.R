#----------------------------------------------------------------------------
# Function   : enrichK
# Aim        : Call clusterProfiler's enrichKEGG on multiple comparisons
# Input      : 
#     isect  : Object resulting from running function 'venn' on a data frame. Contains intersections as gene lists
#     ix     : Index which is of interest. This can either be a number eg : 280 or a name eg : `150mJ.1:150mJ.2:150mJ.3`
#     data   : Object that contains heavy and light abundance (or protein abundance) values for each protein. Shoudl contain columns "light.log","heavy.log" and     #
#               "norm.abundance.ratio". If your columns are names differently, change this code. 
#     pval   : pval cut-off for KEGG enrichment
#   res.dir  : directory in which you want the results written. 
# Output     : Text files with most DE genes for each comparison as well as volcano plots, heatmaps
#----------------------------------------------------------------------------

enrichK <- function(isect,ix,data,pval,res.dir){
  
  genes = isect[ix][[1]]
  comp.name = ""
  
  # Getting the name of the intersection
  if(is.numeric(ix)){
    comp.name = names(isect)[ix]
  }else{
    comp.name = ix
  }
  
  gene.data = data[which(data$accessions %in% genes),]
  plot(gene.data$light.log,gene.data$heavy.log)
  abline(0,1,col="red")
  hist(gene.data$norm.abundance.ratio)
  
  # Separate out the up and down regulated genes
  down.in.crosslink = unique(gene.data[which(gene.data$norm.abundance.ratio<=0),"accessions"])
  up.in.crosslink = unique(gene.data[which(gene.data$norm.abundance.ratio>0),"accessions"])
  
  # Do some KEGG enrichment analysis
  down.names <- bitr(down.in.crosslink, fromType="UNIPROT", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  kk.down <- as.data.frame(enrichKEGG(gene = down.names$ENTREZID, organism = 'hsa',pvalueCutoff = pval))
  kk.down$dir = "Down"
  kk.down$comparison = comp.name
  head(kk.down)
  
  up.names <- bitr(up.in.crosslink, fromType="UNIPROT", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  kk.up <- as.data.frame(enrichKEGG(gene = up.names$ENTREZID, organism = 'hsa',pvalueCutoff = pval))
  kk.up$dir = "Up"
  kk.up$comparison = comp.name
  head(kk.up)
  
  kk = rbind(kk.up, kk.down)
  
  for(i in 1:nrow(kk)){
    kk$protein[i] = paste(bitr(unlist(strsplit(kk$geneID[i],"/")),fromType="ENTREZID", toType=c("UNIPROT"), OrgDb="org.Hs.eg.db")$UNIPROT,collapse=";")
    kk$symbol[i] = paste(bitr(unlist(strsplit(kk$geneID[i],"/")),fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")$SYMBOL,collapse=";")
  }
  return(kk)
}
