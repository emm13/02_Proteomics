# The rownames of samp.dat have to be the same as column names in the expression data matrix
filt.6 = cbind(apply(filt.5[,1:2],1,mean),apply(filt.5[,3:4],1,mean),apply(filt.5[,5:6],1,mean),apply(filt.5[,7:8],1,mean))
colnames(filt.6) = c("NCL","Unstraved","Starved","Insulin")

samp.dat = colnames(filt.6)
samp.dat = cbind(samp.dat,Tag=paste("TMT",1:4,sep=""))
colnames(samp.dat)[1] = "Sample"
samp.dat = data.frame(samp.dat,stringsAsFactors = F)
rownames(samp.dat) = samp.dat$Sample

# Create an MSnSet object
filt.6 = filt.6[order(rownames(filt.6)),]
rownames(comprots.ann) = comprots.ann$Uniprot.name
res <- MSnSet(exprs = as.matrix(filt.6),fData=comprots.ann[10:17],pData = samp.dat)
fData(res)$gene = sapply(strsplit(rownames(fData(res)),"_"),"[[",1)

plot2D(res, fcol = NULL, col = "black")
heatmap.2(exprs(res))
boxplot(exprs(res))

res.tmp = res[which(fData(res)$gene %in% names(km.res$cluster[which(km.res$cluster == 3)]))]
plotDist(res,pcol="lightblue")


# Clustering methods
res.dist <- get_dist(t(exprs(res)), stand = TRUE, method = "pearson")
fviz_dist(res.dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

fviz_nbclust(exprs(res), kmeans, "silhouette")

set.seed(123)

tmp = filt.5
rownames(tmp) = sapply(strsplit(rownames(tmp),"_"),"[[",1)

km.res <- kmeans(exprs(res), 4, nstart = 25)

# Visualize
library("factoextra")
fviz_cluster(km.res, data = exprs(res),
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())


k = kmeans(exprs(res),4,iter.max = 20)
par(mfrow=c(2,2))
for(i in 1:4){
  n = names(k$cluster[which(k$cluster==i)])
  temp = res[which(rownames(exprs(res)) %in% n),]
  plotDist(temp,pcol=rainbow(4)[i])
}


res.hc <- exprs(res) %>%
  scale() %>%                    # Scale the data
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

# Visualize using factoextra
# Cut in 4 groups and color by groups
fviz_dend(res.hc, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)


set.seed(123)
# Compute
library("NbClust")
res.nbclust <- exprs(res)[,2:4] %>%
  scale() %>%
  NbClust(distance = "euclidean",
          min.nc = 2, max.nc = 10, 
          method = "complete", index ="all") 

res.hc <- exprs(res)[,2:4] %>%
  scale() %>%
  eclust("hclust", k = 4, graph = FALSE)

fviz_silhouette(res.hc)


# Working out right number of clusters
pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Vero-RIC_Silhouette-kmeans-based.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
par(mfrow=c(2,4))
col=rainbow(8)

km2 = kmeans(exprs(res),2)

km.res <- kmeans(exprs(res)[,2:4], 4, nstart = 25)
fviz_cluster(km.res, exprs(res)[,2:4],  geom = "point", 
             ellipse= FALSE, show.clust.cent = FALSE,
             palette = "jco", ggtheme = theme_classic())

# For remaining cluster options
for (i in 3:10){
  # Calculate k-means for a cluster value 'i'
  km = kmeans(exprs(res),i)                    
  
  # Calculate dissimilarity matrix for samples in prot.nmf
  dissE = daisy(exprs(res))
  
  # Calculate silhouette scores to judge how well the clustering defines the data 
  sk = silhouette(km$cl,dissE)
  
  # Draw silhouette plot
  plot(sk,col=col[i-2],main = paste("Silhouette plot for k = ",i,sep=""))
  
}

dev.off()

#---------------------------------------------------------------------------------------------------------------
# Comparing this list to Mark/David's list to see if this can be an orthogonal validation
#---------------------------------------------------------------------------------------------------------------
md = read.delim("Input/MarkDavid-RBPlist.txt",sep="\t",header=T)

md.vero = md[which(md$Uniprot.ID %in% unique(common.prots$Accession)),]
table(md.vero$RNAPI)

# Merge md, comm.melt and groupings data
gp = read.table("Input/Common-110-by-group.txt",sep="\t",header=T)
gp.md = merge(gp,md,by.x="Protein",by.y = "Uniprot.ID")
gp.md.up = merge(gp.md,upcom[,1:2],by.x="Protein",by.y="Entry.name")

# Check for any difference in groupings
# There is a larger number of Group1 proteins in RNAPI.dep vs RNAPI.indep samples
y = as.data.frame.matrix(t(table(gp.md$RNAPI,gp.md$Group)))
y$Group = rownames(y)
y$prop.test.pval = round(sapply(1:nrow(y), function(z) prop.test(c(y[z,1],y[z,2]),c(37,64),correct=F, alternative="two.sided")$p.value),4)
y$total.dep = 37
y$dep.prop = round(y$dependent/y$total.dep,2)
y$total.indep = 64
y$indep.prop = round(y$independent/y$total.indep,2)
y = y[,c(3,1,5,6,2,7,8,4)]

c = compareCluster(Entry~Group,data = gp.md.up,fun='enricher', universe = fData(qnt.prot.no.imp)$Master.Protein.Accessions,TERM2GENE = univ.cat.go[,c(2,1)],qvalueCutoff = 0.05)
c@compareClusterResult$Description = goterms[as.data.frame(c)$ID]
c@compareClusterResult$Ontology = ont[as.data.frame(c)$ID]

d = compareCluster(Entry~Group,data = gp.md.up,fun='enricher', universe = fData(qnt.prot.no.imp)$Master.Protein.Accessions,TERM2GENE = univ.cat.doms[,c(2,1)],qvalueCutoff = 0.05)
d.df = data.frame(d)
d.df$symbol = sapply(d.df$geneID, function(x) paste(bitr(unlist(strsplit(x,"/")),fromType="UNIPROT", toType=c("SYMBOL"), OrgDb=org.Hs.eg.db)$SYMBOL,collapse="/"))


ego = enricherPlot(data=c,N=5,suf="GO-Across-7-RBP-profile-groups",trunc.len=50,y.size=14,all.size=16)
epro = enricherPlot(data=d,N=5,suf="GO-Across-7-RBP-profile-groups",trunc.len=50,y.size=14,all.size=16)

c.df = data.frame(c)
c.df$symbol = sapply(c.df$geneID, function(x) paste(bitr(unlist(strsplit(x,"/")),fromType="UNIPROT", toType=c("SYMBOL"), OrgDb=org.Hs.eg.db)$SYMBOL,collapse="/"))


# Split into indep and dependent
md.vero.dep = as.character(md.vero$Uniprot.ID[which(md.vero$RNAPI == "dependent")])
md.vero.indep = as.character(md.vero$Uniprot.ID[which(md.vero$RNAPI == "independent")])

# Get the protein spectral counts
md.vero.dep.rbp = com.melt[which(com.melt$Accession %in% md.vero.dep),]
md.vero.indep.rbp = com.melt[which(com.melt$Accession %in% md.vero.indep),]

# Plot
pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"RNAPI-dependent-common-genes-across-callers-violinplot.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
gg2 = ggplot(data=md.vero.dep.rbp,aes(x=Tmt, y=Exp, group = Tmt,fill=Tmt)) + geom_violin() +stat_summary(aes(group=1),fun.y=median, geom="line", color="black", size=0.5,lty = 2) + theme(legend.position="top")+labs(fill = "Treatment")+ scale_fill_manual(values=c("#E2D200","#46ACC8","#E58601","#B40F20"))+facet_wrap(~Accession)
gg10 <- facet_multiple(plot=gg2, facets="Accession", ncol = 4, nrow = 4, scales = "free_y")
dev.off()

pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"RNAPI-independent-common-genes-across-callers-violinplot.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
gg2 = ggplot(data=md.vero.indep.rbp,aes(x=Tmt, y=Exp, group = Tmt,fill=Tmt)) + geom_violin() +stat_summary(aes(group=1),fun.y=median, geom="line", color="black", size=0.5,lty = 2) + theme(legend.position="top")+labs(fill = "Treatment")+ scale_fill_manual(values=c("#E2D200","#46ACC8","#E58601","#B40F20"))+facet_wrap(~Accession)
gg10 <- facet_multiple(plot=gg2, facets="Accession", ncol = 4, nrow = 4, scales = "free_y")
dev.off()

```
