#---------------------------------------------------------------------------
# Author 	: Manasa Ramakrishna, mr325@le.ac.uk
# Date started 	: 7th March, 2017
# Last modified : 28th March, 2017
# Aim 		: To analyse proteomics data using NMF 
# Notes   : To use pRolocdata on ALICE, need to 'module load netcdf/gcc/4.4.1' and 'module load R' 
#---------------------------------------------------------------------------

# Invoking libraries
library(ggplot2)
library(goseq)
library(GOstats)
library(limma)
library(NMF)
library(pRoloc)
library(pRolocdata)
library(ggbiplot)
require(gridExtra)


#Setting working directories
wd = "/home/m/mr325/Manasa/01_Code/R/01_nmf-for-dmel"
setwd(wd)
getwd()

indir = paste(wd,"Input",sep="/")
outdir = paste(wd,paste(Sys.Date(),"Output",sep = "_"),sep = "/")

if (exists(outdir)){
	print("Outdir exists")
}else{
	dir.create(outdir)
}

# Sourcing functions and files
source(paste(wd,"nmfFunctions.R",sep="/"))
infiles = list.files(indir,full.names=T)
print(infiles)

#===================================
# Processing proteomics data
#===================================

# Reading protein mass spec data for fly embryos
# The data used in this tutorial has been published in Denise J Tan, Heidi Dvinge, Andy Christoforou, Paul Bertone, Alfonso A Martinez, and Kathryn S
# Lilley. Mapping organelle proteins and protein complexes in drosophila melanogaster. J Proteome Res, 8(6):2667–78, Jun 2009. doi:10.1021/pr800866n. 
# The LOPIT technique [1] is used to localise integral and associated membrane proteins in Drosophila melanogaster embryos.The original localisation analysis was 
# performed using partial least square discriminant analysis (PLS-DA). Relative quantitation data was retrieved from the supplementary file pr800866n si 004.xls3
# and replicate 1 is used below.

prot <- read.delim(grep("*.csv",infiles,value=T),sep=",",header=T)
frac <- read.delim(grep("fractions.txt",infiles,value=T),sep="\t",header=T)

prot.nmf = prot[,c(7:10)]
rownames(prot.nmf) = prot$FBgn
prot.meta = prot[,c(1:6,11)]

# So prot.meta data has identifiers in the UniProt format. However, need to map them to FLYBASE markers
# To do this, I exported to file, uploaded to UniProt's Retrieve ID/mapping tool (http://www.uniprot.org/uploadlists/)
# Downloaded Flybase identifier mapping to UniProt, manually curated markets.

prolocMark = as.data.frame(pRolocmarkers("dmel"))
prolocMark = cbind(rownames(prolocMark),prolocMark)
colnames(prolocMark) = c("UniProt","Location")

# Reading in FLybase-Uniprot mappings
prolocFly = read.table("Input/ProlocData-manually-curated-markers-with-Flybase.txt",sep="\t",header=T)
colnames(prolocFly) = c("UniProt","FBgn")
allLabels = merge(prolocFly,prolocMark,all.x=T,all.y=T)

# Merging proloc metadata with marker information (only available for 157/888 markers)
prot.meta.mark = merge(prot.meta,allLabels, all.x=T,all.y=F)
prot.meta.mark$PLS.DA.classification = gsub("^$","Unknown",prot.meta.mark$PLS.DA.classification)

#==============================================================================================
# Clustering by k-means
# Would pick the 'k' with highest silhouette score and sensible number of samples per group
#==============================================================================================

pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"DMel_Silhouette-kmeans-based.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
par(mfrow=c(2,4))
col=rainbow(8)

km2 = kmeans(prot.nmf,2)
cl.km = data.frame(cbind(FBgn=names(km2$cluster),"k.2"=km2$cluster))
m2 = merge(prot.meta,cl.km) # contains PLSDA and new Locations


# For remaining cluster options
for (i in 3:10){
  
  # Calculate k-means for a cluster value 'i'
	km = kmeans(prot.nmf,i)                    
	
	# Convert result to a data frame with sample and cluster assignment
	cl.km = data.frame(cbind(FBgn=names(km$cluster),km$cluster))
	
	# Add column name "k.i" eg: k.3
	colnames(cl.km)[2] = paste("k",i,sep=".")    
	
	# Merge with existing cluster data matrix
	m2 = merge(m2,cl.km)
	
	# Calculate dissimilarity matrix for samples in prot.nmf
	dissE = daisy(prot.nmf)
	
	# Calculate silhouette scores to judge how well the clustering defines the data 
	sk = silhouette(km$cl,dissE)
	
	# Draw silhouette plot
	plot(sk,col=col[i-2],main = paste("Silhouette plot for k = ",i,sep=""))
	
}

dev.off()

# Optimal clusters are either 3 or 7.
# Would like to work out what is in each cluster in the optimal cluster list
# Need to map proteins to gene ontology terms. 

# First, cluster size of 3
km3 = kmeans(prot.nmf,3)
#cl.km3 = data.frame(cbind(FBgn=names(km3$cluster),"k.3"=km3$cluster))
#m.km3 = merge(prot.meta.mark,cl.km3)

# Mapping to GO terms
dmel.ensembl = useEnsembl(biomart="ensembl", dataset="dmelanogaster_gene_ensembl")
listFilters(dmel.ensembl)[60:75,]
#                      name                                        description
# 60   flybase_annotation_id        Flybase Annotation ID(s) [e.g. FBan0011023]
# 61         flybase_gene_id              Flybase Gene ID(s) [e.g. FBgn0031208]
# 62   flybase_transcript_id        Flybase Transcript ID(s) [e.g. FBtr0300689]
#.....
# 72                   go_id             GO Term Accession(s) [e.g. GO:0005515]
# 73           go_to_gene_id                 Go to gene ID(s) [e.g. GO:0016020]
# 74    goslim_goa_accession         GOSlim GOA Accessions(s) [e.g. GO:0005623]
.....

# All proteins in list
fbgn.to.go = getBM(attributes=c('flybase_gene_id','go_id'),filters='flybase_gene_id',values=m.km3$FBgn,mart=dmel.ensembl)

# Proteins in cluster 3
fbgn.to.go.k3.cl1 = getBM(attributes=c('flybase_gene_id','go_id'),filters='flybase_gene_id',values=m.km3$FBgn[which(m.km3$k.3==1)],mart=dmel.ensembl)
fbgn.to.go.k3.cl2 = getBM(attributes=c('flybase_gene_id','go_id'),filters='flybase_gene_id',values=m.km3$FBgn[which(m.km3$k.3==2)],mart=dmel.ensembl)
fbgn.to.go.k3.cl3 = getBM(attributes=c('flybase_gene_id','go_id'),filters='flybase_gene_id',values=m.km3$FBgn[which(m.km3$k.3==3)],mart=dmel.ensembl)

write.table(cl.km3[which(cl.km3$k.3==1),],"2017-03-07_Output/Flybase-k3-cluster1.txt",sep="\t",quote=F,row.names=F)
write.table(cl.km3[which(cl.km3$k.3==2),],"2017-03-07_Output/Flybase-k3-cluster2.txt",sep="\t",quote=F,row.names=F)
write.table(cl.km3[which(cl.km3$k.3==3),],"2017-03-07_Output/Flybase-k3-cluster3.txt",sep="\t",quote=F,row.names=F)


#-----------------------------
# PCA visualisation
#-----------------------------

# PCA for overall data
pca.all = drawPCA(prot.nmf,prot.meta.mark,"All")

#Subsetting first cluster from kmeans clustering with k=3
prot.nmf.k3.cl1 = prot.nmf[!is.na(match(rownames(prot.nmf),names(km3$cluster[which(km3$cluster==1)]))),]
pca.k3.cl1 = drawPCA(prot.nmf.k3.cl1,prot.meta.mark,"k3-cluster1")

#Subsetting second cluster from kmeans clustering with k=3
prot.nmf.k3.cl2 = prot.nmf[!is.na(match(rownames(prot.nmf),m.km3$FBgn[which(m.km3$k.3==2)])),]
pca.k3.cl2 = drawPCA(prot.nmf.k3.cl2,prot.meta.mark,"k3-cluster2")

#Subsetting third cluster from kmeans clustering with k=3
prot.nmf.k3.cl3 = prot.nmf[!is.na(match(rownames(prot.nmf),m.km3$FBgn[which(m.km3$k.3==3)])),]
pca.k3.cl3 = drawPCA(prot.nmf.k3.cl3,prot.meta.mark,"k3-cluster3")

#Subsetting clusters 2&3 from kmeans clustering with k=3
prot.nmf.k3.cl2.3 = prot.nmf[!is.na(match(rownames(prot.nmf),m.km3$FBgn[which(m.km3$k.3!=1)])),]
pca.k3.cl2.3 = drawPCA(prot.nmf.k3.cl2.3,prot.meta.mark,"k3-clusters2-and-3")

# PCA on all original "Unknown"
prot.nmf.unknown = prot.nmf[which(rownames(prot.nmf) %in% prot.meta.mark$FBgn[which(prot.meta.mark$PLS.DA.classification == "Unknown")]),]
pca.unk = drawPCA(prot.nmf.unknown,prot.meta.mark,"Unknown-in-PLSDA")


#================
# Functions
#================

#---------------------------------------------------------------------------------------
# drawPCA: Function to run a Pricipal Components Analysis(PCA) and plot the results
#---------------------------------------------------------------------------------------
drawPCA <- function(data,meta,suf){

	# Scaled and centred
	dmel.pca = prcomp(data,center=T,scale=T)

	# Open PDF
	pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"DMel_PCA",suf,"with-known-markers.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
	
	# Plot based on old PLSDA classification	
	g <- ggbiplot(dmel.pca, obs.scale = 1, var.scale = 1, var.axes=F, groups = meta$PLS.DA.classification[match(rownames(data), meta$FBgn)], ellipse = F, circle = F)
	g <- g + scale_color_discrete(name = '')
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	
	# Plot based on improved location classification
	f <- ggbiplot(dmel.pca, obs.scale = 1, var.scale = 1, var.axes=F, groups = meta$Location[match(rownames(data), meta$FBgn)], ellipse = F, circle = F)
	f <- f + scale_color_discrete(name = '')
	f <- f + theme(legend.direction = 'horizontal', legend.position = 'top')	

	# Print side by side plots
	grid.arrange(g,f,ncol=2)
	
	# Close PDF
	dev.off()
	
	# Return PCA solution
	return(dmel.pca)
}


#---------------------------------------------------------------------------------------
# wssplot: Function to plot within group sums of squares from k-means clustering data
#---------------------------------------------------------------------------------------
wssplot <- function(data, nc=15, seed=1234){
	wss <- (nrow(data)-1)*sum(apply(data,2,var))
        for (i in 2:nc){
	        set.seed(seed)
                wss[i] <- sum(kmeans(data, centers=i)$withinss)
	}
        plot(1:nc, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
}
