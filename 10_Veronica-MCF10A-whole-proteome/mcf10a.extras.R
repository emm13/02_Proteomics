
# Background
all.prot.names = unlist(strsplit(rownames(exprs(qnt.prot.no.imp)),"; "))
all.prot.entrez = unique(bitr(geneID = exc.2hr.names,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)$ENTREZID)

# Test

# Exclusive to insulin starved cells
exc.starved.names = unique(unlist(strsplit(rownames(exc.starved),"; ")))
exc.starved.entrez = na.omit(unique(bitr(geneID = exc.starved.names,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)$ENTREZID))

# Exclusive to 30mins post insulin stimulation
exc.30min.names = unique(unlist(strsplit(rownames(exc.30min),"; ")))
exc.30min.entrez = na.omit(unique(bitr(geneID = exc.30min.names,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)$ENTREZID))

# Exclusive to 2hr post insulin stimulation
exc.2hr.names = unique(unlist(strsplit(rownames(exc.2hr),"; ")))
exc.2hr.entrez = na.omit(unique(bitr(geneID = exc.2hr.names,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)$ENTREZID))

# GO enrichment

go.comp = compareCluster(list(starved = exc.starved.names,ins.30min = exc.30min.names, ins.2hr = exc.2hr.names),keytype = "UNIPROT",fun="enrichGO",OrgDb = org.Hs.eg.db)
kegg.comp = compareCluster(list(starved = exc.starved.entrez,ins.30min = exc.30min.entrez, ins.2hr = exc.2hr.entrez),fun="enrichKEGG")

go.2hr = enrichGO(gene = exc.2hr.entrez,OrgDb = 'org.Hs.eg.db',keytype = "ENTREZID")
go.30min = enrichGO(gene = exc.30min.entrez,OrgDb = 'org.Hs.eg.db',keytype = "ENTREZID")
go.starved = enrichGO(gene = exc.starved.entrez,OrgDb = 'org.Hs.eg.db',keytype = "ENTREZID")


# Calculate means/variance
exp.com = exprs(combat.qnt.prot)
vars <- apply(exp.com,1,var,na.rm=T)
means <- rowMeans(exp.com,na.rm=T)
cv2 <- vars/means^2

# Plot data distribution
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(means),log(cv2))

# Pick the mean cut-off for top 5% of most variable genes
minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
useForFit <- means >= minMeanForFit # & spikeins

# Fit a generalised linear regression
fit <- glm.fit(cbind( a0 = 1, a1tilde = 1/means[useForFit]),cv2[useForFit])
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

# Repeating previous plot with fit lines
#----------------------------------------
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(means),log(cv2));
xg <- exp(seq( min(log(means[means>0]),na.rm=T), max(log(means),na.rm=T), length.out=1000 ))
vfit <- a1/xg + a0

# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(data) - 1

# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")

# Rank genes
afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- exp.com[varorder,]
oed[is.na(oed)] = 0
# save for the next exercise
save(oed,file="oed.RData")

# repeat previous plot and add dots
#----------------------------------------
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 100 genes
points(log(means[varorder[1:100]]),log(cv2[varorder[1:100]]),col=2)

# Heatmap of most variable genes
m <- data.frame(oed[1:50,])
heatmap(as.matrix(m)/apply(m,1,max),zlim=c(0,1),col=gray.colors(100),Rowv=NA,labRow=rownames(m),scale="none",ColSideColors=rep(c("#edf8b1","#7fcdbb","#2c7fb8"),times=c(3,3,3)))

# Heatmap of least variable genes
n <- data.frame(oed[4900:4962,])
heatmap(as.matrix(n)/apply(n,1,max),zlim=c(0,1),col=gray.colors(100),Rowv=NA,labRow=rownames(n),scale="none",ColSideColors=rep(c("#edf8b1","#7fcdbb","#2c7fb8"),times=c(3,3,3)))


# Vero's list
vero = c("Q09666","O75534","Q96GQ7","Q12906","P46013","Q9Y520","P42696","Q8IY81","Q9BVJ6")
data[vero,]

#---------------------------------------------------------------------------------------------------
# Looking for those proteins that are exclusively present in the starved but absent from stimulated.
#---------------------------------------------------------------------------------------------------

# Looking for those proteins that are exclusively present in one group but not in the other two
exc.starved = exp[which(exp$zero.starved == 0 & exp$zero.2hr > 1 & exp$zero.30min > 1),]
exc.30min = exp[which(exp$zero.starved > 1 & exp$zero.2hr > 1 & exp$zero.30min == 0),]
exc.2hr = exp[which(exp$zero.starved > 1 & exp$zero.2hr == 0 & exp$zero.30min > 1),]

# Merging list to annotate
excl = rbind(cbind(Condition = "Starved.only",exc.starved),cbind(Condition = "30min.only",exc.30min),cbind(Condition="2hr.only",exc.2hr))
excl$Uniprot = rownames(excl)
df2=data.frame(cSplit(excl,"Uniprot","; ", direction="long"))

# Annotate this list
ann = data.frame(queryMany(df2$Uniprot,scopes="uniprot",fields=c("ensembl.gene","uniprot.Swiss-Prot","entrezgene","symbol","go.MF.term","go.CC.term","go.BP.term","interpro.desc"),species=9606))

# Modify some annotations
for(j in c("ensembl","interpro","go.BP","go.CC","go.MF")){
  name = paste(j,"mod",sep=".")
  ann[,name] = sapply(ann[,j],function(x) paste(unlist(x),collapse="; "))
}

# Merge it with original data
excl.ann = merge(ann,df2,by.x="query",by.y = "Uniprot",all=T)

# Trim down data frame and add extra annotations to make it compatible with the RBP capture data
final.excl = excl.ann[order(excl.ann$Condition,excl.ann$symbol),c(17,7,1,11,5,12,18:29,13:16)]
extra.ann = read.delim("Files-to-send/Whole-proteome-highly-variable-id-mapping.tab",sep="\t",header=T)

f = merge(final.excl,extra.ann,by.x = "uniprot.Swiss.Prot",by.y = "Query",all.x=T,all.y=F)
f = f[,c(1:6,23:30,7:22)]
f = f[order(f$Condition,f$Entry.name),]

write.table(f,paste(outdir,"Exclusive-proteins-with-data-and-annotations.txt",sep="/"),sep="\t",row.names=F,quote=F)
f[grep("RNA-binding",f$interpro.mod),]

# Trying to corect for batch with missing values
# Combat doesn't work as within each batch, all could be missing. 
# Modified filter to be more inclusive trying to use this reference (https://support.bioconductor.org/p/75965/)
exp = data[which(data$zero.starved<3 & data$zero.30min <3 & data$zero.2hr <3),1:9] # 5719 proteins where there is atleast 1 value in each condition

# Comparing data spread before and after filtering
apply(exp,2,summary) # Mean and median a lot more similar
apply(exprs(res),2,summary)

# Doing some batch correction now that we have excluded all proteins where one condition is completely missing.
pheno = pData(qnt.prot.no.imp)
edata = exp
feat = featureData(qnt.prot.no.imp)[rownames(exp)]

# Combat model
batch = pheno$Rep
modcombat = model.matrix(~1,data=pheno)

# Batch corrected data. Visit later as with NAs this doesn't work!
combat_edata = sva::ComBat(dat=na.omit(edata), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=T)

# Draw a pca to check batch correction
prot.pca = prcomp(t(na.omit(exprs(qnt.prot.no.imp))))
summary(prot.pca)
combat.pca = prcomp(t(combat_edata),scale=T)
summary(combat.pca)

j <- ggbiplot(prot.pca,choices=c(1,2), var.axes=F, groups = factor(rep(c("starved","Ins.30min","Ins.2hr"),times=c(3,3,3))), circle = T,obs.scale=1,labels=rownames(prot.pca$x))
print(j)

y <- ggbiplot(combat.pca, choices = c(1,2),var.axes=F, groups = factor(rep(c("starved","Ins.30min","Ins.2hr"),times=c(3,3,3))), circle = T,obs.scale=1,labels=rownames(combat.pca$x))
print(y)

# Make a new object to contain combat results
combat.qnt.prot = MSnSet(exprs = as.matrix(combat_edata),pData = pheno, fData = feat[which(feat$Master.Protein.Accessions %in% rownames(combat_edata))])

#----------------------------------------------------
# Finding most variable genes after batch correction
# Function : vargenes
#----------------------------------------------------

var = vargenes(combat.qnt.prot)


#---------------------------------------------------------------------------------------------------------------------
# Filter data for highly expressing - lose 85 proteins
#---------------------------------------------------------------------------------------------------------------------
keep.exprs.noimp <- rowSums(exp>10)>=3
exp.filt <- exp[keep.exprs.noimp,]

#------------------------------------------------------------
# Histogram to check where an expression cut-off can be drawn
#------------------------------------------------------------

# Set up colours
nsamples <- 10
coul <- brewer.pal(10, "Paired")
col = colorRampPalette(coul)(nsamples)

par(mfrow=c(1,2))

r = exp
r[is.na(r)] = 0
plot(density(r[,1]), col=col[1],ylim = c(0,0.05),xlim = c(0,200), lwd=2, las=2, main="", xlab="")
title(main=paste("A. Raw data (n = 0, ",nrow(r)," proteins)", sep=""), xlab="Relative protein abundance")
abline(v=0, lty=3)
for (i in 2:9){
  den <- density(r[,i],na.rm=T)
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright",colnames(r)[1:9], text.col=col, bty="n",cex=1)

rfilt = exp.filt
rfilt[is.na(rfilt)] = 0
plot(density(rfilt[,1]), col=col[1],ylim = c(0,0.05),xlim = c(0,200), lwd=2, las=2, main="", xlab="")
title(main=paste("B. Filtered data (n = 3, ",nrow(rfilt)," proteins)", sep=""), xlab="Relative protein abundance")
abline(v=0, lty=3)
for (i in 2:9){
  den <- density(rfilt[,i],na.rm=T)
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright",colnames(rfilt)[1:9], text.col=col, bty="n",cex=1)


#-------------------------------------------------------
# Function : ComBat
#-------------------------------------------------------

function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE, 
          mean.only = FALSE) 
{
  if (mean.only == TRUE) {
    cat("Using the 'mean only' version of ComBat\n")
  }
  if (length(dim(batch)) > 1) {
    stop("This version of ComBat only allows one batch variable")
  }
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1 + batch)
  cat("Found", nlevels(batch), "batches\n")
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  }
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) {
    mean.only = TRUE
    cat("Note: one batch has only one sample, setting mean.only=TRUE\n")
  }
  n.array <- sum(n.batches)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n.batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if (ncol(design) > (n.batch + 1)) {
      if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[, 
                                                          -c(1:n.batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      }
      else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  NAs = any(is.na(dat))
  if (NAs) {
    cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"), 
        sep = " ")
  }
  cat("Standardizing Data across genes\n")
  if (!NAs) {
    B.hat <- solve(t(design) %*% design) %*% t(design) %*% 
      t(as.matrix(dat))
  }
  else {
    B.hat = apply(dat, 1, Beta.NA, design)
  }
  grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
  if (!NAs) {
    var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
                                                          n.array)
  }
  else {
    var.pooled <- apply(dat - t(design %*% B.hat), 1, var, 
                        na.rm = T)
  }
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
                                                           n.array)))
  cat("Fitting L/S model and finding priors\n")
  batch.design <- design[, 1:n.batch]
  if (!NAs) {
    gamma.hat <- solve(t(batch.design) %*% batch.design) %*% 
      t(batch.design) %*% t(as.matrix(s.data))
  }
  else {
    gamma.hat = apply(s.data, 1, Beta.NA, batch.design)
  }
  delta.hat <- NULL
  for (i in batches) {
    if (mean.only == TRUE) {
      delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
    }
    else {
      delta.hat <- rbind(delta.hat, apply(s.data[, i], 
                                          1, var, na.rm = T))
    }
  }
  gamma.bar <- apply(gamma.hat, 1, mean)
  t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  if (prior.plots & par.prior) {
    par(mfrow = c(2, 2))
    tmp <- density(gamma.hat[1, ])
    plot(tmp, type = "l", main = "Density Plot")
    xx <- seq(min(tmp$x), max(tmp$x), length = 100)
    lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
    qqnorm(gamma.hat[1, ])
    qqline(gamma.hat[1, ], col = 2)
    tmp <- density(delta.hat[1, ])
    invgam <- 1/rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
    tmp1 <- density(invgam)
    plot(tmp, typ = "l", main = "Density Plot", ylim = c(0, 
                                                         max(tmp$y, tmp1$y)))
    lines(tmp1, col = 2)
    qqplot(delta.hat[1, ], invgam, xlab = "Sample Quantiles", 
           ylab = "Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
    title("Q-Q Plot")
  }
  gamma.star <- delta.star <- NULL
  if (par.prior) {
    cat("Finding parametric adjustments\n")
    for (i in 1:n.batch) {
      if (mean.only) {
        gamma.star <- rbind(gamma.star, postmean(gamma.hat[i, 
                                                           ], gamma.bar[i], 1, 1, t2[i]))
        delta.star <- rbind(delta.star, rep(1, nrow(s.data)))
      }
      else {
        temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, 
                                                         ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i], 
                       b.prior[i])
        gamma.star <- rbind(gamma.star, temp[1, ])
        delta.star <- rbind(delta.star, temp[2, ])
      }
    }
  }
  else {
    cat("Finding nonparametric adjustments\n")
    for (i in 1:n.batch) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]), 
                         gamma.hat[i, ], delta.hat[i, ])
      gamma.star <- rbind(gamma.star, temp[1, ])
      delta.star <- rbind(delta.star, temp[2, ])
    }
  }
  cat("Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
                                                       ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
                                                                                                           n.batches[j])))
    j <- j + 1
  }
  bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
                                                        n.array)))) + stand.mean
  return(bayesdata)
}

# Finding highly variable genes


cv.cut = summary(cv2)[5]

minMeanForFit <- unname( quantile( means[ which( cv2 > cv.cut ) ], .95 ) )
useForFit <- means >= minMeanForFit # & spikeins

# Fit a generalised linear regression 
fit <- glm.fit(cbind( a0 = 1, a1tilde = 1/means[useForFit]),cv2[useForFit])
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

# Repeating previous plot
smoothScatter(log(means),log(cv2),main="Colour density representation of mean and variance with fit and conf intervals")
xg <- exp(seq( min(log(means[means>0]),na.rm=T), max(log(means),na.rm=T), length.out=1000 ))
vfit <- a1/xg + a0

# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(data) - 1

# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")

# Rank genes
afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- exp[varorder,]
oed[is.na(oed)] = 0

# repeat previous plot and add dots
smoothScatter(log(means),log(cv2),main="Colour density representation of mean and variance showing most variable (50) proteins")
lines( log(xg), log(vfit), col="black",lwd=3 )
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")




# Checking the effects of normalisation on covariance

# Non-normalised data
sd1 <- apply(log2(exprs(impute.res)),1,sd)
mn1 <- apply(log2(exprs(impute.res)),1,mean)
cv1 <- sd1/mn1

# Normalised data
sd2 <- apply(exprs(qnt.vsn),1,sd)
mn2 <- apply(exprs(qnt.vsn),1,mean)
cv2 <- sd2/mn2

# Merge the two datasets
dfr <- rbind(data.frame(rank=order(mn1),cv=cv1,norm="raw"),data.frame(rank=order(mn2),cv=cv2,norm="vsn"))

rmed1 <- rollapply(cv1,7,function(x) median(x,na.rm=TRUE))
rmed2 <- rollapply(cv2,7,function(x) median(x,na.rm=TRUE))
dfr2 <- rbind(data.frame(x=seq(0,30,by=30/length(rmed1))[-1],y=rmed1,norm="raw"),data.frame(x=seq(0,30,by=30/length(rmed2))[-1],y=rmed2,norm="vsn"))

p <- ggplot()+geom_line(data=dfr2,aes(x=x,y=y,col=norm)) + theme_gray(7)
plot(p)



#clusdat = comprots.ann[which(common.prots$Source != "DIA_ISOQUANT"),5:16]
#rownames(clusdat) = paste(comprots.ann$Accession[which(common.prots$Source != "DIA_ISOQUANT")],comprots.ann$Source[which(common.prots$Source != "DIA_ISOQUANT")],sep=".")

#commonvar = vargenes(clusdat,"Common-to-callers",rep(c(1,2,3,4),each=3),"symbol")
# Plotting a trend plot for each gene and for all callers
# Includes Isoquant though this is in fmol rather than unitless
#pdf(paste(outdir,paste(format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"Filtered-common-genes-across-callers.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
# PLot of gene-wise trends for all 4 callers
#gg1 = ggplot(data=com.melt,aes(x=Tmt, y=Exp, group = Tmt,fill=Tmt)) + geom_violin() + theme(axis.text.x=element_text(angle=90, hjust=1),legend.position="top")+facet_wrap(~Accession)
#gg10 <- facet_multiple(plot=gg1, facets="Accession", ncol = 4, nrow = 4, scales = "free_y")
#dev.off()


#---------------------------------------------------------------------------------------------------------------
# Sneaky use of maSigPro to detect trends in the data
#---------------------------------------------------------------------------------------------------------------

seqdat = filt.4[which(filt.4$Source == "DDA_SCAFFOLD"),5:16]
rownames(seqdat) = filt.4$Accession[which(filt.4$Source == "DDA_SCAFFOLD")]

des = as.matrix(cbind(Time = rep(c(0,3,6,9),each=3),Rep = rep(c(1,2,3,4),each=3),NCL = rep(c(1,0),times=c(3,9)),U = rep(c(0,1,0),times=c(3,3,6)),S = rep(c(0,1,0),times=c(6,3,3)),IT = rep(c(0,1),times=c(9,3))))
rownames(des) = colnames(seqdat)

e = make.design.matrix(des)

# Run maSigPro and return lists of genes
library(maSigPro)
ma.res = runMaSigPro(seqdat,e,"NCL")
ma.res

# For see.genes, clusters = 4, cluster.method = hclust, agglo.method = "complete", distance = "cor"
sg = see.genes(seqdat[ma.res$summary,], edesign = des.series,k=4, cluster.method="hclust", main = "Time Course", show.fit = F, dis = e$dis, groups.ve