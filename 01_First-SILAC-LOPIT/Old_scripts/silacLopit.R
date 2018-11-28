#---------------------------------------------------------------------------
# Author 	      : Manasa Ramakrishna, mr325@le.ac.uk
# Date started 	: 1st June, 2017
# Last modified : 7th June, 2017
# Aim 		      : To take a look at first SILAC labelled LOPIT data on Trizol
# Depends       : On 'silacFunctions.R'. Make sure they are in the same directory
# Notes         : Works on data from Rayner's first experiments
#---------------------------------------------------------------------------

# Invoking libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite("Mus.musculus")
library(outliers)
library(reshape2)
library(ggplot2)


#Setting working directories
wd = "/Users/manasa/Documents/Work/TTT/02_Proteomics/01_First-SILAC-LOPIT/"
setwd(wd)
getwd()

indir = paste(wd,"Input",sep="/")
outdir = paste(wd,paste(Sys.Date(),"Output",sep = "_"),sep = "/")

if (exists(outdir)){
  print("Outdir exists")
}else{
  dir.create(outdir)
}

# Sourcing function file
source("RRollupMod.R")


# Read data

# File of contaminants - proteins to exclude from analysis as are things like keratin, alcohol dehydrogenase etc....
contam = read.delim("Input/Common contaminant_all.csv",sep=",",header=T)

# Data files
infiles = grep("Trizol",list.files("Input/",full.names = T),value=T)
prot.data = NULL
for (i in infiles){
  in.dat = read.delim(i,sep="\t",comment.char="",as.is=T,header=F)
  in.dat$sample = strsplit(i,"//")[[1]][2]
  print(names(in.dat))
  prot.data = rbind(prot.data,in.dat)
}

colnames(prot.data) = prot.data[1,]
dim(prot.data)
head(prot.data)

# Remove header lines as they differ in one of the columns
remove.head = which(prot.data[,1]=="Checked")
prot.data = prot.data[-(remove.head),]
colnames(prot.data)
dim(prot.data)

# Change header names a little to make them neutral and remove space, special characters
colnames(prot.data) = tolower(colnames(prot.data))
colnames(prot.data) = gsub(" ",".",colnames(prot.data))
colnames(prot.data) = gsub("#","no",colnames(prot.data))
colnames(prot.data) = gsub("\\:\\.f2\\:","n",colnames(prot.data))
colnames(prot.data)[12] = "theoretical.mass"
colnames(prot.data)[13] = "light.sample"
colnames(prot.data)[14] = "heavy.sample"
colnames(prot.data)[15] = "abundance.ratio.heavy.to.light"
colnames(prot.data)[16] = "abundance.light"
colnames(prot.data)[17] = "abundance.heavy"
colnames(prot.data)[29] = "sample"

# Add rep, reagent and UV amount columns
prot.data$uv = sapply(strsplit(prot.data$sample,"_"),"[[",2)
prot.data$rep = gsub(".txt","",gsub("rep","",sapply(strsplit(prot.data$sample,"_"),"[[",3)))
prot.data$rep = paste(prot.data$uv,prot.data$rep,sep=".")
prot.data$reagent = sapply(strsplit(prot.data$sample,"_"),"[[",1)
head(prot.data)

# ------------------
# Step 1 : Filter 
# ------------------

# Step 1a : Filter only for those peptides that have a unique master protein. Done using column "quan.info"
dim(prot.data)
filt.1a = prot.data[which(prot.data$quan.info == "Unique"),]
dim(filt.1a)
peptide.stats = cbind("total"= table(prot.data$sample), "unique" = table(filt.1a$sample), "non.unique"=table(prot.data$sample)-table(filt.1a$sample))
peptide.stats


head(filt.1a)
table(filt.1a$sequence,filt.1a$sample)

# This table is very odd. Rayner had an explanation but I'm still muddled
table(filt.1a$heavy.sample,filt.1a$light.sample)

# Making two new columns "heavy.mod" and "light.mod" where "High" is converted to "Found"
# Table makes a bit more sense now
# Light = "Cross linked"; Heavy = "Non-crosslinked"
# 15810 peptides across all reps found in Light and not found in Heavy
# 1672 peptides across all reps found in Heavy and not found in Light - contaminants ?
# 16897 peptides across all reps found in both - RBPs ?

filt.1a$heavy.mod = filt.1a$heavy.sample
filt.1a$light.mod = filt.1a$light.sample
filt.1a$light.mod = gsub("Not Found","not.found",filt.1a$light.mod)
filt.1a$light.mod = gsub("Peak Found","peak.found",filt.1a$light.mod)
filt.1a$heavy.mod = gsub("Not Found","not.found",filt.1a$heavy.mod)
filt.1a$heavy.mod = gsub("Peak Found","peak.found",filt.1a$heavy.mod)
filt.1a$heavy.mod[which(filt.1a$heavy.mod == "High")] = "peak.found"
filt.1a$light.mod[which(filt.1a$light.mod == "High")] = "peak.found"
table(light=filt.1a$light.mod,heavy=filt.1a$heavy.mod)

# Moving on......
head(filt.1a)

# Step 1b : Filter out those proteins that are contaminants from the contaminants list. 
filt.1b = filt.1a[-which(filt.1a$master.protein.accessions %in% contam$Protein.Group.Accessions),]
num.contams = length(which(filt.1a$master.protein.accessions %in% contam$Protein.Group.Accessions))
dim(filt.1a)
dim(filt.1b)

# ------------------
# Step 2 : Normalise 
# heavy = non-crosslinked
# light = crosslinked
# ------------------

# Annotate which peptides are missing heavy, light or both, abundance values
filt.1b$missing.val = rowSums(is.na(filt.1b[,c("abundance.heavy", "abundance.light")])) > 0

# log transform data 
filt.1b$abundance.heavy = as.numeric(filt.1b$abundance.heavy)
filt.1b$abundance.light = as.numeric(filt.1b$abundance.light)

filt.1b$heavy.log = log(filt.1b$abundance.heavy,2)
filt.1b$light.log = log(filt.1b$abundance.light,2)

# Generate an abundance ratio which for log transformed data is a subtraction
filt.1b$norm.abundance.ratio = filt.1b$light.log - filt.1b$heavy.log

# Data is normalised
norm.data = filt.1b
head(norm.data)

# ------------------
# Step 3 : Roll-up 
# Technique where peptides are collapsed into proteins
# ------------------


# Filter potential Rbps - i.e where peak.found in light and not in heavy
pot.rbps = filt.1a[which(filt.1a$light.mod =="peak.found" & filt.1a$heavy.mod == "not.found"),]
dim(pot.rbps)
dim(table(pot.rbps$sequence,pot.rbps$rep)) # 7071 peptide seqs
dim(table(pot.rbps$master.protein.accessions,pot.rbps$rep)) # 1500
rbp.prot.tab = table(pot.rbps$sequence,pot.rbps$rep)
write.table(rbp.prot.tab,"2017-06-01_Output/Potential-RBP-by-replicate.txt",sep="\t",quote=F,row.names=T)


library(limma)
vennDiagram(rbp.prot.tab[,1:3])
vennDiagram(rbp.prot.tab[,4:6])
vennDiagram(rbp.prot.tab[,7:9])





