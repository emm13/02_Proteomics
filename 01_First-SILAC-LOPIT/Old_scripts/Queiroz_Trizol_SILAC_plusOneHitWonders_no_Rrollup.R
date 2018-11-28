
#############Trizol enriched RNA-BPs_ SILAC quan#############
##################Welcome Trust - TTT group##################
#################Rayner Queiroz / June 2017##################

setwd("C:/Users/rq214/Desktop/R_Analysis_Trizol")
source('C:/Users/rq214/Desktop/R_Analysis_Trizol/scripts/RRollupMod.R', encoding = 'UTF-8')
source('C:/Users/rq214/Desktop/R_Analysis_Trizol/scripts/LimmaRankProd.R', encoding = 'UTF-8')
library(Hmisc)


############## Preparing data
## Normalization and Rollup, uses RRollup from DanteR
par(mfrow=c(1,1))

NumReps = 3 # Number of biological replicates 
numNAcond = 0 # maximum number of missing values (missing condition) permited per replicate
#numNAmax = 4 # maximum number of missing values permited (2 donors missing)
#numNAmax_II = 6 # maximum number of missing values permited (3 donors missing)
PSM=uniquePSM=PSMsM=tDat=PSM_D=Pep_D=Prot_D=NAs=OneHitWond=OnlyLight=OnlyHeavy=OnlyHeavy_Prot=OnlyLight_Prot=OnlyLightTotData=OnlyHeavyTotData=NULL

##Normalization and removal of contaminating proteins 
Contaminant = read.csv("Common contaminant_all.csv")

for ( i in 1:NumReps ) {
  cat( "Processing replicate ", i, "\n")
  #PSM[[i]] = read.csv(paste("Trizol_150mJ_rep",i,".csv",sep=""))
  #PSM[[i]] = read.csv(paste("Trizol_275mJ_rep",i,".csv",sep=""))
  #PSM[[i]] = read.csv(paste("Trizol_400mJ_rep",i,".csv",sep=""))

  #filter unique PSMs 
  uniquePSM[[i]] = PSM[[i]][PSM[[i]]$Quan.Info=="Unique",] ##RQ: look here the importance of pep rollup! quans from heavy are not quite the same from quans from light when each of them is id`ed`
  
  uniquePSM[[i]]$Sequence = toupper(uniquePSM[[i]]$Sequence) # RQ: redundant for PD v.2.1
  PSMsM[[i]] = uniquePSM[[i]][,c(3,10,16,17)]  # RQ: check column name in the PD version used
  colnames(PSMsM[[i]]) = c("Sequence","Master.Protein.Accessions","Light","Heavy")
  
  #Remove PSMs from Common Contaminant Proteins 
  PSMsM[[i]] = PSMsM[[i]][! PSMsM[[i]][,2]%in% Contaminant[,1] ,] 
  
  #log2/median normalization 
  tDat[[i]]=log (PSMsM[[i]][,3:4],2)
  tDat[[i]]<-t(t(tDat[[i]])-apply(tDat[[i]],2,median,na.rm=T))
  boxplot(tDat[[i]],ylab="log(Intensity)",xlab="Label")
  PSMsM[[i]][,3:4]=tDat[[i]]
}


#Split by donor
PSM_D[[1]] = PSMsM[[1]] 
PSM_D[[2]] = PSMsM[[2]]
PSM_D[[3]] = PSMsM[[3]] 

#Count NAs per replicate and remove PSMs below threshold.
for (i in 1:NumReps) {
  cat( "Processing replicate ", i, "\n")
  PSM_D[[i]] = PSM_D[[i]] [rowSums(is.na(PSM_D[[i]][,3:4]))<= numNAcond , ]
}
for (i in 1:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  NAs[[i]] = PSMsM[[i]] [rowSums(is.na(PSMsM[[i]][,3:4]))> numNAcond , ]
}

#Histograms
par(mfrow=c(3,2))
hist(PSM_D[[1]][,3],100,main="Rep1 Light",xlab="Intensity",cex.main=0.8)
hist(PSM_D[[1]][,4],100,main="Rep1 Heavy",xlab="Intensity",cex.main=0.8)

hist(PSM_D[[2]][,3],100,main="Rep2 Light",xlab="Intensity",cex.main=0.8)
hist(PSM_D[[2]][,4],100,main="Rep2 Heavy",xlab="Intensity",cex.main=0.8)

hist(PSM_D[[3]][,3],100,main="Rep3 Light",xlab="Intensity",cex.main=0.8)
hist(PSM_D[[3]][,4],100,main="Rep3 Heavy",xlab="Intensity",cex.main=0.8)
par(mfrow=c(1,1))


#Peptide to protein 
for  (i in 1:NumReps){
  cat( "Processing replicate ", i, "\n")
  ##average PSMs from same sequences
  #Pep_D[[i]] = RRollup(PSM_D[[i]][,3:4],PSM_D[[i]][,1],Mode="mean",oneHitWonders=T,center=F,minPep=1)$RolledUp
  #Pep_D[[i]] = merge(Pep_D[[i]],PSM_D[[i]][,1:2],by.x=0,by.y=1)
  #Pep_D[[i]] = unique(Pep_D[[i]])
  #Pep_D[[i]] = Pep_D[[i]][!duplicated(Pep_D[[i]][,1]),]
  
  ##proper rollup to proteins
  #Prot_D[[i]] = RRollup(Pep_D[[i]][,2:3],Pep_D[[i]][,4],Mode="mean",oneHitWonders=T,center=F,minPep=1)$RolledUp
  #Prot_D[[i]] = Prot_D[[i]]-rowMeans(Prot_D[[i]],na.rm=T) # another normalization?!
  
  Pep_D[[i]] = PSM_D[[i]][,2:4]
  Prot_D[[i]] = aggregate(.~Master.Protein.Accessions, data = Pep_D[[i]], FUN = mean )
}

# merge replicates
#Prot_D[[1]] = cbind(Row.Names=rownames(Prot_D[[1]]),Prot_D[[1]])
TotData = Prot_D[[1]]

for (i in 2:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  #Prot_D[[i]] = cbind(Row.Names=rownames(Prot_D[[i]]),Prot_D[[i]])
  rownames(Prot_D[[i]]) <- NULL
  TotData= merge(TotData,Prot_D[[i]],by="Master.Protein.Accessions",all=T)
  
}
rownames(TotData) = TotData[,1]
TotData = TotData[,2:(dim(TotData)[2])]  

#########################################################
#Dealing with One Hit Wonders

for (i in 1:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  OnlyLight[[i]] =  NAs[[i]][is.na(NAs[[i]]$Heavy),]
  OnlyHeavy[[i]] =  NAs[[i]][is.na(NAs[[i]]$Light),]
  OnlyLight[[i]] =  OnlyLight[[i]][,c(2,3)]
  OnlyHeavy[[i]] =  OnlyHeavy[[i]][,c(2,4)]
  
  
  #Impute minimal detected abandance in missing conditions 
  ##RQ: imputation only in missing condition, but not for missing replicate, as for a SILAC experiment a missing condition is in fact abscence or too low abundance 
  OnlyHeavy[[i]]$Light = min(PSM_D[[i]]$Light)
  OnlyLight[[i]]$Heavy = min(PSM_D[[i]]$Heavy)

  OnlyLight_Prot[[i]] = aggregate(Light  ~ Master.Protein.Accessions , data = OnlyLight[[i]], FUN = mean )
  OnlyLight_Prot[[i]]$Heavy = min(OnlyLight[[i]]$Heavy)
 
  OnlyHeavy_Prot[[i]] = aggregate(Heavy  ~ Master.Protein.Accessions , data = OnlyHeavy[[i]], FUN = mean )
  OnlyHeavy_Prot[[i]]$Light = min(OnlyHeavy[[i]]$Light)
}
 
  # merge replicates OnlyLight
  OnlyLightTotData = OnlyLight_Prot[[1]]
  rownames(OnlyLightTotData) = OnlyLightTotData[,1]
  
for (i in 2:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  OnlyLightTotData= merge(OnlyLightTotData,OnlyLight_Prot[[i]],by="Master.Protein.Accessions",all=T)
}

  # merge replicates OnlyHeavy
  OnlyHeavyTotData = OnlyHeavy_Prot[[1]]
  rownames(OnlyHeavyTotData) = OnlyHeavyTotData[,1]
  
for (i in 2:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  OnlyHeavyTotData= merge(OnlyHeavyTotData,OnlyHeavy_Prot[[i]],by="Master.Protein.Accessions",all=T)
}
  rownames(OnlyHeavyTotData) = OnlyHeavyTotData[,1]
  OnlyHeavyTotData = OnlyHeavyTotData[,2:(dim(OnlyHeavyTotData)[2])]
  rownames(OnlyLightTotData) = OnlyLightTotData[,1]
  OnlyLightTotData = OnlyLightTotData[,2:(dim(OnlyLightTotData)[2])]

#Combine tables and average with imputed data
CombData=NULL
CombData=rbind(TotData,OnlyHeavyTotData,OnlyLightTotData)
CombData = aggregate(.~rownames(CombData), data = CombData, FUN = mean ,na.action = na.pass)
rownames(CombData) = CombData[,1]
CombData = CombData[,2:(dim(CombData)[2])]

#Calculate percentage of missing values
#sum(is.na(TotData))/prod(dim(TotData))*100
sum(is.na(TotData))/prod(dim(TotData))*100
sum(is.na(OnlyHeavyTotData))/prod(dim(OnlyHeavyTotData))*100
sum(is.na(OnlyLightTotData))/prod(dim(OnlyLightTotData))*100
sum(is.na(CombData))/prod(dim(CombData))*100

#create ratios
CombData$"L/H-Rep1" = CombData[,1]-CombData[,2]
CombData$"L/H-Rep2" = CombData[,3]-CombData[,4]
CombData$"L/H-Rep3" = CombData[,5]-CombData[,6]


par(mfrow=c(3,1))
hist(CombData[,7],100,main="Rep1 Ratio Light/Heavy",xlab="Intensity (log2/median)",cex.main=0.8)
hist(CombData[,8],100,main="Rep2 imputed Ratio Light/Heavy",xlab="Intensity (log2/median)",cex.main=0.8)
hist(CombData[,9],100,main="Rep3 Ratio Light/Heavy",xlab="Intensity (log2/median)",cex.main=0.8)
par(mfrow=c(1,1))

##Assign RNA-BPs: Proteins with log2 L/H ratio above 0  
RegProt=RegData= NULL
for (i in 1:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  RegProt[[i]] = CombData[,c(i,i+1,i+6)]
  RegProt[[i]] = RegProt[[i]] [RegProt[[i]][,3]> 0 , ]
  RegProt[[i]] = RegProt[[i]] [rowSums(!is.na(RegProt[[i]][,1:3])) == 3 , ]
}

# merge replicates
RegProt[[1]] = cbind(Row.Names=rownames(RegProt[[1]]),RegProt[[1]])
RegData = RegProt[[1]]

for (i in 2:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  RegProt[[i]] = cbind(Row.Names=rownames(RegProt[[i]]),RegProt[[i]])
  #rownames(RegProt[[i]]) <- NULL
  RegData= merge(RegData,RegProt[[i]],by="Row.Names",all=T)
}

##Assign Trizol noise: Proteins with log2 L/H ratio below or equal to 0  
NoiseProt=NoiseData= NULL
for (i in 1:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  NoiseProt[[i]] = CombData[,c(i,i+1,i+6)]
  NoiseProt[[i]] = NoiseProt[[i]] [NoiseProt[[i]][,3]<= 0 , ]
  NoiseProt[[i]] = NoiseProt[[i]] [rowSums(!is.na(NoiseProt[[i]][,1:3])) == 3 , ]
}

# merge replicates
NoiseProt[[1]] = cbind(Row.Names=rownames(NoiseProt[[1]]),NoiseProt[[1]])
NoiseData = NoiseProt[[1]]

for (i in 2:NumReps) { 
  cat( "Processing replicate ", i, "\n")
  NoiseProt[[i]] = cbind(Row.Names=rownames(NoiseProt[[i]]),NoiseProt[[i]])
  #rownames(RegProt[[i]]) <- NULL
  NoiseData= merge(NoiseData,NoiseProt[[i]],by="Row.Names",all=T)
}

#Store data
#write.csv(CombData,file=paste("150mJ_Trizol_normData",numNAcond,"missing_Cond.csv",sep=""))
#write.csv(RegData,file=paste("150mJ_Trizol_RNA-BPs",numNAcond,"missing_Cond.csv",sep=""))
#write.csv(NoiseData,file=paste("150mJ_Trizol_noise",numNAcond,"missing_Cond.csv",sep=""))
#write.csv( OnlyHeavyTotData ,file=paste("150mJ_Trizol_onlyHeavy.csv"))
#write.csv( OnlyLightTotData ,file=paste("150mJ_Trizol_onlyLight.csv"))

#write.csv(CombData,file=paste("275mJ_Trizol_normData",numNAcond,"missing_Cond.csv",sep=""))
#write.csv(RegData,file=paste("275mJ_Trizol_RNA-BPs",numNAcond,"missing_Cond.csv",sep=""))
#write.csv(NoiseData,file=paste("275mJ_Trizol_noise",numNAcond,"missing_Cond.csv",sep=""))
#write.csv( OnlyHeavyTotData ,file=paste("275mJ_Trizol_onlyHeavy.csv"))
#write.csv( OnlyLightTotData ,file=paste("275mJ_Trizol_onlyLight.csv"))

#write.csv(CombData,file=paste("400mJ_Trizol_normData",numNAcond,"missing_Cond.csv",sep=""))
#write.csv(RegData,file=paste("400mJ_Trizol_RNA-BPs",numNAcond,"missing_Cond.csv",sep=""))
#write.csv(NoiseData,file=paste("400mJ_Trizol_noise",numNAcond,"missing_Cond.csv",sep=""))
#write.csv( OnlyHeavyTotData ,file=paste("400mJ_Trizol_onlyHeavy.csv"))
#write.csv( OnlyLightTotData ,file=paste("400mJ_Trizol_onlyLight.csv"))


################### Statistics###################################################
#For details, see Schw?mmle, V.; Leon, I. R. and Jensen, O. N. Assessment and
#improvement of statistical tools for comparative proteomics analysis of sparse
#data sets with few experimental replicates J Proteome Res, 2013, 12, 3874-3883.
#################################################################################

### Set parameters
qlim = 0.05 #q-value stringency
NumCond = 2 #Number of conditions
NumReps = 3 #Set number of replicates
#ChosenTable =  # Compiled table of normalized data
#ChosenTable = read.csv(paste(".csv"))
#RefCond = 1 # Condition to which the data is compared (150mJ=1 / 275mJ=2 / 400mJ=3)
####

##prep data tables
MAData =Against=DataSet=RR= NULL
DataSet = ChosenTable #[,c(2:11)]
Against = seq(RefCond,NumCond*NumReps,NumCond) 
RR = 1:(NumReps*NumCond)
RR = RR[-Against]
RR = rbind(RR,rep(Against,each=NumCond-1))
MAData = NULL
for (i in 1:ncol(RR))
  MAData = cbind(MAData,DataSet[,RR[1,i]]-DataSet[,RR[2,i]])
rownames(MAData) = rownames(ChosenTable)
MAData = MAData[,1:(dim(MAData)[2])]
MAReps<-rep(1:(NumCond-1),NumReps)
#calculating q-values and vulcano plots
qvalues<-LimmaRankProd(MAData,MAReps)

##STOP HERE!!!!

#########################################4


HiLoList = AllHiLo = NULL

for (i in 1:(NumCond-1)) {
  HiLoList[[i]]<-(qvalues[[i]])[!is.na(qvalues[[i]][,2]) & !is.na(qvalues[[i]][,3]) & (qvalues[[i]][,2]<qlim | qvalues[[i]][,3] < qlim),]
  colnames(HiLoList[[i]])<-c("t-test","limma","rank products","log2 change")
  print(paste(i,":",dim(HiLoList[[i]])))
}
#Compiling statistics table
AllHiLo = HiLoList[[1]]



colnames(AllHiLo)<-c("Protein Accession",
                     paste("Cond",rep(c(2:5) #It must be corrected if it is not against Condition 1!!!
                                      ,each=4),"vs Cond",RefCond,c("t-test","limma","rank products","log2 change")))


#write.csv(AllHiLo,file=paste("min1pep_Reg_Proteins_RefCondition",RefCond,"_atLeast2of5Donors.csv",sep=""))


write.csv(AllHiLo,file=paste("Reg_TotalProteins_",numNAcond,"missingCond_24h.csv",sep=""))

