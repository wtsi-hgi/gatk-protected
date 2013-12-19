getFDR<- function(rawTable,id, lenient = F) {
  
  
  maxAC<- max(as.numeric(rawTable$resource.oneKGAC),na.rm=T)
  incr<- 0.05
  bins <-unique(round(10^(seq(from=0.0,to=log10(maxAC),by=incr))))
  
  TPs = vector(length=length(bins),mode="numeric")
  Points=FPs=TPs
  
  totalTP = nrow(subset(rawTable,FILTER=="PASS" & AC > 0))
  totalFP = nrow(subset(rawTable,FILTER=="LowQual" & AC==0))
  for (ac in c(1:maxAC)) {
    acSlice = subset(rawTable,resource.oneKGAC == ac )
   # acSlice = acSlice[grep(",",acSlice$ALT,invert=T),]
    binIdx = max(which(ac>=bins))

    tpSet = (subset(acSlice, FILTER=="PASS" & AC > 0))
    fpSet = (subset(acSlice, FILTER=="LowQual" & AC==0 ))
    
    tp = nrow(tpSet)
    fp = nrow(fpSet)
    if (lenient == F) {
      
    }
    
    TPs[binIdx] = TPs[binIdx]+ tp
    FPs[binIdx] = FPs[binIdx] + fp 
    Points[binIdx] =  Points[binIdx]+tp+fp
    
  }
  allTPs = sum(TPs)
  allFPs = sum(FPs)
  allFDR = allFPs/(allFPs+allTPs)
  
  fdrData<-data.frame(TP=TPs,FP=FPs,
                      FDR=FPs/(TPs+TPs),
                      AFBin = bins/maxAC,
                      ACBin=bins,
                      Points=Points,
                      dataSet=id,
                      oneKGFDR = allFDR, 
                      oneKGTPs = allTPs,
                      oneKGFPs = allFPs,
                      oneKGCalls = allTPs+allFPs,
                      totalTPs = totalTP,
                      totalFPs = totalFP,
                      totalCalls = totalTP+totalFP,
                      totalFDR = totalFP/(totalTP+totalFP))
  fdrData
}


library(ggplot2)

pbase <- "/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/finalRunGGANew."

# dafsnp_noref <- read.table(paste(pbase,"afSNPs.noRef.filtered.annotated.table",sep=""),header=T)
# dafindel_noref <- read.table(paste(pbase,"afIndels.noRef.filtered.annotated.table",sep=""),header=T)
# dunifsnp_noref <- read.table(paste(pbase,"unifSNPs.noRef.filtered.annotated.table",sep=""),header=T)
# dunifindel_noref <- read.table(paste(pbase,"unifIndels.noRef.filtered.annotated.table",sep=""),header=T)
# fdrAFSNP_noref <- getFDR(dafsnp_noref, "AFSNP")
# fdrAFINDEL_noref <- getFDR(dafindel_noref, "AFINDEL")
# fdrUNIFSNP_noref <- getFDR(dunifsnp_noref, "UNIFSNP")
# fdrUNIFINDEL_noref <- getFDR(dunifindel_noref, "UNIFINDEL")
# 
# dataFDRnoRef <- rbind(fdrAFSNP_noref,fdrAFINDEL_noref)
# dataFDRnoRef <- rbind(dataFDRnoRef,fdrUNIFSNP_noref)
# dataFDRnoRef <- rbind(dataFDRnoRef,fdrUNIFINDEL_noref)

# with ref sample
dafsnp_withRef <- read.table(paste(pbase,"afSNPs.withRef.filtered.annotated.table",sep=""),header=T)
dafindel_withRef <- read.table(paste(pbase,"afIndels.withRef.filtered.annotated.table",sep=""),header=T)
dunifsnp_withRef <- read.table(paste(pbase,"unifSNPs.withRef.filtered.annotated.table",sep=""),header=T)
dunifindel_withRef <- read.table(paste(pbase,"unifIndels.withRef.filtered.annotated.table",sep=""),header=T)
#dlostimp_withRef <- read.table(paste(pbase,"lostToImputation.withRef.filtered.annotated.table",sep=""),header=T)

fdrUNIFSNP_withRef <- getFDR(dunifsnp_withRef, "SNPs, Uniformly Distributed",F)
fdrAFSNP_withRef <- getFDR(dafsnp_withRef, "SNPs, AF Distributed",F)
fdrAFINDEL_withRef <- getFDR(dafindel_withRef, "Indels, AF Distributed",F)
fdrUNIFINDEL_withRef <- getFDR(dunifindel_withRef, "Indels, Uniformly Distributed",F)
#fdrLOSTIMP_withRef = getFDR(dlostimp_withRef, "SNPs lost to imputation",F)

dataIndelFDRwithRef <- rbind(fdrUNIFINDEL_withRef,fdrAFINDEL_withRef)
dataSNPFDRwithRef <- rbind(fdrAFSNP_withRef,fdrUNIFSNP_withRef)
dataUnifFDRwithRef <- rbind(fdrUNIFINDEL_withRef,fdrUNIFSNP_withRef)
dataAFFDRwithRef <- rbind(fdrAFSNP_withRef,fdrAFINDEL_withRef)

dataAFFDRwithRef$log10TP = log10(dataAFFDRwithRef$TP)
dataUnifFDRwithRef$log10TP = log10(dataUnifFDRwithRef$TP)

da<-ggplot(subset(dataAFFDRwithRef,is.finite(FDR)),
           aes(x=ACBin,y=FDR,color=dataSet,weight=TP, size = TP))
da<-da+geom_point()
da <- da + scale_x_log10("1000 Genomes Phase 1 Allele Count")
da <- da + stat_smooth(method="lm", formula=y ~ poly(x,3))
da <- da + opts(title="1000 Genomes Lowpass False Discovery rate")
png("fdr_lsv_vs_af.png",width=11,height=8.5,units="in",res=600)
print(da)
dev.off()


#da<-ggplot(subset(dataUnifFDRwithRef,is.finite(FDR)),
#           aes(x=ACBin,y=FDR,color=dataSet,weight=TP, size = TP))
#da<-da+geom_point()
#da <- da + scale_x_log10("1000 Genomes Phase 1 Allele Count")
#da <- da + stat_smooth(method="lm", formula=y ~ poly(x,3))
#da <- da + opts(title="1000 Genomes Lowpass False Discovery rate")
#png("fdr_lsv_vs_unif.png",width=640,height=480)
#print(da)
#dev.off()

#LOF
#dlofsnp_withRef <- read.table(paste(pbase,"LOFSNP.withRef.filtered.annotated.LOF.table",sep=""),header=T)
#dlofindel_withRef <- read.table(paste(pbase,"LOFINDEL.withRef.filtered.annotated.LOF.table",sep=""),header=T)
fdrLOFSNP_withRef <- getFDR(dlofsnp_withRef, "LOFSNP")
fdrLOFINDEL_withRef <- getFDR(dlofindel_withRef, "LOFINDEL")

# # controls
# domnipoly_withRef <- read.table(paste(pbase,"omniPoly.withRef.filtered.annotated.table",sep=""),header=T)
# domnimono_withRef <- read.table(paste(pbase,"omniMono.withRef.filtered.annotated.table",sep=""),header=T)
# dexome_withRef <- read.table(paste(pbase,"exomeChip.withRef.filtered.annotated.table",sep=""),header=T)
# dmillspoly_withRef <- read.table(paste(pbase,"millsPoly.withRef.filtered.annotated.table",sep=""),header=T)
# fdrOMNIPOLY_withRef <- getFDR(domnipoly_withRef, "OMNIPOLY")
# fdrOMNIMONO_withRef <- getFDR(domnimono_withRef, "OMNIMONO")
# fdrEXOME_withRef <- getFDR(dexome_withRef, "EXOMECHIP")
# fdrMILLSPOLY_withRef <- getFDR(dmillspoly_withRef, "MILLSPOLY")
# 

## specificity stuff
dsnpsatBaits <- read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/dataAnalysis/results/UGDiscoveryAtAllBaits.ref.snp.recalibrated.annotatedWithPhase1AC.table",header=T)

# show, for each AC, fracion of called sites which are poly in 1000G
maxAC = max(dsnpsatBaits$AC,na.rm=T)

incr<- 0.05
bins <-unique(round(10^(seq(from=0.0,to=log10(maxAC),by=incr))))
numVariants = numIn1000G = vector(length=length(bins),mode="numeric")

for (ac in 1:maxAC) {
  acSlice = subset(dsnpsatBaits,AC == ac & FILTER=="PASS" & TYPE == "SNP")
  binIdx = max(which(ac>=bins))
  numIn1000G[binIdx] = numIn1000G[binIdx] + dim(subset(acSlice,is.finite(resource.AC)))[1] 
  numVariants[binIdx] = numVariants[binIdx] + dim(acSlice)[1]
  
}

data1000G = data.frame(acBins = bins, afBins = bins/maxAC, numVariants = numVariants, numIn1000G = numIn1000G,
                       fracIn1000G = numIn1000G/numVariants)

p <- ggplot(data1000G, aes(x=afBins, y = fracIn1000G))+geom_point() + geom_smooth()
p <- p + scale_x_log10("Allele Frequency")
p <- p + scale_y_continuous("Fraction on variants in 1000 Genomes")

png("frac_1kg_vs_af.png",width=11,height=8.5,units="in",res=600)
print(p)
dev.off()

