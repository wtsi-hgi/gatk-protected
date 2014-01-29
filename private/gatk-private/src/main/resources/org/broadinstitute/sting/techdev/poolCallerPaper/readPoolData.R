rm(list=ls())
writeData <- T
processData <-function(input,dataSet,pool,inPoolCaller=F,capNA=T) {
  da<-data.frame(CHROM=input$CHROM,POS=input$POS,AC=input$AC,REF=input$REF,ALT=input$ALT,TYPE=input$TYPE,DataSet=dataSet) #, EVENTLENGTH=input$EVENTLENGTH) 
  if (inPoolCaller) {
    da$goodSite = ((input$FILTER == "PASS" | input$FILTER=="LowQual") & !is.na(input$AC) )
    da$inPoolCaller = inPoolCaller
    da$FILTER=input$FILTER
    da$LOF = input$resource.LOF
    da$NCALLED=input$NCALLED
    da$QUAL=input$QUAL
    names(da)[3] = "AC.PoolCaller"
  }
  else {
#    names(da)[3] = paste("AC.",dataSet,sep="")     
  }
  if (max(input$AC,na.rm=T) > 24 & capNA==T) da$AC = NA
  da$poolNumber=pool
  da
}

rmdupPos <- function(allData) {
  # remove duplicates in first
  merged <- allData
  badPos = subset(count(merged,vars="POS"),freq>1)
  dups <- subset(merged,is.element(POS,badPos$POS))
  merged = merged[!is.element(merged$POS,dups$POS),]
  
}
mergePC <-function(allData, newData,inPoolCaller = F) {
  if (inPoolCaller) {
#    cols <- c("CHROM","POS","REF","ALT","TYPE","DataSet","AC.PoolCaller","FILTER","goodSite","inPoolCaller","LOF","NCALLED","QUAL","poolNumber")
    merged <- merge(allData,newData,all=T)
  } else {
    nd<-data.frame(CHROM=newData$CHROM,POS=newData$POS,AC=newData$AC,REF=newData$REF,ALT=newData$ALT)
    dataSet = unique(newData$DataSet)[1]
    colnames(nd)[3] = paste("AC.",dataSet,sep="")
    merged <-  merge(allData,nd,all=T)
#     # remove duplicates in first
#     merged <- rmdupPos(allData)
#     newData <- rmdupPos(newData)
# 
#     commonPOS = intersect(merged$POS, newData$POS)    
#     indsAll = (merged$POS %in% commonPOS)
#     indsNew = (newData$POS %in% commonPOS)
#     
#     merged$acnew = NA
#     #id <- which(colnames(newData) == "AC")
#     merged$acnew[indsAll==T] = newData$AC[indsNew==T]
#     dataSet = unique(newData$DataSet)[1]
#     
#     idn <- which(colnames(merged) == "acnew")
#     colnames(merged)[idn] = paste("AC.",dataSet,sep="") 
  }
  merged
}

getStats <-function(rinput,pool) {
  inputSNP<-subset(rinput,goodSite==T & TYPE=="SNP" )
  inputINDEL <- subset(rinput,goodSite==T & TYPE=="INDEL")
  da<-data.frame()
  #sensitivity = nrow(subset(input, PoolCallerAC>0 & TruthAC>0))/nrow(subset( input, TruthAC>0))
  #ss = nrow(subset(input, PoolCallerAC>0 & TruthAC==1))/nrow(subset( input, TruthAC==1))
  #sp = nrow(subset(input, PoolCallerAC==0 & TruthAC==0))/nrow(subset( input, TruthAC==0))
  #da<-data.frame(sensitivity=sensitivity,singletonSensitivity=ss, specificity=sp,DataSet=dataSet,
  
  if (nrow(inputSNP) > 0) {
    pcorrS = cor(inputSNP$AC.PoolCaller,inputSNP$AC.OMNI,use="pairwise.complete.obs")
    da<-rbind(da,data.frame(corrValue=pcorrS,TYPE="SNP",DataSet="OMNI"))
    
    axiomS = cor(inputSNP$AC.PoolCaller,inputSNP$AC.Axiom,use="pairwise.complete.obs")
    da<-rbind(da,data.frame(corrValue=axiomS,TYPE="SNP",DataSet="Axiom"))
    eccorrS = cor(inputSNP$AC.PoolCaller,inputSNP$AC.ExomeChip,use="pairwise.complete.obs")
    da<-rbind(da,data.frame(corrValue=eccorrS,TYPE="SNP",DataSet="ExomeChip"))
    ogcorrS = cor(inputSNP$AC.PoolCaller,inputSNP$AC.OneKG,use="pairwise.complete.obs")
    da<-rbind(da,data.frame(corrValue=ogcorrS,TYPE="SNP",DataSet="OneKG"))
    
  }
  
  if (nrow(inputINDEL) > 0) {
#    pcorrI = cor(inputINDEL$AC.PoolCaller,inputINDEL$AC.OMNI,use="pairwise.complete.obs")
#    da<-rbind(da,data.frame(corrValue=pcorrI,TYPE="INDEL",DataSet="OMNI"))
    axiomI = cor(inputINDEL$AC.PoolCaller,inputINDEL$AC.Axiom,use="pairwise.complete.obs")
    da<-rbind(da,data.frame(corrValue=axiomI,TYPE="INDEL",DataSet="Axiom"))
    eccorrI = cor(inputINDEL$AC.PoolCaller,inputINDEL$AC.ExomeChip,use="pairwise.complete.obs")
    da<-rbind(da,data.frame(corrValue=eccorrI,TYPE="INDEL",DataSet="ExomeChip"))
    ogcorrI = cor(inputINDEL$AC.PoolCaller,inputINDEL$AC.OneKG,use="pairwise.complete.obs")
    da<-rbind(da,data.frame(corrValue=ogcorrI,TYPE="INDEL",DataSet="OneKG"))
  }
  da$poolNumber = pool
  da
}

allPoolData<-data.frame()
allPoolStats<-data.frame()

# get global tables first
runName<-"finalRunGGANew"
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.LOFSNP.withRef.filtered.annotated.LOF.table",runName),header=T)
allData<-processData(da,"LOF",pool=0,inPoolCaller=T,capNA=F)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.LOFINDEL.withRef.filtered.annotated.LOF.table",runName),header=T)
ad2<-processData(da,"LOF",pool=0,inPoolCaller=T,capNA=F)
allData<-mergePC(allData,ad2,inPoolCaller=T)

da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.afSNPs.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,dataSet="afSNPs",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.afIndels.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,"afIndels",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.unifSNPs.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,dataSet="unifSNPs",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.unifIndels.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,"unifIndels",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.exomeChip.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,"ExomeChipControl",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.millsPoly.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,"MillsPoly",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.omniPoly.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,"OmniPoly",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.omniMono.withRef.filtered.annotated.table",runName),header=T)
allData<-mergePC(allData,processData(da,"OmniMono",pool=0,inPoolCaller=T,capNA=F),inPoolCaller=T)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/1000G.exomechip.20121009.subsetToGoodPoolSamples.genotypes.table",header=T)
allData<-mergePC(allData,processData(da,"ExomeChip",pool=0,inPoolCaller=F,capNA=F),inPoolCaller=F)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/ALL.wex.axiom.20120206.snps_and_indels.subsetToGoodPoolSamples.genotypes.table",header=T)
allData<-mergePC(allData,processData(da,"Axiom",pool=0,inPoolCaller=F,capNA=F),inPoolCaller=F)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/Omni25_genotypes_2141_samples.b37.samplesInGoodPools.table",header=T)
allData<-mergePC(allData,processData(da,"OMNI",pool=0,inPoolCaller=F,capNA=F),inPoolCaller=F)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/ALL.wgs.phase1_release_v3.20101123.snps_indels_svs.genotypes.InGoodSamples.table",header=T)
allData<-mergePC(allData,processData(da,"OneKG",pool=0,inPoolCaller=F,capNA=F),inPoolCaller=F)

# read na12878-only data to substract from pools 1-10
dna12878<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/NA12878.GGA.single_sample_at_all_validationSites.table",header=T)
da2 <- processData(dna12878,"NA12878",0,inPoolCaller=F,capNA=F)
allData<-mergePC(allData,da2,inPoolCaller=F)


# extra filtering for QUAL
allData$goodSite[allData$QUAL < 100 & allData$AC.PoolCaller > 0] = F
allData <- subset(allData,!is.na(DataSet))
allStats <- getStats(allData,0)
maxPool <- 92

for (pool in c(1:maxPool)) {
  if (pool != 12 && pool != 24 && pool != 72) {
    pstr = sprintf("%d",pool)
    if (pool < 10) {
      pstr = sprintf("0%d",pool)
    }
    # lof pool caller data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.LOFSNP.withRef.filtered.annotated.LOF.atPool%s.table",pstr,runName,pstr),header=T)  
    da2a<-processData(da,dataSet="LOF",pool=pool,inPoolCaller=T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.LOFINDEL.withRef.filtered.annotated.LOF.atPool%s.table",pstr,runName,pstr),header=T)  
    da2b<-processData(da,dataSet="LOF",pool=pool,inPoolCaller=T)
    poolData<-mergePC(da2a,da2b,inPoolCaller=T)

    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.afSNPs.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)  
    poolData<-mergePC(poolData,processData(da,"afSNPs",pool,inPoolCaller=T),inPoolCaller=T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.afIndels.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)  
    poolData<-mergePC(poolData,processData(da,"afIndels",pool,inPoolCaller=T),inPoolCaller=T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.unifSNPs.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)  
    poolData<-mergePC(poolData,processData(da,"unifSNPs",pool,inPoolCaller=T),inPoolCaller=T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.unifIndels.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)  
    poolData<-mergePC(poolData,processData(da,"unifIndels",pool,inPoolCaller=T),inPoolCaller=T)
    # control omni
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.millsPoly.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)  
    poolData<-mergePC(poolData,processData(da,"MillsPoly",pool,inPoolCaller=T),inPoolCaller=T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.omniPoly.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)  
    poolData<-mergePC(poolData,processData(da,"OmniPoly",pool,inPoolCaller=T),inPoolCaller=T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.omniMono.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)  
    poolData<-mergePC(poolData,processData(da,"OmniMono",pool,inPoolCaller=T),inPoolCaller=T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.exomeChip.withRef.filtered.annotated.atPool%s.table",pstr,runName,pstr),header=T)
    poolData<-mergePC(poolData,processData(da,"ExomeChipControl",pool,T),inPoolCaller=T)
    
    #omni data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/Omni25_genotypes_2141_samples.b37.samplesInGoodPools.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"OMNI",pool,inPoolCaller=F)
    poolData<-mergePC(poolData,da2,inPoolCaller=F)

    #axiom data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/ALL.wex.axiom.20120206.snps_and_indels.subsetToGoodPoolSamples.genotypes.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"Axiom",pool,inPoolCaller=F)
    poolData<-mergePC(poolData,da2,inPoolCaller=F)
    
    #exome chip data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/1000G.exomechip.20121009.subsetToGoodPoolSamples.genotypes.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"ExomeChip",pool,inPoolCaller=F)
    poolData<-mergePC(poolData,da2,inPoolCaller=F)
    
    # 1000g data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/ALL.wgs.phase1_release_v3.20101123.snps_indels_svs.genotypes.InGoodSamples.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"OneKG",pool,inPoolCaller=F)
    poolData<-mergePC(poolData,da2,inPoolCaller=F)
    
    # NA12878 data
    da2 <- processData(subset(dna12878,is.element(POS,poolData$POS)),"NA12878",pool,inPoolCaller=F)
    poolData<-mergePC(poolData,da2,inPoolCaller=F)
    poolData <- subset(poolData,!is.na(DataSet))
    
    
#    badPos = subset(count(poolData,vars="POS"),freq>1)
#    dups <- subset(poolData,is.element(POS,badPos$POS))
#    poolData$DUPLICATE = F
#    poolData[is.element(poolData$POS,dups$POS),]$DUPLICATE = T
    
    if (pool <= 10) {
      poolData$AC.PoolCaller = poolData$AC.PoolCaller - poolData$AC.NA12878      
    }

    poolData$goodSite[poolData$QUAL < 100 & poolData$AC.PoolCaller > 0] = F
    
    allPoolData<-rbind(allPoolData,poolData)
    poolStats <- getStats(allPoolData,pool)
    allPoolStats <- rbind(allPoolStats, poolStats)
    cat(pstr)
  }
  
}

#pdf("corr_per_pool.pdf",width=11,height=8.5)
ggplot(allPoolStats,aes(x=poolNumber,y=corrValue,color=DataSet))+geom_point() + facet_grid(TYPE~.)+opts(title="Per-pool Cross-correlation of estimated AC counts vs. chip genotype data or 1000 G genotypes")
#dev.off()

if (writeData) {
  runName = paste(runName,"Qual100",sep="")
  outputDir="/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/processedTables"
  write.table(x=allPoolStats,file=sprintf("%s/%s.allPoolStats.table",outputDir,runName),sep="\t",row.names=F,quote=F)
  write.table(x=allPoolData,file=sprintf("%s/%s.allPoolData.table",outputDir,runName),sep="\t",row.names=F,quote=F)
  write.table(x=allStats,file=sprintf("%s/%s.allStats.table",outputDir,runName),sep="\t",row.names=F,quote=F)
  write.table(x=allData,file=sprintf("%s/%s.allData.table",outputDir,runName),sep="\t",row.names=F,quote=F)
  
}
# 
# snpMatrix = matrix(nrow=nrow(subset(allPoolData,TYPE=="SNP" & poolNumber == 2 & goodSite == T)),ncol=92)
# indelMatrix = matrix(nrow=nrow(subset(allPoolData,TYPE=="INDEL" & poolNumber == 2 & goodSite == T)),ncol=92)
# for (pool  in c(1:92)) {
#   if (pool != 12 && pool != 24 && pool != 72) {
#     snpMatrix[,pool] = subset(allPoolData,TYPE=="SNP" & poolNumber == pool & goodSite == T)$AC.PoolCaller
#     indelMatrix[,pool] = subset(allPoolData,TYPE=="INDEL" & poolNumber == pool & goodSite == T)$AC.PoolCaller
#   }
#   else {
#     snpMatrix[,pool] = -1
#     indelMatrix[,pool] = -1
#   }
#   cat(pool)
# }

# png("PCAC_vs_omni_lof_by_chrom.png",width=640,height=480)
# ggplot(allPoolData,aes(x=TruthAC,y=PoolCallerAC,color=CHROM))+geom_point(alpha=0.3)+facet_grid(refSample~DataSet)+geom_jitter()
# dev.off()
# 
# png("sens_vs_omni_lof_by_chrom.png",width=640,height=480)
# ggplot(subset(allPoolStats,refSample==T),aes(x=poolNumber,y=sensitivity,color=DataSet))+geom_point(size=3,alpha=0.5)
# dev.off()
# 
# png("singletons_sens_vs_omni_lof_by_chrom.png",width=640,height=480)
# ggplot(subset(allPoolStats,refSample==T),aes(x=poolNumber,y=singletonSensitivity,color=DataSet))+geom_point(size=3,alpha=0.5)
# dev.off()
# 
# png("singletons_spec_vs_omni_lof_by_chrom.png",width=640,height=480)
# ggplot(subset(allPoolStats,refSample==T),aes(x=poolNumber,y=specificity,color=DataSet))+geom_point(size=3,alpha=0.5)
# dev.off()
# 
# totalTruthSensOmniRef = nrow(subset(allPoolData,refSample==T & DataSet=="OMNI" & TruthAC>0 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==T & DataSet=="OMNI" & TruthAC>0))
# totalTruthSensOmniNoRef = nrow(subset(allPoolData,refSample==F & DataSet=="OMNI" & TruthAC>0 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==F & DataSet=="OMNI" & TruthAC>0))
# totalTruthSensLOFRef = nrow(subset(allPoolData,refSample==T & DataSet=="LOF" & TruthAC>0 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==T & DataSet=="LOF" & TruthAC>0))
# totalTruthSensLOFNoRef = nrow(subset(allPoolData,refSample==F & DataSet=="LOF" & TruthAC>0 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==F & DataSet=="LOF" & TruthAC>0))
# 
# totalSingletonTruthSensOmniRef = nrow(subset(allPoolData,refSample==T & DataSet=="OMNI" & TruthAC==1 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==T & DataSet=="OMNI" & TruthAC==1))
# totalSingletonTruthSensOmniNoRef = nrow(subset(allPoolData,refSample==F & DataSet=="OMNI" & TruthAC==1 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==F & DataSet=="OMNI" & TruthAC==1))
# totalSingletonTruthSensLOFRef = nrow(subset(allPoolData,refSample==T & DataSet=="LOF" & TruthAC==1 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==T & DataSet=="LOF" & TruthAC==1))
# totalSingletonTruthSensLOFNoRef = nrow(subset(allPoolData,refSample==F & DataSet=="LOF" & TruthAC==1 & PoolCallerAC > 0))/nrow(subset(allPoolData,refSample==F & DataSet=="LOF" & TruthAC==1))
# 
# cat(totalTruthSensOmniRef)
# cat(totalTruthSensOmniNoRef)
# cat(totalTruthSensLOFRef)
# cat(totalTruthSensLOFNoRef)
# cat(totalSingletonTruthSensOmniRef)
# cat(totalSingletonTruthSensOmniNoRef)
# cat(totalSingletonTruthSensLOFRef)
# cat(totalSingletonTruthSensLOFNoRef)
