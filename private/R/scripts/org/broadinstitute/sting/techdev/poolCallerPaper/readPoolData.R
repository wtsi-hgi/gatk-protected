 
processData <-function(input,dataSet,pool,doFilter,capNA=T) {
  da<-data.frame(CHROM=input$CHROM,POS=input$POS,AC=input$AC,REF=input$REF,ALT=input$ALT,TYPE=input$TYPE, EVENTLENGTH=input$EVENTLENGTH) 
  if (doFilter) {
    da$goodSite = ((input$FILTER == "PASS" | input$FILTER=="LowQual") & !is.na(input$AC) )
    da$inLOF = !is.na(input$AC)
    da$FILTER=input$FILTER
    da$LOF = input$resource.LOF
    da$NCALLED=input$NCALLED
    da$QUAL=input$QUAL
  }
  
  if (max(input$AC,na.rm=T) > 24 & capNA==T) da$AC = NA
  names(da)[3] = paste("AC.",dataSet,sep="") 
  da$poolNumber=pool
  da
}

getStats <-function(rinput,pool) {
  inputSNP<-subset(rinput,goodSite==T & TYPE=="SNP" )
  inputINDEL <- subset(rinput,goodSite==T & TYPE=="INDEL")
  da<-data.frame()
  #sensitivity = nrow(subset(input, PoolCallerAC>0 & TruthAC>0))/nrow(subset( input, TruthAC>0))
  #ss = nrow(subset(input, PoolCallerAC>0 & TruthAC==1))/nrow(subset( input, TruthAC==1))
  #sp = nrow(subset(input, PoolCallerAC==0 & TruthAC==0))/nrow(subset( input, TruthAC==0))
  #da<-data.frame(sensitivity=sensitivity,singletonSensitivity=ss, specificity=sp,DataSet=dataSet,
  pcorrS = cor(inputSNP$AC.LOF,inputSNP$AC.OMNI,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=pcorrS,TYPE="SNP",DataSet="OMNI"))
  
  axiomS = cor(inputSNP$AC.LOF,inputSNP$AC.Axiom,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=axiomS,TYPE="SNP",DataSet="Axiom"))
  eccorrS = cor(inputSNP$AC.LOF,inputSNP$AC.ExomeChip,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=eccorrS,TYPE="SNP",DataSet="ExomeChip"))
  ogcorrS = cor(inputSNP$AC.LOF,inputSNP$AC.OneKG,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=ogcorrS,TYPE="SNP",DataSet="OneKG"))
  
  pcorrI = cor(inputINDEL$AC.LOF,inputINDEL$AC.OMNI,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=pcorrI,TYPE="INDEL",DataSet="OMNI"))
  axiomI = cor(inputINDEL$AC.LOF,inputINDEL$AC.Axiom,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=axiomI,TYPE="INDEL",DataSet="Axiom"))
  eccorrI = cor(inputINDEL$AC.LOF,inputINDEL$AC.ExomeChip,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=eccorrI,TYPE="INDEL",DataSet="ExomeChip"))
  ogcorrI = cor(inputINDEL$AC.LOF,inputINDEL$AC.OneKG,use="pairwise.complete.obs")
  da<-rbind(da,data.frame(corrValue=ogcorrI,TYPE="INDEL",DataSet="OneKG"))
  da$poolNumber = pool
  da
}

allPoolData<-data.frame()
allPoolStats<-data.frame()

# get global tables first
runName<-"finalRunGGA"
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.LOFSNP.withRef.filtered.annotated.LOF.table",runName),header=T)
allData<-processData(da,"LOF",0,T,F)
da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/%s.LOFINDEL.withRef.filtered.annotated.LOF.table",runName),header=T)
ad2<-processData(da,"LOF",0,T,F)
allData<-merge(allData,ad2,all=T)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/1000G.exomechip.20121009.subsetToGoodPoolSamples.genotypes.table",header=T)
allData<-merge(allData,processData(da,"ExomeChip",0,F,F),all=T)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/ALL.wex.axiom.20120206.snps_and_indels.subsetToGoodPoolSamples.genotypes.table",header=T)
allData<-merge(allData,processData(da,"Axiom",0,F,F),all=T)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/Omni25_genotypes_2141_samples.b37.samplesInGoodPools.table",header=T)
allData<-merge(allData,processData(da,"OMNI",0,F,F),all=T)

da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/ALL.wgs.phase1_release_v3.20101123.snps_indels_svs.genotypes.InGoodSamples.table",header=T)
allData<-merge(allData,processData(da,"OneKG",0,F,F),all=T)

# read na12878-only data to substract from pools 1-10
dna12878<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/NA12878.GGA.indels.singleSample.table",header=T)
#da$AC[da$NCALLED==1 &is.na(da$AC) ]=0
#da$AC[da$QUAL<0] = "LowCoverage"
#dm<-subset(da,is.element(POS,allData$POS))
da2 <- processData(dna12878,"NA12878",0,doFilter=F,capNA=F)
da<-read.table("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/NA12878.GGA.snps.singleSample.table",header=T)
dna12878=merge(dna12878,da,all=T)
#da$AC[da$NCALLED==1 &is.na(da$AC) ]=0
#da$AC[da$QUAL<0] = "LowCoverage"
#dm<-subset(da,is.element(POS,allData$POS))
da2 <- processData(dna12878,"NA12878",0,doFilter=F,capNA=F)
allData<-merge(allData,da2,all=T)

# mark duplicate positions
badPos = subset(count(allData,vars="POS"),freq>1)
dups <- subset(allData,is.element(POS,badPos$POS))
allData$DUPLICATE = F
allData[is.element(allData$POS,dups$POS),]$DUPLICATE = T



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
    da2a<-processData(da,"LOF",pool,T)
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/%s.LOFINDEL.withRef.filtered.annotated.LOF.atPool%s.table",pstr,runName,pstr),header=T)  
    da2b<-processData(da,"LOF",pool,T)
    poolData<-merge(da2a,da2b,suffixes=c(".SNP",".INDEL"),all=T)

    # control omni
    # TODO
    #control exomechip
    #TODO
    #omni data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/Omni25_genotypes_2141_samples.b37.samplesInGoodPools.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"OMNI",pool,F)
    poolData<-merge(poolData,da2,all=T)

    #axiom data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/ALL.wex.axiom.20120206.snps_and_indels.subsetToGoodPoolSamples.genotypes.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"Axiom",pool,F)
    poolData<-merge(poolData,da2,all=T)
    
    #exome chip data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/1000G.exomechip.20121009.subsetToGoodPoolSamples.genotypes.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"ExomeChip",pool,F)
    poolData<-merge(poolData,da2,all=T)
    
    # 1000g data
    da<-read.table(sprintf("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/allPoolData/poolData%s/ALL.wgs.phase1_release_v3.20101123.snps_indels_svs.genotypes.InGoodSamples.atPool%s.table",pstr,pstr),header=T)  
    da2<-processData(da,"OneKG",pool,F)
    poolData<-merge(poolData,da2,all=T)
    
    # NA12878 data
    da2 <- processData(subset(dna12878,is.element(POS,poolData$POS)),"NA12878",pool,F)
    poolData<-merge(poolData,da2,all=T)
    
    
    badPos = subset(count(poolData,vars="POS"),freq>1)
    dups <- subset(poolData,is.element(POS,badPos$POS))
    poolData$DUPLICATE = F
    poolData[is.element(poolData$POS,dups$POS),]$DUPLICATE = T
    
    allPoolData<-rbind(allPoolData,poolData)
    poolStats <- getStats(allPoolData,pool)
    allPoolStats <- rbind(allPoolStats, poolStats)
    cat(pstr)
  }
  
}

pdf("corr_per_pool.pdf",width=11,height=8.5)
ggplot(allPoolStats,aes(x=poolNumber,y=corrValue,color=DataSet))+geom_point() + facet_grid(TYPE~.)+opts(title="Per-pool Cross-correlation of estimated AC counts vs. chip genotype data or 1000 G genotypes")
dev.off()

outputDir="/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/processedTables"
write.table(x=allPoolStats,file=sprintf("%s/%s.allPoolStats.table",outputDir,runName),sep="\t",row.names=F,quote=F)
write.table(x=allPoolData,file=sprintf("%s/%s.allPoolData.table",outputDir,runName),sep="\t",row.names=F,quote=F)
write.table(x=allStats,file=sprintf("%s/%s.allStats.table",outputDir,runName),sep="\t",row.names=F,quote=F)
write.table(x=allData,file=sprintf("%s/%s.allData.table",outputDir,runName),sep="\t",row.names=F,quote=F)

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
