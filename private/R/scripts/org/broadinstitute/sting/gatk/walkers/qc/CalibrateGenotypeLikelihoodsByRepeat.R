#!/bin/env Rscript

library("lattice")
library("ggplot2")
library("splines")

##########################################
### Accessory functions
##########################################

addEmpiricalPofG <- function(d) {
  r = c()
  #
  # TODO -- this is a really naive estimate of the accuracy, as it assumes the comp
  # track is perfect.  In reality the chip is at best Q30 accurate (replicate samples have
  # level than this level of concordance).  At low incoming confidence, we can effectively
  # ignore this term but when the incoming Q is near or above Q30 this approximation clearly
  # breaks down.
  #
  for ( i in 1:dim(d)[1] ) {
    row = d[i,]
    if ( row$pGGivenDType == "QofAAGivenD" ) v = row$HOM_REF
    if ( row$pGGivenDType == "QofABGivenD" ) v = row$HET
    if ( row$pGGivenDType == "QofBBGivenD" ) v = row$HOM_VAR
    ratio = (v+1) / (row$Sum+2)
    print(list(v=v, sum = row$Sum, ratio=ratio))
    r = c(r, ratio)
  }

  #print(length(r))
  d$EmpiricalPofG = r
  d$EmpiricalPofGQ = round(sapply(-10*log10(1-d$EmpiricalPofG), function(x) min(x, 93)))
  return(d)
}

genotypeCounts <- function(x) {
  type = unique(x$variable)[1]
  t = addmargins(table(x$comp))
  return(t)
}

rmseFunction <- function(data) {
  rmse = sqrt(sum((data$pGGivenD-data$EmpiricalPofGQ)^2*data$Sum)/sum(data$Sum))
  c(eqn=paste("RMSE=",round(rmse,digits=3),sep=""), Sum=mean(data$Sum))
}

digestTable <- function(inputDataFile, doCompAll, doRepeats) {
  d = read.table(inputDataFile, header=T)
  allOnlyr = subset(d, rg == "ALL" & isRepeat == "true")
  byRG = subset(d, rg != "ALL")
  byRGr = subset(byRG, isRepeat == "true")
  eByCompRGr = addEmpiricalPofG(ddply(byRGr, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
  eByCompRGr$rg = factor(eByCompRGr$rg)
  eByCompRGr$isRepeat = "true"
  if ( doCompAll == 1 ) {
    eByCompAllr = addEmpiricalPofG(ddply(allOnlyr, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
    eByCompAllr$rg = factor(eByCompAllr$rg)
    eByCompAllr$isRepeat= "true"
  } else {
    eByCompAllr <- eByCompRGr
  }
  
  if (doRepeats == 0) {
    return(list(d=d, eByCompRG = eByCompRGr, eByCompAll = eByCompAllr))
    
  } else {
    allOnlyn = subset(d, rg == "ALL" & isRepeat == "false")
    byRGn = subset(byRG, isRepeat == "false")
    eByCompRGn = addEmpiricalPofG(ddply(byRGn, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
    eByCompRGn$rg = factor(eByCompRGn$rg)
    eByCompRGn$isRepeat = "false"
    if ( doCompAll == 1 ) {
      eByCompAlln = addEmpiricalPofG(ddply(allOnlyn, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
      eByCompAlln$rg = factor(eByCompAlln$rg)
      eByCompAlln$isRepeat= "false"
    } else {
      eByCompAlln <- eByCompRGn
    }
    
    return(list(d=d, eByCompRG = rbind(eByCompRGr,eByCompRGn), eByCompAll = rbind(eByCompAllr,eByCompRGn)))
    
  }
}

##########################################
### The script
##########################################
require(ggplot2)
require(plyr)
args <- commandArgs(TRUE)
inputDataFile = args[1]
onCmdLine = ! is.na(inputDataFile)
if ( onCmdLine ) {
  doCompAll <- as.numeric(args[2])
  doRepeats <- as.numeric(args[3])
  digested <- digestTable(inputDataFile, doCompAll, doRepeats)
} else {
  doCompAll <- 0
  doRepeats <- 1
  digested <- digestTable("/humgen/gsa-scr1/delangel/GATK/Sting_unstable_ping/du6rg.out", doCompAll, doRepeats)
}
if ( onCmdLine ) pdf(paste(inputDataFile, ".pdf", sep=""))

plotMe <- function(eByComp, includeByReadGroup, title, doRepeats) {
  for ( xmax in c(30, 40) ) { # loop over just meaningful subset and all 
    ymax = xmax
    if(xmax > 70) {
      labelPlacement = xmax-10
    } else {
      labelPlacement = xmax-5
    }
    goodEByComp = subset(eByComp, Sum > 10 & EmpiricalPofGQ < Inf & pGGivenD <= xmax)
    rmseLabel <- ddply(goodEByComp, .(pGGivenDType), rmseFunction)
    rmseLabel$Sum = as.numeric(max(rmseLabel$Sum))

    if (doRepeats == 1) {
      #First graph, overall likelihoods 
      tryCatch(
        print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ isRepeat, size=log10(Sum), color=pGGivenDType, geom=c("jitter"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_smooth(se=T, aes(weight=Sum)) + geom_abline(slope=1, linetype=2) + geom_text(data=rmseLabel, aes(x=0, y=labelPlacement, label=eqn, hjust=0, vjust=0), fontface="bold"), title=title),
        error = function(e) {
          print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ isRepeat, size=log10(Sum), color=pGGivenDType, geom=c("jitter"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_abline(slope=1, linetype=2) + geom_text(data=rmseLabel, aes(x=0, y=labelPlacement, label=eqn, hjust=0, vjust=0), fontface="bold"), title=title)
        }
      )  
      
      if ( includeByReadGroup ) {
        #Second graph, likelihoods by read group  
        tryCatch(
          print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ isRepeat, size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_abline(slope=1, linetype=2) + geom_smooth(se=F, aes(weight=Sum)), title=title),
          error = function(e) {
            print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ isRepeat, size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_abline(slope=1, linetype=2), title=title)
          }
        )
      }
      
      #Third graph, likelihoods difference from empirical    
      tryCatch(
        print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ isRepeat, size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-ymax/2,ymax/2), xlab="Reported Genotype Quality", ylab="Genotype Quality Accuracy (reported - empirical)") + geom_abline(slope=0, linetype=2) + geom_smooth(se=F, method=lm, formula = y ~ ns(x,1), aes(weight=Sum)), title=title),
        error = function(e) {
          print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ isRepeat, size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-ymax/2,ymax/2), xlab="Reported Genotype Quality", ylab="Genotype Quality Accuracy (reported - empirical)") + geom_abline(slope=0, linetype=2), title=title)
        }
      )
      
    }
    else {
      #First graph, overall likelihoods 
      tryCatch(
        print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=pGGivenDType, geom=c("jitter"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_smooth(se=T, aes(weight=Sum)) + geom_abline(slope=1, linetype=2) + geom_text(data=rmseLabel, aes(x=0, y=labelPlacement, label=eqn, hjust=0, vjust=0), fontface="bold"), title=title),
        error = function(e) {
          print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=pGGivenDType, geom=c("jitter"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_abline(slope=1, linetype=2) + geom_text(data=rmseLabel, aes(x=0, y=labelPlacement, label=eqn, hjust=0, vjust=0), fontface="bold"), title=title)
        }
      )  
      
      if ( includeByReadGroup ) {
        #Second graph, likelihoods by read group  
        tryCatch(
          print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_abline(slope=1, linetype=2) + geom_smooth(se=F, aes(weight=Sum)), title=title),
          error = function(e) {
            print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Reported Genotype Quality", ylab="Empirical Genotype Quality") + geom_abline(slope=1, linetype=2), title=title)
          }
        )
      }
      
      #Third graph, likelihoods difference from empirical    
      tryCatch(
        print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-ymax/2,ymax/2), xlab="Reported Genotype Quality", ylab="Genotype Quality Accuracy (reported - empirical)") + geom_abline(slope=0, linetype=2) + geom_smooth(se=F, method=lm, formula = y ~ ns(x,1), aes(weight=Sum)), title=title),
        error = function(e) {
          print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-ymax/2,ymax/2), xlab="Reported Genotype Quality", ylab="Genotype Quality Accuracy (reported - empirical)") + geom_abline(slope=0, linetype=2), title=title)
        }
      )
      
    }
  }
}
  
plotMe(digested$eByCompRG, T, "GLs within read groups", doRepeats)
if ( doCompAll == 1 ) {
    plotMe(digested$eByCompAll, F, "GLs across read groups", doRepeats)
}

if ( onCmdLine ) dev.off()
  
