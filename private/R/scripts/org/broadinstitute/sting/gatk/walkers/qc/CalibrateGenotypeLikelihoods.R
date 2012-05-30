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

digestTable <- function(inputDataFile, doCompAll) {
  d = read.table(inputDataFile, header=T)
  byRG = subset(d, rg != "ALL")
  allOnly = subset(d, rg == "ALL")
  eByCompRG = addEmpiricalPofG(ddply(byRG, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
  eByCompRG$rg = factor(eByCompRG$rg)
  if ( doCompAll == 1 ) {
    eByCompAll = addEmpiricalPofG(ddply(allOnly, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
    eByCompAll$rg = factor(eByCompAll$rg)
  } else {
    eByCompAll <- eByCompRG
  }
  
  return(list(d=d, eByCompRG = eByCompRG, eByCompAll = eByCompAll))
}

##########################################
### The script
##########################################

args <- commandArgs(TRUE)
inputDataFile = args[1]
onCmdLine = ! is.na(inputDataFile)
if ( onCmdLine ) {
  doCompAll <- as.numeric(args[2])
  digested <- digestTable(inputDataFile, doCompAll)
} else {
  doCompAll <- 1
  digested <- digestTable("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/QuantizedQuals/levels.64.qq.reduced.cgl", doCompAll)
}
if ( onCmdLine ) pdf(paste(inputDataFile, ".pdf", sep=""))

plotMe <- function(eByComp, includeByReadGroup, title) {
  for ( xmax in c(30, 99) ) { # loop over just meaningful subset and all 
    ymax = xmax
    goodEByComp = subset(eByComp, Sum > 10 & EmpiricalPofGQ < Inf)
    
    #First graph, overall likelihoods 
    tryCatch(
      print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=pGGivenDType, geom=c("jitter"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_smooth(se=T, aes(weight=Sum)) + geom_abline(slope=1, linetype=2), title=title),
      error = function(e) {
        print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=pGGivenDType, geom=c("jitter"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2), title=title)
        }
    )  
    
    if ( includeByReadGroup ) {
      #Second graph, likelihoods by read group  
      tryCatch(
        print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2) + geom_smooth(se=F, aes(weight=Sum)), title=title),
        error = function(e) {
          print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2), title=title)    
        }
      )
    }
    
    #Third graph, likelihoods difference from empirical    
    tryCatch(
      print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-ymax/2,ymax/2)) + geom_abline(slope=0, linetype=2) + geom_smooth(se=F, method=lm, formula = y ~ ns(x,1), aes(weight=Sum)), title=title),
      error = function(e) {
        print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-ymax/2,ymax/2)) + geom_abline(slope=0, linetype=2), title=title)
      }
    )
  }
}
  
plotMe(digested$eByCompRG, F, "GLs within read groups")
if ( doCompAll == 1 ) {
    plotMe(digested$eByCompAll, F, "GLs across read groups")
}

if ( onCmdLine ) dev.off()
  