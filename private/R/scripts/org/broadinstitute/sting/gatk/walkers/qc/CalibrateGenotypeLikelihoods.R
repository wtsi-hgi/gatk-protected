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
    r = c(r, v / row$Sum)
  }

  #print(length(r))
  d$EmpiricalPofG = r
  d$EmpiricalPofGQ = round(-10*log10(1-r))
  return(d)
}

genotypeCounts <- function(x) {
  type = unique(x$variable)[1]
  t = addmargins(table(x$comp))
  return(t)
}

digestTable <- function(inputDataFile) {
  d = subset(read.table(inputDataFile, header=T), rg != "ALL")
  eByComp = addEmpiricalPofG(ddply(d, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
  
  return(list(d=d, eByComp = eByComp))
}

##########################################
### The script
##########################################

args <- commandArgs(TRUE)
inputDataFile = args[1]
onCmdLine = ! is.na(inputDataFile)
if ( onCmdLine ) {
  eByComp <- digestTable(inputDataFile)$eByComp
}
pdf(paste(inputDataFile, ".pdf", sep=""))

for ( xmax in c(30, 99) ) { # loop over just meaningful subset and all 
ymax = 30
goodEByComp = subset(eByComp, Sum > 10 & EmpiricalPofGQ < Inf)

#First graph, overall likelihoods 
tryCatch(
  print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=pGGivenDType, geom=c("jitter", "smooth"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2)),
  error = function(e) {
    print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=pGGivenDType, geom=c("jitter"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2))
    }
)  

#Second graph, likelihoods by read group  
tryCatch(
  print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2) + geom_smooth(se=F, aes(weight=Sum))),
  error = function(e) {
    print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2))    
  }
)

#Third graph, likelihoods difference from empirical    
tryCatch(
  print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-10,10)) + geom_abline(slope=0, linetype=2) + geom_smooth(se=F, method=lm, formula = y ~ ns(x,1), aes(weight=Sum))),
  error = function(e) {
    print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., size=log10(Sum), color=rg, geom=c("jitter"), group=rg, xlim=c(0,xmax), ylim=c(-10,10)) + geom_abline(slope=0, linetype=2))
  }
)
}
