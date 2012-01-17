#!/bin/env Rscript

require("lattice")
require("ggplot2")
require("splines")

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
  d$technology <- factor(1, levels=c("HiSeq-paper", "GA2-1000G", "HiSeq-recent"))
  d$technology[grepl("ERR.*", d$rg)] <- "GA2-1000G"
  d$technology[grepl("20.*", d$rg)] <- "HiSeq-paper"
  d$technology[grepl("B00EG.*", d$rg)] <- "HiSeq-recent"
  print(summary(d$technology))

  eByComp = addEmpiricalPofG(ddply(d, .(rg, technology, pGGivenDType, pGGivenD), genotypeCounts))
  return(list(d=d, eByComp = eByComp))
  #countsByTech = addEmpiricalPofG(ddply(d, .(technology, pGGivenDType, pGGivenD), genotypeCounts))
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

ymax = xmax = 30
goodEByComp = subset(eByComp, Sum > 10 & EmpiricalPofGQ < Inf)

print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, size=log10(Sum), facets = pGGivenDType ~ technology, color=pGGivenDType, geom=c("point", "smooth"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2))
print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ technology, color=rg, geom=c("blank"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2) + geom_smooth(se=F, aes(weight=Sum)))
print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ technology, color=rg, geom=c("blank"), group=rg, xlim=c(0,xmax), ylim=c(-10,10)) + geom_abline(slope=0, linetype=2) + geom_smooth(se=F, method=lm, formula = y ~ ns(x,1), aes(weight=Sum)))
print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., color=technology, geom=c("blank"), group=technology, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2) + geom_smooth(se=T, size=1.5, aes(weight=Sum)))