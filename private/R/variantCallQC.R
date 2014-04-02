library(scales)
library(gplots)
library(gsalib)
library(ggplot2)
library(reshape)
library(tools)


# Utilities for displaying multiple plots per page
# -------------------------------------------------------

distributeGraphRows <- function(graphs, heights = c()) {
  # Viewport layout 2 graphs top to bottom with given relative heights
  #
  #
  if (length(heights) == 0) {
    heights <- rep.int(1, length(graphs))
  }
  heights <- heights[!is.na(graphs)]
  graphs <- graphs[!is.na(graphs)]
  numGraphs <- length(graphs)
  Layout <- grid.layout(nrow = numGraphs, ncol = 1, heights=heights)
  grid.newpage()
  pushViewport(viewport(layout = Layout))
  subplot <- function(x) viewport(layout.pos.row = x, layout.pos.col = 1)
  for (i in 1:numGraphs) {
    print(graphs[[i]], vp = subplot(i))
  }
}

distributeLogGraph <- function(graph, xName) {
  continuousGraph <- graph + scale_x_continuous(xName)
  logGraph <- graph + scale_x_log10(xName) + opts(title="")
  distributeGraphRows(list(continuousGraph, logGraph))
}

distributePerSampleGraph <- function(perSampleGraph, distGraph, ratio=c(2,1)) {
  distributeGraphRows(list(perSampleGraph, distGraph), ratio)
}

removeExtraStrats <- function(variantEvalDataFrame, moreToRemove=c()) {
  # Remove the standard extra stratification columns FunctionalClass, Novelty, and others in moreToRemove from the variantEvalDataFrame
  #
  # Only keeps the column marked with "all" for each removed column
  #
  for ( toRemove in c("FunctionalClass", "Novelty", moreToRemove) ) {
    if (toRemove %in% colnames(variantEvalDataFrame)) {
      variantEvalDataFrame <- variantEvalDataFrame[variantEvalDataFrame[[toRemove]] == "all",]
    }
  }
  variantEvalDataFrame    
}

openPDF <- function(outputPDF) {
  # Open the outputPDF file with standard dimensions, if outputPDF is not NA
  if ( ! is.na(outputPDF) ) {
    pdf(outputPDF, height=8.5, width=11)
  }
}

closePDF <- function(outputPDF) {
  # close the outputPDF file if not NA, and try to compact the PDF if possible
  if ( ! is.na(outputPDF) ) {
    dev.off()
    if (exists("compactPDF")) {
      compactPDF(outputPDF)
    }
  }
}

makeRatioDataFrame <- function(ACs, num, denom, widths = NULL) {
  if ( is.null(widths) ) widths <- rep(1, length(ACs))
  
  value = NULL
  titv <- data.frame(AC=ACs, width = widths, num=num, denom = denom, ratio = num / denom)
}

.reduceACs <- function(binWidthForAC, ACs) {
  # computes data structures necessary to reduce the full range of ACs
  #
  # binWidthForAC returns the number of upcoming bins that should be merged into 
  # that AC bin.  ACs is a vector of all AC values from 0 to 2N that should be 
  # merged together
  #
  # Returns a list containing the reduced ACs starts, their corresponding widths,
  # and a map from original ACs to their new ones (1 -> 1, 2 -> 2, 3 -> 2, etc)
  maxAC <- max(ACs)
  newACs <- c()
  widths <- c()
  newACMap <- c()
  ac <- 0
  while ( ac < maxAC ) {
    newACs <- c(newACs, ac)
    width <- binWidthForAC(ac)
    widths <- c(widths, width)
    newACMap <- c(newACMap, rep(ac, width))
    ac <- ac + width
  }
  list(ACs = newACs, widths=widths, newACMap = newACMap)
}

# geometricACs <- function(k, ACs) {
#   nBins <- round(k * log10(max(ACs)))
#   
#   binWidthForAC <- function(ac) {
#     max(ceiling(ac / nBins), 1)
#   }
#   
#   return(reduceACs(binWidthForAC, ACs))
# }

reduce.AC.on.LogLinear.intervals <- function(scaleFactor, ACs) {
  # map the full range of AC values onto a log linear scale
  #
  # Reduce the full AC range onto one where the width of each new AC increases at a rate of
  # 10^scaleFactor in size with growing AC values.  This is primarily useful for accurately
  # computing ratios or other quantities by AC that aren't well determined when the AC 
  # values are very large
  #
  # Returns a list containing the reduced ACs starts, their corresponding widths,
  # and a map from original ACs to their new ones (1 -> 1, 2 -> 2, 3 -> 2, etc)
  maxAC <- max(ACs)
  afs <- ACs / maxAC
  breaks <- 10^(seq(-4, -1, scaleFactor))
  widths <- c()
  lastBreak <- 1
  for ( i in length(breaks):1 ) {
    b <- breaks[i]
    width <- sum(afs < lastBreak & afs >= b)
    widths <- c(widths, width)
    lastBreak <- b
  }
  widths <- rev(widths)
  
  binWidthForAC <- function(ac) {
    af <- ac / maxAC
    value = 1
    for ( i in length(breaks):1 )
      if ( af >= breaks[i] ) {
        value = widths[i]
        break
      }
    
    return(value)
  }
  
  return(.reduceACs(binWidthForAC, ACs))
}

.remapACs <- function(remapper, k, df) {
  newACs <- remapper(k, df$AC)
  
  n = length(newACs$ACs)
  num = rep(0, n)
  denom = rep(0, n)
  for ( i in 1:dim(df)[1] ) {
    rowI = df$AC == i
    row = df[rowI,]
    newAC = newACs$newACMap[row$AC]
    newRowI = newACs$ACs == newAC
    num[newRowI] = num[newRowI] + df$num[rowI]
    denom[newRowI] = denom[newRowI] + df$denom[rowI]
  }
  
  newdf <- makeRatioDataFrame(newACs$ACs, num, denom, newACs$widths )
  newdf
}

compute.ratio.on.LogLinear.AC.intervals <- function(ACs, num, denom, scaleFactor = 0.1) {
  df = makeRatioDataFrame(ACs, num, denom, 1)
  return(.remapACs(reduce.AC.on.LogLinear.intervals, scaleFactor, df))
}

plotVariantQC <- function(metrics, measures, requestedStrat = "Sample", 
                          fixHistogramX=F, anotherStrat = NULL, nObsField = "n_indels", 
                          onSamePage=F, facetVariableOnXPerSample = F, facetVariableOnXForDist = T, 
                          moreTitle="", note = NULL) {
  metrics$strat = metrics[[requestedStrat]]
  
  otherFacet = "."
  id.vars = c("strat", "nobs")
  metrics$nobs <- metrics[[nObsField]]
  
  # keep track of the other strat and it's implied facet value
  if (! is.null(anotherStrat)) { 
    id.vars = c(id.vars, anotherStrat)
    otherFacet = anotherStrat
  }
  
  molten <- melt(metrics, id.vars=id.vars, measure.vars=c(measures))
  perSampleGraph <- ggplot(data=molten, aes(x=strat, y=value, group=variable, color=variable, fill=variable))

  # create the title
  titleText=paste(paste(paste(measures, collapse=", "), "by", requestedStrat), moreTitle)
  if ( !is.null(note) ) {
    titleText=paste(titleText, note, sep="\n")
  }
  paste(titleText)
  title <- opts(title=titleText)
  
  determineFacet <- function(onX) {
    if ( onX ) { 
      paste(otherFacet, "~ variable")
    } else {
      paste("variable ~", otherFacet)
    }
  }
  
  sampleFacet = determineFacet(facetVariableOnXPerSample)
  distFacet   = determineFacet(facetVariableOnXForDist)
  
  if ( requestedStrat == "Sample" ) {
    perSampleGraph <- perSampleGraph + geom_text(aes(label=strat), size=1.5) + geom_blank() # don't display a scale
    perSampleGraph <- perSampleGraph + scale_x_discrete("Sample (ordered by nSNPs)", labels=function(x) "")
  } else { # by AlleleCount
    perSampleGraph <- perSampleGraph + geom_point(aes(size=log10(nobs))) #+ geom_smooth(aes(weight=log10(nobs)))
    perSampleGraph <- perSampleGraph + scale_x_log10("AlleleCount")
  }    
  perSampleGraph <- perSampleGraph + ylab("Variable value") + title
  perSampleGraph <- perSampleGraph + facet_grid(sampleFacet, scales="free")
  
  nValues = length(unique(molten$value))
  if (nValues > 2) {
    if ( requestedStrat == "Sample" ) {
      distGraph <- ggplot(data=molten, aes(x=value, group=variable, fill=variable))
    } else {
      distGraph <- ggplot(data=molten, aes(x=value, group=variable, fill=variable, weight=nobs))
    }
    distGraph <- distGraph + geom_histogram(aes(y=..ndensity..))
    distGraph <- distGraph + geom_density(alpha=0.5, aes(y=..scaled..))
    distGraph <- distGraph + geom_rug(aes(y=NULL, color=variable, position="jitter"))
    scale = "free"
    if ( fixHistogramX ) scale = "fixed"
    distGraph <- distGraph + facet_grid(distFacet, scales=scale)
    distGraph <- distGraph + ylab("Relative frequency")
    distGraph <- distGraph + xlab("Variable value (see facet for variable by color)")
    distGraph <- distGraph + opts(axis.text.x=theme_text(angle=-45)) # , legend.position="none")
  } else {
    distGraph <- NA
  }
  
  if ( onSamePage ) {
    suppressMessages(distributePerSampleGraph(perSampleGraph, distGraph))
  } else {
    suppressMessages(print(perSampleGraph))
    suppressMessages(print(distGraph + title))
  }
}


# for testing only
#source("~/Desktop/broadLocal/GATK/unstable/public/R/src/org/broadinstitute/sting/utils/R/gsalib/R/gsa.variantqc.utils.R")

# TODOs:
#  Assumes you have indels in your call set.  If not you will get errors

args <- commandArgs(TRUE)

onCMDLine <- ! is.na(args[1])
LOAD_DATA <- T

if ( onCMDLine ) {
  projectName <- args[1]
  bySampleEval <- args[2]
  byACEval <- args[3]
  IndelQCEval <- args[4]
  outputPDF <- args[5]
} else {
  projectName <- "InDevelopmentInR"

  #root <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/indelQC/C783_277_826_calling8Mar2012"
  root <- "~/Desktop/broadLocal/tmp.data/esp.all.snps.unfiltered.indels.filtered.bcf"
  #root <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/indelQC/esp.all.unannotated.chr1"
  #root <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/indelQC/ALL.wex.broad.illumina.20110521.snps.indels.genotypes"
  #root <- "/humgen/gsa-hpprojects/ESP/calls/broadOnly_chr1_v2/esp.all.unannotated.chr1"
  
  bySampleEval <- paste(root, "bySample.eval", sep=".")
  byACEval <- paste(root, "byAC.eval", sep=".")
  IndelQCEval <- paste(root, "indelQC.eval", sep=".")
  
  # comment out to stop printing to PDF
  outputPDF <- NA
  #outputPDF <- "variantCallQC.pdf"
}

theme_set(theme_bw())

print("Report")
print(paste("Project          :", projectName))
print(paste("bySampleEval     :", bySampleEval))
print(paste("byACEval         :", byACEval))
print(paste("IndelQC          :", IndelQCEval))
print(paste("outputPDF        :", outputPDF))

expandVEReport <- function(d) {
  d$TiTvVariantEvaluator$tiTvRatio <- round(d$TiTvVariantEvaluator$tiTvRatio,2) 
  d$CountVariants$insertionDeletionRatio <- round(d$CountVariants$insertionDeletionRatio,2) 
  d$CountVariants$nIndels <- d$CountVariants$nInsertions + d$CountVariants$nDeletions
  return(d)
}

# Filters out any metric that isnt cumulative for the specified column name
selectCumulativeMetrics <- function(data,column) {
  lapply(data,function(x) { subset(x,x[column]=='all') })
}

createMetricsBySites <- function(bySampleReport, byACReport) {
  # Metrics by sites:
  #  bySample -> counts of SNPs and Indels by novelty, with expectations
  #  byAC -> snps and indels (known / novel)
  bySample <- expandVEReport(selectCumulativeMetrics(bySampleReport,'Sample'))
  byAC <- byACReport
  r <- list(bySample = bySample,byAC = byAC)
  r$byAC$CountVariants$nIndels <- r$byAC$CountVariants$nInsertions + r$byAC$CountVariants$nDeletions
  r$byAC$TiTvVariantEvaluator$nSNPs <- r$byAC$TiTvVariantEvaluator$nTi + r$byAC$TiTvVariantEvaluator$nTv
  r$byAC$CountVariants$AC <- r$byAC$CountVariants$AlleleCount
  r$byAC$TiTvVariantEvaluator$AC <- r$byAC$TiTvVariantEvaluator$AlleleCount
  return(r)
}

summaryTable <- function(metricsBySites, metricsBySample) {
  # SNP summary statistics
  merged <- merge(metricsBySites$bySample$CountVariants, metricsBySites$bySample$TiTvVariantEvaluator)
  sub <- subset(merged, FunctionalClass=="all")
  raw <- melt(sub, id.vars=c("Novelty"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "insertionDeletionRatio"))

    # Create a data.frame with a column for missense/silent ratio by Novelty
  countMissenseSilent <- subset(metricsBySites$bySample$CountVariants, FunctionalClass %in% c("silent", "missense"))
  noveltyMissenseSilent <- recast(countMissenseSilent, Novelty ~ FunctionalClass, id.var=c("Novelty","FunctionalClass"), measure.var=c("nCalledLoci"))
  noveltyMissenseSilent$missenseSilentRatio <- round(noveltyMissenseSilent$missense / noveltyMissenseSilent$silent, 2)

  raw <- rbind(data.frame(raw), melt.data.frame(noveltyMissenseSilent, id.vars=c("Novelty"), measure.vars=c("missenseSilentRatio")))
  table <- cast(raw, Novelty ~ ...)
  # doesn't work with textplot
  colnames(table) <- c("Novelty", "Target Size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "Ins/Del Ratio", "Ka/Ks")
  return(table)
}

sampleSummaryTable <- function(metricsBySamples, missenseSilentSummary) {
  # SNP summary statistics
  raw <- melt(metricsBySamples, id.vars=c("Novelty", "Sample"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "insertionDeletionRatio"))
  table <- cast(rbind(raw, missenseSilentSummary), Novelty ~ variable, mean)
  table$nSNPs <- round(table$nSNPs, 0)
  table$nIndels <- round(table$nIndels, 0)
  table$tiTvRatio <- round(table$tiTvRatio, 2)
  table$insertionDeletionRatio <- round(table$insertionDeletionRatio, 2)
  table$missenseSilentRatio <- round(table$missenseSilentRatio, 2)
  colnames(table) <- c("Novelty", "Target Size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "Ins/Del Ratio", "Ka/Ks")
  return(table)
}

overallSummaryTable <- function(metricsBySites, metricsBySamples, missenseSilentSummary) {
  sitesSummary <- as.data.frame(summaryTable(metricsBySites, metricsBySamples))
  sitesSummary$"Metric Type" <- "Sites"
  sampleSummary <- as.data.frame(sampleSummaryTable(metricsBySamples, missenseSilentSummary))
  sampleSummary$"Metric Type" <- "Per-sample avg."
  # that last item puts the metric.type second in the list
  table <- rbind(sitesSummary, sampleSummary)[, c(1,8,2,3,4,7,5,6)]
  # remove columns with all NA/NaN
  table <- table <- table[,colSums(is.na(table))<nrow(table)]
  return(table)
}

summaryPlots <- function(metricsBySites) {
  # See note about ggplot with less than two points in perSamplePlots.
  numAC <- max(metricsBySites$byAC$TiTvVariantEvaluator$AC)
  countVariants <- subset(metricsBySites$byAC$CountVariants, FunctionalClass == "all")
  tiTvVariantEvaluator <- subset(metricsBySites$byAC$TiTvVariantEvaluator, FunctionalClass == "all")

  name <- "SNP and Indel count by novelty and allele count"
  molten <- melt(subset(countVariants, Novelty != "all" & AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(c("nSNPs", "nIndels")))
  acGraph <- ggplot(data=molten, aes(x=AC, y=value+1, color=Novelty, fill=Novelty), group=variable)
  acGraph <- acGraph + opts(title = name)
  acGraph <- acGraph + scale_y_log10("Number of variants")
  acGraph <- acGraph + geom_point(alpha=0.5, size=3)
  acGraph <- acGraph + geom_line(size=1)
  acGraph <- acGraph + facet_grid(variable ~ ., scales="free")
  distributeLogGraph(acGraph, "Allele count (AC)")
  
  name <- "Transition / transversion ratio by allele count"
  # nSNPs > 0 => requires that we have some data here, otherwise Ti/Tv is zero from VE  
  minSNPsToInclude <- 0
  if (sum(tiTvVariantEvaluator$nSNPs) > 0) {
    byACNoAll <- subset(tiTvVariantEvaluator, Novelty != "all" & AC > 0 & nSNPs > minSNPsToInclude)
    byACNoAll$weightSNPs <- byACNoAll$nSNPs
    byACNoAll$weightSNPs[byACNoAll$weightSNPs < 10] <- 0 # no weight given to bins with small counts
    acGraph <- ggplot(data=byACNoAll, aes(x=AC, y=(nTi+1)/(nTv+1), color=Novelty))
    acGraph <- acGraph + scale_y_continuous("Transition / transversion ratio", limits=c(0,4))
    acGraph <- acGraph + opts(title = name)
    acGraph <- acGraph + geom_point(aes(size=log10(nSNPs), weight=weightSNPs), alpha=0.5) + geom_smooth()
    distributeLogGraph(acGraph, "Allele count (AC)")
  }
  
  countFunctional <- subset(metricsBySites$byAC$CountVariants, FunctionalClass != "all" & AlleleCount > 0)

  name <- "SNP counts by functional class" 
  if (sum(countVariants$nSNPs) > 0) {
    molten <- melt(subset(countFunctional, Novelty != "all"), id.vars=c("Novelty", "FunctionalClass"), measure.vars=c(c("nSNPs")))
    if ( sum(molten$value) > 0 ) {
      p <- ggplot(data=cast(molten, Novelty + FunctionalClass ~ ..., sum), aes(x=factor(FunctionalClass), y=nSNPs, fill=Novelty), group=FunctionalClass)
      p <- p + opts(title = name)
      p <- p + scale_y_log10("No. of SNPs")
      p <- p + geom_bar(position="dodge")
      p <- p + xlab("FunctionalClass")
      print(p)

      countNoNonsense = subset(countFunctional, Novelty == "all")
      p <- ggplot(data=countNoNonsense, aes(x=AlleleCount, y=nCalledLoci, color=factor(FunctionalClass)))
      p <- p + opts(title = "Functional Class by allele count")
      p <- p + geom_point(alpha=0.5, size=3)
      p <- p + geom_line(size=1)
      p <- p + scale_y_log10("No. of SNPs")
      distributeLogGraph(p, "Allele count (AC)")
    }
  }
}

addSection <- function(name) {
    print(paste("Running section", name))
    par("mar", c(5, 4, 4, 2))
    frame()
    title(name, cex=2)
}
 
# -------------------------------------------------------
# read functions
# -------------------------------------------------------

createMetricsBySamples <- function(bySampleReport) {
  bySample <- expandVEReport(selectCumulativeMetrics(bySampleReport,'FunctionalClass'))
  r <- merge(bySample$TiTvVariantEvaluator, bySample$CountVariants)
  r <- merge(r, bySample$CompOverlap)
  # order the samples by nSNPs -- it's the natural ordering.
  x <- subset(r, Novelty=="all")
  r$Sample <- factor(x$Sample, levels=x$Sample[order(x$nSNPs)])

  return(subset(r, Sample != "all"))
}

createMissenseSilentSummary <- function(bySampleReport) {
  raw <- melt(subset(bySampleReport$CountVariants, FunctionalClass %in% c("silent", "missense") & Sample != "all"), id.vars=c("FunctionalClass", "Novelty", "Sample"), measure.vars=c("nVariantLoci"))
  table <- data.frame(cast(raw, Novelty + Sample ~ FunctionalClass))
  table$missenseSilentRatio <- table$missense / table$silent
  return(melt(table, id.vars=c("Novelty", "Sample"), measure.vars=c("missenseSilentRatio")))
}

# -------------------------------------------------------
# Per sample plots
# -------------------------------------------------------

perSamplePlots <- function(metricsBySamples) {
  # TODO - add MNPs
  #measuresSets <- list(c("nSNPs", "nIndels"), c("tiTvRatio", "insertionDeletionRatio")) # for testing
  measuresSets <- list(c("nSNPs", "nIndels", "nSingletons"), c("tiTvRatio", "insertionDeletionRatio", "hetHomRatio"))
  for ( measures in measuresSets ) {
    plotVariantQC(metricsBySamples, measures, requestedStrat = "Sample", fixHistogramX=F, anotherStrat = "Novelty", nObsField = "nSNPs", 
                  onSamePage=T, facetVariableOnXPerSample=F, facetVariableOnXForDist=T)
    for ( measure in measures ) {
      plotVariantQC(metricsBySamples, measure, requestedStrat = "Sample", fixHistogramX=F, anotherStrat = "Novelty", nObsField = "nSNPs", 
                    onSamePage=T, facetVariableOnXPerSample=T, facetVariableOnXForDist=F)
    }
      
    # known / novel ratio by sample
    # TODO -- would ideally not conflate SNPs and Indels
    d <- subset(metricsBySamples, Novelty == "all" & CompRod == "dbsnp")
    noveltyRateName <- "Percent of variants in dbSNP"
    d[[noveltyRateName]] <- d$compRate
    plotVariantQC(d, c(noveltyRateName), requestedStrat = "Sample", fixHistogramX=F, nObsField = "nSNPs", onSamePage=T)
  
    for ( novelty in c("all", "known", "novel") ) {
      subd = subset(metricsBySamples, Novelty==novelty)
      plotVariantQC(subd, measures, requestedStrat = "Sample", fixHistogramX=F, nObsField = "nSNPs", onSamePage=T, 
                    moreTitle=paste(novelty, "sites"), facetVariableOnXPerSample=F, facetVariableOnXForDist=T)
    }
  }
}

# -------------------------------------------------------
# Detailed indel QC statistics 
# -------------------------------------------------------

plotRatioByAlleleCount <- function(AC, num, denom, name, expectedValue, keepFullAC=F) {
  df <- data.frame(AlleleCount = AC, nobs = denom, num = num, denom = denom, ratio = num / denom)
  df <- subset(df, num > 0 & denom > 0)
  
  maxAC = max(AC)
  note = NULL
  if ( maxAC > 100 ) {
    highAC = maxAC * 0.1 # all values with MAF > 10%
    highACValues = subset(df, AlleleCount > highAC)
    highACRatio = sum(highACValues$num) / sum(highACValues$denom)
    note=paste("Note ratios at high MAF with many samples may appear artificially low.  Ratio for MAF > 10% is", sprintf("%.2f", highACRatio))
  }    
  
  df[[name]] <- df$ratio
  plotVariantQC(df, c(name), requestedStrat = "AlleleCount", nObsField="nobs", note = note)
}

indelLengthDistribution <- function(indelHistogram) {
  indelHistogram$callset <- "per sample"
  numSamples = length(unique(indelHistogram$Sample))
  indelHistogram$callset[indelHistogram$Sample == "all"] = "overall"
  indelHistogram$callset = factor(indelHistogram$callset, levels=c("per sample", "overall"), ordered=T)
  p <- ggplot(data=subset(indelHistogram, Sample != "all"), aes(x=Length, y=Freq, group=interaction(Length, TandemRepeat), fill=TandemRepeat))
  p <- p + geom_vline(x=c(-9,-6,-3,3,6,9), linetype="dashed", color="grey",size=2)
  p <- p + geom_boxplot()
  p <- p + scale_x_continuous(breaks=unique(sort(indelHistogram$Length)))
  p <- p + xlab("Indel length (negative is deletion)") + ylab("Relative frequency")

  p2 <- p + facet_grid(TandemRepeat ~ .)
  #p2 <- p2 + geom_point(position="jitter")
  p2 <- p2 + geom_line(aes(group=Sample, color=TandemRepeat), data=subset(indelHistogram, Sample == "all"), size=2)
  
  print(p)
  print(p2)
  #distributePerSampleGraph(p, p2, ratio=c(1,1))
}
#indelLengthDistribution(lengthHistogram)

indelPlots <- function(IndelQCReport, byACReport) {
  IndelSummaryBySample <- IndelQCReport$IndelSummary
  IndelSummaryBySampleWithoutAll = subset(IndelSummaryBySample, Sample != "all")
  IndelSummaryBySampleNoOtherStrats <- removeExtraStrats(IndelSummaryBySampleWithoutAll, c("OneBPIndel", "TandemRepeat"))
  IndelSummaryByAC = removeExtraStrats(subset(byACReport$IndelSummary, AlleleCount > 0)) 
  
  # write out the values for all IndelSummary report fields transposed for easy viewing
  all = subset(removeExtraStrats(IndelSummaryBySample, c("OneBPIndel", "TandemRepeat")), Sample == "all")
  all = all[, !(colnames(all) %in% c("IndelSummary", "Sample", "CompRod","EvalRod") )]
  rownames(all) <- "Combined callset"
  textplot(t(all))
  
  plotVariantQCWithAllStrats <- function(values, ...) {
    plotVariantQC(IndelSummaryBySampleNoOtherStrats, values, ...)
    plotVariantQC(removeExtraStrats(IndelSummaryBySampleWithoutAll, c("OneBPIndel")), values, anotherStrat="TandemRepeat", ...)
    plotVariantQC(removeExtraStrats(IndelSummaryBySampleWithoutAll, c("TandemRepeat")), values, anotherStrat="OneBPIndel",...)
  }

  plotVariantQCWithTandemStrat <- function(values, ...) {
    plotVariantQC(IndelSummaryBySampleNoOtherStrats, values, ...)
    plotVariantQC(removeExtraStrats(IndelSummaryBySampleWithoutAll, c("OneBPIndel")), values, anotherStrat="TandemRepeat",...)
  }
  
  if ( F ) {
    plotVariantQC(IndelSummaryBySampleNoOtherStrats, c("n_SNPs", "n_indels", "SNP_to_indel_ratio"))
    plotVariantQC(IndelSummaryBySampleNoOtherStrats, c("n_singleton_SNPs", "n_singleton_indels", "SNP_to_indel_ratio_for_singletons"))
    plotVariantQCWithAllStrats(c("n_indels", "n_singleton_indels"))
    plotVariantQCWithAllStrats(c("n_indels_matching_gold_standard", "gold_standard_matching_rate"))
    
    if ( all$n_multiallelic_indel_sites > 0 ) { # if there are at least some multi-allelic sites
      plotVariantQCWithAllStrats(c("n_multiallelic_indel_sites", "percent_of_sites_with_more_than_2_alleles"))
    }
	
    plotVariantQCWithAllStrats(c("indel_novelty_rate"))
    
    # insertion to deletion information
    plotVariantQCWithAllStrats(c("insertion_to_deletion_ratio"), fixHistogramX=T)
    
    # special handling of these ratios, as OneBP strat doesn't make sense
    plotVariantQCWithTandemStrat(c("ratio_of_1_and_2_to_3_bp_insertions", "ratio_of_1_and_2_to_3_bp_deletions"), fixHistogramX=T)
	
    # het : hom ratios information
    plotVariantQCWithAllStrats(c("SNP_het_to_hom_ratio", "indel_het_to_hom_ratio"), fixHistogramX=T)
    
    plotRatioByAlleleCount(IndelSummaryByAC$AlleleCount, IndelSummaryByAC$n_insertions, IndelSummaryByAC$n_deletions, "Insertion to deletion ratio", 1)
    plotRatioByAlleleCount(IndelSummaryByAC$AlleleCount, IndelSummaryByAC$n_novel_indels*100, IndelSummaryByAC$n_indels, "% novel indels", 50)
    plotRatioByAlleleCount(IndelSummaryByAC$AlleleCount, IndelSummaryByAC$n_SNPs, IndelSummaryByAC$n_indels, "SNP to indel ratio", 50)
    # optional frameshift counts
    hasFrameShift = ! is.na(max(IndelSummaryBySampleNoOtherStrats$frameshift_rate_for_coding_indels))
    if ( hasFrameShift ) {
      plotVariantQC(IndelSummaryBySampleNoOtherStrats, c("frameshift_rate_for_coding_indels"))
      plotRatioByAlleleCount(IndelSummaryByAC$AlleleCount, IndelSummaryByAC$n_coding_indels_frameshifting * 100, IndelSummaryByAC$n_coding_indels_frameshifting + IndelSummaryByAC$n_coding_indels_in_frame, "Indel frameshift rate (%)", 0.1)
    } 
    lengthHistogram <- removeExtraStrats(IndelQCReport$IndelLengthHistogram, c("OneBPIndel"))
    indelLengthDistribution(lengthHistogram)
   } else {
      # for testing only
     plotRatioByAlleleCount(IndelSummaryByAC$AlleleCount, IndelSummaryByAC$n_insertions, IndelSummaryByAC$n_deletions, "Insertion to deletion ratio", 1)
   }
}

# -------------------------------------------------------
# Actually invoke the above plotting functions 
# -------------------------------------------------------

# load the data.
if ( onCMDLine || LOAD_DATA ) {
  bySampleReport <- gsa.read.gatkreport(bySampleEval)
  byACReport <- gsa.read.gatkreport(byACEval)
  IndelQCReport <- gsa.read.gatkreport(IndelQCEval)
  metricsBySites <- createMetricsBySites(bySampleReport, byACReport)
  metricsBySamples <- createMetricsBySamples(bySampleReport)
  missenseSilentSummary <- createMissenseSilentSummary(bySampleReport)
}

pdf(outputPDF, height=8.5, width=11)

# Table of overall counts and quality
textplot(overallSummaryTable(metricsBySites, metricsBySamples, missenseSilentSummary), show.rownames=F)
title(paste("Summary metrics for project", projectName), cex=3)

addSection("Summary plots")
summaryPlots(metricsBySites)
addSection("Properties per sample")
perSamplePlots(metricsBySamples)
addSection("Detailed indel QC metrics")
indelPlots(IndelQCReport, metricsBySites$byAC)

dev.off()
compactPDF(outputPDF)
