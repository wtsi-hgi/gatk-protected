library(gplots)
library(gsalib)
library(ggplot2)
#library(tools)

# TODOs:
#  Assumes you have indels in your call set.  If not you will get errors

args <- commandArgs(TRUE)

onCMDLine <- ! is.na(args[1])
LOAD_DATA <- T

# creates an array of c(sampleName1, ..., sampleNameN)
parseHighlightSamples <- function(s) {
  return(unlist(strsplit(s, ",", fixed=T)))
}

if ( onCMDLine ) {
  projectName <- args[1]
  bySampleEval <- args[2]
  byACEval <- args[3]
  IndelQCEval <- args[4]
  outputPDF <- args[5]
  if ( ! is.na(args[6]) )
    highlightSamples <- parseHighlightSamples(args[6])
  else
    highlightSamples <- c()
} else {
  projectName <- "InDevelopmentInR"

  root <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/indelQC/C783_277_826_calling8Mar2012"
  #root <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/indelQC/esp.all.unannotated.chr1"
  #root <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/indelQC/ALL.wex.broad.illumina.20110521.snps.indels.genotypes"
  #root <- "/humgen/gsa-hpprojects/ESP/calls/broadOnly_chr1_v2/esp.all.unannotated.chr1"
  
  bySampleEval <- paste(root, "bySample.eval", sep=".")
  byACEval <- paste(root, "byAC.eval", sep=".")
  IndelQCEval <- paste(root, "indelQC.eval", sep=".")
  
  # comment out to stop printing to PDF
  outputPDF <- NA
  #outputPDF <- "variantCallQC.pdf"
  highlightSamples <- c() # parseHighlightSamples("29029,47243")
}

theme_set(theme_bw())

print("Report")
print(paste("Project          :", projectName))
print(paste("bySampleEval     :", bySampleEval))
print(paste("byACEval         :", byACEval))
print(paste("IndelQC          :", IndelQCEval))
print(paste("outputPDF        :", outputPDF))
print(paste("highlightSamples :", highlightSamples))

# -------------------------------------------------------
# Utilities for displaying multiple plots per page
# -------------------------------------------------------

# Viewport (layout 2 graphs top to bottom)
distributeGraphRows <- function(graphs, heights = c()) {
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
  colnames(table) <- c("Novelty", "Target Size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "Del/Ins Ratio", "Ka/Ks")
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
  colnames(table) <- c("Novelty", "Target Size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "Del/Ins Ratio", "Ka/Ks")
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

  # Counts vs. Allele frequency 
  name <- "Variant counts by allele count"
  for ( measure in c("nSNPs", "nIndels")) {
    molten <- melt(subset(countVariants, AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(measure))
    if ( sum(molten$value > 0) ) {
      p <- ggplot(data=molten, aes(x=AC, y=value+1, color=Novelty), group=variable)
      p <- p + opts(title = paste(name, ":", measure))
      p <- p + scale_y_log10("Number of variants")
      p <- p + scale_x_log10("Allele count (AC)")
      if (numAC > 2) {
        p <- p + geom_smooth(aes(weight=value), size=1, method="lm", formula = y ~ x)
      } else {
        p <- p + geom_line(size=1)
      }
      p <- p + geom_point(alpha=0.5, size=4)
      p <- p + facet_grid(Novelty ~ ., scales="free")
      print(p)
    }
  }
  
  name <- "Transition / transversion ratio by allele count"
  # nSNPs > 0 => requires that we have some data here, otherwise Ti/Tv is zero from VE  
  minSNPsToInclude <- 0
  if (sum(tiTvVariantEvaluator$nSNPs) > 0) {
    byACNoAll <- subset(tiTvVariantEvaluator, Novelty != "all" & AC > 0 & nSNPs > minSNPsToInclude)
    acGraph <- ggplot(data=byACNoAll, aes(x=AC, y=tiTvRatio, color=Novelty))
    acGraph <- acGraph + scale_y_continuous("Transition / transversion ratio", limits=c(0,4))
    acGraph <- acGraph + opts(title = name)
    if (numAC > 2) {
      acGraph <- acGraph + geom_smooth(size=2)
    } else {
      acGraph <- acGraph + geom_line(size=2)
    }
    acGraph <- acGraph + geom_point(aes(size=log10(nSNPs), weight=nSNPs), alpha=0.5)
    distributeLogGraph(acGraph, "Allele count (AC)")
  }
 
  # SNPs to indels ratio by allele frequency
  name <- "SNPs to indels ratio by allele count"
  if ( sum(countVariants$nIndels) > 0 & sum(countVariants$nSNPs) > 0 ) {
    countVariants$SNP.Indel.Ratio <- countVariants$nSNPs / countVariants$nIndels
    countVariants$SNP.Indel.Ratio[countVariants$nIndels == 0] <- NaN
    p <- ggplot(data=subset(countVariants, Novelty == "all" & nSNPs > 0), aes(x=AC, y=SNP.Indel.Ratio))
    p <- p + opts(title = name)
    p <- p + scale_y_continuous("SNP to indel ratio")
    if (numAC > 2) {
      p <- p + geom_smooth(size=2, aes(weight=nIndels))
    } else {
      p <- p + geom_line(size=2)
    }
    p <- p + geom_point(alpha=0.5, aes(size=log10(nIndels)))
    distributeLogGraph(p, "Allele count (AC)")
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

  # add highlight info
  r$highlight <- r$Sample %in% highlightSamples

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
  # Some ggplot function do not work with two or less samples, for example
  # PASS: ggplot(data=data.frame(Grouping=c('a','a','a','c','c','c'),Values=1:6),aes(x=Values,group=Grouping)) + geom_density()
  # FAIL: ggplot(data=data.frame(Grouping=c('a','a','b','b','c','c'),Values=1:6),aes(x=Values,group=Grouping)) + geom_density()
  #   Error: attempt to apply non-function
  numSamples <- length(unique(metricsBySamples$Sample))
  metricsBySamples$highlightTextSizes <- c(1,2)[metricsBySamples$highlight+1]
  if (TRUE %in% metricsBySamples$highlight) {
    # TODO: How do you scale the relative aes and the geom_text at the same time?
    # In R 2.13 the font size is much bigger and not sure the syntax to size it correctly
    sampleTextLabel <- geom_text(aes(label=Sample, size=highlightTextSizes))
    sampleTextLabelScale <- scale_size("Highlighted samples", to=c(3,5), breaks=c(1,2), labels=c("regular", "highlighted"))
  } else {
    # Nothing to highlight
    # If the size was in the aes then R 2.13 wasn't scaling the text
    # https://groups.google.com/d/msg/ggplot2/uj38mwCyY0Q/dOHhLyLEW8YJ
    sampleTextLabel <- geom_text(aes(label=Sample), size=1.5)
    sampleTextLabelScale <- geom_blank() # don't display a scale
  }
  xAxis <- scale_x_discrete("Sample (ordered by nSNPs)", formatter=function(x) "")

  measures <- c("nSNPs", "tiTvRatio", "nSingletons", "nIndels", "insertionDeletionRatio")
  name <- "by sample"
  for ( measure in measures ) {
    molten <- melt(metricsBySamples, id.vars=c("Novelty", "Sample", "highlightTextSizes"), measure.vars=c(measure))

    perSampleGraph <- ggplot(data=molten, aes(x=Sample, y=value, group=Novelty, color=Novelty), y=value)
    perSampleGraph <- perSampleGraph + opts(title = paste(measure, name))
    if (numSamples > 2) {
      perSampleGraph <- perSampleGraph + geom_smooth(alpha=0.5, aes(group=Novelty))
    } else {
      perSampleGraph <- perSampleGraph + geom_line(alpha=0.5, size=1)
    }
    perSampleGraph <- perSampleGraph + sampleTextLabel + sampleTextLabelScale
    perSampleGraph <- perSampleGraph + facet_grid(Novelty ~ ., scales="free")
    perSampleGraph <- perSampleGraph + xAxis
    
    if (numSamples > 2) {
      distGraph <- ggplot(data=molten, aes(x=value, group=Novelty, fill=Novelty))
      #distGraph <- distGraph + geom_density(alpha=0.5)
      distGraph <- distGraph + geom_histogram(aes(y=..ndensity..))
      distGraph <- distGraph + geom_density(alpha=0.5, aes(y=..scaled..))
      distGraph <- distGraph + geom_rug(aes(y=NULL, color=Novelty, position="jitter"))
      distGraph <- distGraph + facet_grid(. ~ Novelty, scales="free")
      distGraph <- distGraph + ylab("Relative frequency")
      distGraph <- distGraph + scale_x_continuous(measure)
    } else {
      distGraph <- NA
    }

    distributePerSampleGraph(perSampleGraph, distGraph)
  }
    
  # known / novel ratio by sample
  # TODO -- would ideally not conflate SNPs and Indels
  d <- subset(metricsBySamples, Novelty == "all" & CompRod == "dbsnp")
  title <- opts(title = "Novelty rate by sample")

  perSampleGraph <- ggplot(data=d, aes(x=Sample, y=compRate))
  perSampleGraph <- perSampleGraph + title
  if (numSamples > 2) {
    perSampleGraph <- perSampleGraph + geom_smooth(alpha=0.5, aes(group=Novelty))
  } else {
    perSampleGraph <- perSampleGraph + geom_line(alpha=0.5, size=1)
  }
  perSampleGraph <- perSampleGraph + sampleTextLabel + sampleTextLabelScale
  perSampleGraph <- perSampleGraph + geom_rug(aes(x=NULL, position="jitter"))
  perSampleGraph <- perSampleGraph + xAxis
  perSampleGraph <- perSampleGraph + scale_y_continuous("Percent of variants in dbSNP")

  if (numSamples > 2) {
    distGraph <- ggplot(data=d, aes(x=compRate))
    #distGraph <- distGraph + geom_density(alpha=0.5)
    distGraph <- distGraph + geom_histogram(aes(y=..ndensity..))
    distGraph <- distGraph + geom_density(alpha=0.5, aes(y=..scaled..))
    distGraph <- distGraph + geom_rug(aes(y=NULL, position="jitter"))
    distGraph <- distGraph + ylab("Relative frequency")
    distGraph <- distGraph + scale_x_continuous("Percent of variants in dbSNP")
  } else {
    distGraph <- NA
  }

  distributePerSampleGraph(perSampleGraph, distGraph)

  for ( novelty in c("all", "known", "novel") ) {
    # TODO -- how can I color it as before?
    # TODO -- add marginal distributions?
    molten <- melt(subset(metricsBySamples, Novelty==novelty), id.vars=c("Sample", "highlightTextSizes"), measure.vars=measures)
    p <- ggplot(data=molten, aes(x=Sample, y=value))
    p <- p + opts(title = paste(name, ":", novelty))
    p <- p + sampleTextLabel + sampleTextLabelScale
    p <- p + facet_grid(variable ~ ., scales="free")
    # how do we remove the labels?
    p <- p + xAxis
    print(p)
  }
}

# -------------------------------------------------------
# Detailed indel QC statistics 
# -------------------------------------------------------

removeExtraStrats <- function(df, moreToRemove=c()) {
  for ( toRemove in c("FunctionalClass", "Novelty", moreToRemove) ) {
    if (toRemove %in% colnames(df)) {
      df <- df[df[[toRemove]] == "all",]
    }
  }
  df    
}

indelQCPlot <- function(metrics, measures, requestedStrat = "Sample", fixHistogramX=F, anotherStrat = NULL) {
  numSamples = dim(metrics)[1] - 1
  metrics$strat = metrics[[requestedStrat]]
  
  otherFacet = "."
  id.vars = c("strat", "n_indels")
  
  # keep track of the other strat and it's implied facet value
  if (! is.null(anotherStrat)) { 
    id.vars = c(id.vars, anotherStrat)
    otherFacet = anotherStrat
  }
  
  molten <- melt(metrics, id.vars=id.vars, measure.vars=c(measures))
  perSampleGraph <- ggplot(data=molten, aes(x=strat, y=value, group=variable, color=variable, fill=variable))
  if ( requestedStrat == "Sample" ) {
    perSampleGraph <- perSampleGraph + geom_text(aes(label=strat), size=1.5) + geom_blank() # don't display a scale
    perSampleGraph <- perSampleGraph + scale_x_discrete("Sample (ordered by nSNPs)", formatter=function(x) "")
  } else {
    perSampleGraph <- perSampleGraph + geom_point(aes(size=log10(n_indels))) + geom_smooth(aes(weight=log10(n_indels)))
    perSampleGraph <- perSampleGraph + scale_x_log10("AlleleCount")
  }    
  perSampleGraph <- perSampleGraph + ylab("Variable value")

  perSampleGraph <- perSampleGraph + facet_grid(paste("variable ~ ", otherFacet), scales="free")

    
  if (numSamples > 2) {
    distGraph <- ggplot(data=molten, aes(x=value, group=variable, fill=variable))
    distGraph <- distGraph + geom_histogram(aes(y=..ndensity..))
    distGraph <- distGraph + geom_density(alpha=0.5, aes(y=..scaled..))
    distGraph <- distGraph + geom_rug(aes(y=NULL, color=variable, position="jitter"))
    scale = "free"
    if ( fixHistogramX ) scale = "fixed"
    distGraph <- distGraph + facet_grid(paste(otherFacet, " ~ variable"), scales=scale)
    distGraph <- distGraph + ylab("Relative frequency")
    distGraph <- distGraph + xlab("Variable value (see facet for variable by color)")
    distGraph <- distGraph + opts(axis.text.x=theme_text(angle=-45)) # , legend.position="none")
  } else {
    distGraph <- NA
  }
  
  print(perSampleGraph)
  print(distGraph)
  #distributePerSampleGraph(perSampleGraph, distGraph)
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

  indelQCPlotWithAllStrats <- function(values, ...) {
    indelQCPlot(IndelSummaryBySampleNoOtherStrats, values, ...)
    indelQCPlot(removeExtraStrats(IndelSummaryBySampleWithoutAll, c("OneBPIndel")), values, anotherStrat="TandemRepeat",...)
    indelQCPlot(removeExtraStrats(IndelSummaryBySampleWithoutAll, c("TandemRepeat")), values, anotherStrat="OneBPIndel",...)
  }
  
  if ( T ) {
    indelQCPlot(IndelSummaryBySampleNoOtherStrats, c("n_SNPs", "n_indels", "SNP_to_indel_ratio"))
    indelQCPlot(IndelSummaryBySampleNoOtherStrats, c("n_singleton_SNPs", "n_singleton_indels", "SNP_to_indel_ratio_for_singletons"))
    indelQCPlotWithAllStrats(c("n_indels", "n_singleton_indels"))
    indelQCPlot(IndelSummaryByAC, c("SNP_to_indel_ratio"), "AlleleCount")
    indelQCPlotWithAllStrats(c("n_indels_matching_gold_standard", "gold_standard_matching_rate"))
    
    if ( all$n_multiallelic_indel_sites > 0 ) { # if there are at least some multi-allelic sites
      indelQCPlotWithAllStrats(c("n_multiallelic_indel_sites", "percent_of_sites_with_more_than_2_alleles"))
    }

    indelQCPlotWithAllStrats(c("indel_novelty_rate"))
    indelQCPlot(IndelSummaryByAC, c("indel_novelty_rate"), "AlleleCount")
    
    # insertion to deletion information
    indelQCPlotWithAllStrats(c("insertion_to_deletion_ratio"), fixHistogramX=T)
    indelQCPlotWithAllStrats(c("ratio_of_1_and_2_to_3_bp_insertions", "ratio_of_1_and_2_to_3_bp_deletions"), fixHistogramX=T)

    # het : hom ratios information
    indelQCPlotWithAllStrats(c("SNP_het_to_hom_ratio", "indel_het_to_hom_ratio"), fixHistogramX=T)

    # optional frameshift counts
    hasFrameShift = ! is.na(max(IndelSummaryBySampleNoOtherStrats$frameshift_rate_for_coding_indels))
    if ( hasFrameShift ) {
      indelQCPlot(IndelSummaryBySampleNoOtherStrats, c("frameshift_rate_for_coding_indels"))
      indelQCPlot(IndelSummaryByAC, c("frameshift_rate_for_coding_indels"), "AlleleCount")
    } 

    lengthHistogram <- removeExtraStrats(IndelQCReport$IndelLengthHistogram, c("OneBPIndel"))
    indelLengthDistribution(lengthHistogram)
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

if ( ! is.na(outputPDF) ) {
  pdf(outputPDF, height=8.5, width=11)
}

# Table of overall counts and quality
textplot(overallSummaryTable(metricsBySites, metricsBySamples, missenseSilentSummary), show.rownames=F)
title(paste("Summary metrics for project", projectName), cex=3)

summaryPlots(metricsBySites)
perSamplePlots(metricsBySamples)
indelPlots(IndelQCReport, metricsBySites$byAC)

if ( ! is.na(outputPDF) ) {
  dev.off()
  if (exists("compactPDF")) {
    compactPDF(outputPDF)
  }
}
