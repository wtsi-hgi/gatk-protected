library(gsalib)
library(ggplot2)
library(gplots)
library(tools)

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
  outputPDF <- args[4]
  if ( ! is.na(args[5]) )
    highlightSamples <- parseHighlightSamples(args[5])
  else
    highlightSamples <- c()
} else {
  projectName <- "InDevelopmentInR"
  bySampleEval <- "/humgen/gsa-hpprojects/dev/kshakir/scratch/postqc/1kg_example.bySample.eval"
  byACEval <- "/humgen/gsa-hpprojects/dev/kshakir/scratch/postqc/1kg_example.byAC.eval"
  outputPDF <- "1kg_example.test.pdf"
  highlightSamples <- c() # parseHighlightSamples("29029,47243")
}

print("Report")
print(paste("Project          :", projectName))
print(paste("bySampleEval     :", bySampleEval))
print(paste("byACEval         :", byACEval))
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

distributePerSampleGraph <- function(perSampleGraph, distGraph) {
  distributeGraphRows(list(perSampleGraph, distGraph), c(2,1))
}

expandVEReport <- function(d) {
  d$TiTvVariantEvaluator$tiTvRatio <- round(d$TiTvVariantEvaluator$tiTvRatio,2) 
  d$CountVariants$deletionInsertionRatio <- round(d$CountVariants$deletionInsertionRatio,2) 
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
  raw <- melt(sub, id.vars=c("Novelty"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "deletionInsertionRatio"))

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
  raw <- melt(metricsBySamples, id.vars=c("Novelty", "Sample"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "deletionInsertionRatio"))
  table <- cast(rbind(raw, missenseSilentSummary), Novelty ~ variable, mean)
  table$nSNPs <- round(table$nSNPs, 0)
  table$nIndels <- round(table$nIndels, 0)
  table$tiTvRatio <- round(table$tiTvRatio, 2)
  table$deletionInsertionRatio <- round(table$deletionInsertionRatio, 2)
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
    print(p)
  }
  
  countFunctional <- subset(metricsBySites$byAC$CountVariants, FunctionalClass != "all" & AlleleCount > 0)

  name <- "SNP counts by functional class" 
  if (sum(countVariants$nSNPs) > 0) {
    molten <- melt(subset(countFunctional, Novelty != "all"), id.vars=c("Novelty", "FunctionalClass"), measure.vars=c(c("nSNPs")))
    if ( sum(molten$value) > 0 ) {
      p <- ggplot(data=cast(molten, Novelty + FunctionalClass ~ ..., sum), aes(x=FunctionalClass, y=nSNPs, fill=Novelty), group=FunctionalClass)
      p <- p + opts(title = name)
      p <- p + scale_y_log10("No. of SNPs")
      p <- p + geom_bar(position="dodge")
      print(p)

      countNoNonsense = subset(countFunctional, Novelty == "all" & FunctionalClass != "nonsense")
      p <- ggplot(data=countNoNonsense, aes(x=AlleleCount, y=nCalledLoci, color=FunctionalClass))
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

  measures <- c("nSNPs", "tiTvRatio", "nSingletons", "nIndels", "deletionInsertionRatio")
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
      distGraph <- distGraph + geom_density(alpha=0.5)
      distGraph <- distGraph + geom_rug(aes(y=NULL, color=Novelty, position="jitter"))
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
    distGraph <- distGraph + geom_density(alpha=0.5)
    distGraph <- distGraph + geom_rug(aes(y=NULL, position="jitter"))
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
# Actually invoke the above plotting functions 
# -------------------------------------------------------

# load the data.
if ( onCMDLine || LOAD_DATA ) {
  bySampleReport <- gsa.read.gatkreport(bySampleEval)
  byACReport <- gsa.read.gatkreport(byACEval)
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

if ( ! is.na(outputPDF) ) {
  dev.off()
  if (exists("compactPDF")) {
    compactPDF(outputPDF)
  }
}
