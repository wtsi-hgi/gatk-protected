require("gsalib")
require("ggplot2")
require("gplots")
require("tools")

# TODOs:
#  Assumes you have indels in your call set.  If not you will get errors
#  Create pre/post calling sections
#  Allow conditional use of the preQCFile (where it's not available)

args <- commandArgs(TRUE)

onCMDLine <- ! is.na(args[1])
LOAD_DATA <- F

# creates an array of c(sampleName1, ..., sampleNameN)
parseHighlightSamples <- function(s) {
  return(unlist(strsplit(s, ",", fixed=T)))
}

preQCFile <- NA
if ( onCMDLine ) {
  ProjectName <- args[1]
  VariantEvalRoot <- args[2]
  outputPDF <- args[3]
  if ( ! is.na(args[4]) ) 
    preQCFile <- args[4]
  if ( ! is.na(args[5]) ) 
    highlightSamples <- parseHighlightSamples(args[5])
  else
    highlightSamples <- c()
} else {
  ProjectName <- "InDevelopmentInR"
  VariantEvalRoot <- "t2dChr20Eval/ALL.chr20.freeze20110608_umich.genotypes.vcf.gz"
  outputPDF <- "variantCallQC.pdf"
  preQCFile <- NA # "~/Desktop/broadLocal/GATK/trunk/qcTestData/GoT2D_exomes_batch_005_per_sample_metrics.tsv"
  highlightSamples <- c() # parseHighlightSamples("29029,47243")
}

print("Report")
print(paste("Project          :", ProjectName))
print(paste("VariantEvalRoot  :", VariantEvalRoot))
print(paste("outputPDF        :", outputPDF))
print(paste("preQCFile        :", preQCFile))
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

createMetricsBySites <- function(VariantEvalRoot, PreQCMetrics) {
  # Metrics by sites:
  #  bySite -> counts of SNPs and Indels by novelty, with expectations
  #  byAC -> snps and indels (known / novel)
  bySiteEval <- expandVEReport(selectCumulativeMetrics(gsa.read.gatkreport(paste(VariantEvalRoot, ".bySampleByFunctionalClass.eval", sep="")),'Sample'))
  byACEval <- gsa.read.gatkreport(paste(VariantEvalRoot, ".byAC.eval", sep=""))
  r <- list(bySite = bySiteEval,byAC = byACEval)
  r$byAC$CountVariants$nIndels <- r$byAC$CountVariants$nInsertions + r$byAC$CountVariants$nDeletions
  r$byAC$TiTvVariantEvaluator$nSNPs <- r$byAC$TiTvVariantEvaluator$nTi + r$byAC$TiTvVariantEvaluator$nTv
  r$byAC$CountVariants$AC <- r$byAC$CountVariants$AlleleCount
  r$byAC$TiTvVariantEvaluator$AC <- r$byAC$TiTvVariantEvaluator$AlleleCount
  return(r)
}

summaryTable <- function(metricsBySites, metricsBySample) {
  # SNP summary statistics
  merged <- merge(metricsBySites$bySite$CountVariants, metricsBySites$bySite$TiTvVariantEvaluator)
  sub <- subset(merged, FunctionalClass=="all")
  raw <- melt(sub, id.vars=c("Novelty"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "deletionInsertionRatio"))
  table <- cast(raw, Novelty ~ ...)
  # doesn't work with textplot
  colnames(table) <- c("Novelty", "Target size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "deletion/insertion ratio")
  return(table)
}

sampleSummaryTable <- function(metricsBySamples) {
  # SNP summary statistics
  raw <- melt(metricsBySamples, id.vars=c("Novelty", "Sample"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "deletionInsertionRatio"))
  table <- cast(raw, Novelty ~ variable, mean)
  table$nSNPs <- round(table$nSNPs, 0)
  table$nIndels <- round(table$nIndels, 0)
  table$tiTvRatio <- round(table$tiTvRatio, 2)
  table$deletionInsertionRatio <- round(table$deletionInsertionRatio, 2)
  colnames(table) <- c("Novelty", "Target size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "deletion/insertion ratio")
  return(table)
}

overallSummaryTable <- function(metricsBySites, metricsBySamples) {
  sitesSummary <- as.data.frame(summaryTable(metricsBySites, metricsBySamples))
  sitesSummary$Metric.Type <- "Sites"
  sampleSummary <- as.data.frame(sampleSummaryTable(metricsBySamples))
  sampleSummary$Metric.Type <- "Per-sample avg."
  # that last item puts the metric.type second in the list
  return(rbind(sitesSummary, sampleSummary)[, c(1,7,2,3,4,5,6)])
}

summaryPlots <- function(metricsBySites) {
  # See note about ggplot with less than two points in perSamplePlots.
  numAC <- max(metricsBySites$byAC$TiTvVariantEvaluator$AC)

  name <- "SNP and Indel count by novelty and allele frequency" 
  molten <- melt(subset(metricsBySites$byAC$CountVariants, Novelty != "all" & AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(c("nSNPs", "nIndels")))
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
    molten <- melt(subset(metricsBySites$byAC$CountVariants, AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(measure))
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
  if (sum(metricsBySites$byAC$TiTvVariantEvaluator$nSNPs) > 0) {
    byACNoAll <- subset(metricsBySites$byAC$TiTvVariantEvaluator, Novelty != "all" & AC > 0 & nSNPs > minSNPsToInclude)
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
  name <- "SNPs to indels ratio by allele frequency" 
  if ( sum(metricsBySites$byAC$CountVariants$nIndels) > 0 & sum(metricsBySites$byAC$CountVariants$nSNPs) > 0 ) {
    metricsBySites$byAC$CountVariants$SNP.Indel.Ratio <- metricsBySites$byAC$CountVariants$nSNPs / metricsBySites$byAC$CountVariants$nIndels
    metricsBySites$byAC$CountVariants$SNP.Indel.Ratio[metricsBySites$byAC$CountVariants$nIndels == 0] <- NaN
    p <- ggplot(data=subset(metricsBySites$byAC$CountVariants, Novelty == "all" & nSNPs > 0), aes(x=AC, y=SNP.Indel.Ratio))
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
  
  name <- "SNP counts by functional class" 
  if (sum(metricsBySites$byAC$CountVariants$nSNPs) > 0) {
    molten <- melt(subset(metricsBySites$bySite$CountVariants, Novelty != "all" & FunctionalClass != "all"), id.vars=c("Novelty", "FunctionalClass"), measure.vars=c(c("nSNPs")))
    if ( sum(molten$value) > 0 ) {
      p <- ggplot(data=molten, aes(x=FunctionalClass, y=value, fill=Novelty), group=FunctionalClass)
      p <- p + opts(title = name)
      p <- p + scale_y_log10("No. of SNPs")
      p <- p + geom_bar(position="dodge")
      print(p)
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

createMetricsBySamples <- function(VariantEvalRoot) {
  bySampleEval <- expandVEReport(selectCumulativeMetrics(gsa.read.gatkreport(paste(VariantEvalRoot,".bySampleByFunctionalClass.eval", sep="")),'FunctionalClass'))
  r <- merge(bySampleEval$TiTvVariantEvaluator, bySampleEval$CountVariants)
  r <- merge(r, bySampleEval$CompOverlap)
  if ( ! is.na(preQCFile) ) {
    preQCMetrics <- read.table(preQCFile, header=T)
    r <- merge(r, preQCMetrics)
  }
  # order the samples by nSNPs -- it's the natural ordering.
  x <- subset(r, Novelty=="all")
  r$Sample <- factor(x$Sample, levels=x$Sample[order(x$nSNPs)])

  # add highlight info
  r$highlight <- r$Sample %in% highlightSamples

  #r <- merge(merge(preQCMetrics, byACEval$TiTvVariantEvaluator), byACEval$CountVariants)
  return(subset(r, Sample != "all"))
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
  metricsBySites <- createMetricsBySites(VariantEvalRoot)
  metricsBySamples <- createMetricsBySamples(VariantEvalRoot)
}

if ( ! is.na(outputPDF) ) {
  pdf(outputPDF, height=8.5, width=11)
}

# Table of overall counts and quality
textplot(overallSummaryTable(metricsBySites, metricsBySamples), show.rownames=F)
title(paste("Summary metrics for project", ProjectName), cex=3)

summaryPlots(metricsBySites)
perSamplePlots(metricsBySamples)

if ( ! is.na(outputPDF) ) {
  dev.off()
  if (exists("compactPDF")) {
    compactPDF(outputPDF)
  }
}
