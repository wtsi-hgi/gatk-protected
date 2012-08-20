library(gsalib)
library(ggplot2)
library(reshape)
#library(gplots)
library(tools)

args <- commandArgs(TRUE)

onCMDLine <- ! is.na(args[1])
LOAD_DATA <- ! exists("allReports")

RUNTIME_UNITS = "(hours)"
ORIGINAL_UNITS_TO_RUNTIME_UNITS = 1/1000/60/60
MIN_RUN_TIME_FOR_SUCCESSFUL_JOB = 10^-3

if ( onCMDLine ) {
  file <- args[1]
  outputPDF <- args[2]
} else {
  #file <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/parallelBQSR/GATKPerformanceOverTime.jobreport.txt"
  #file <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/parallelCombineVariants2/GATKPerformanceOverTime.jobreport.txt"
  file <- "~/Desktop/broadLocal/GATK/unstable/GATKPerformanceOverTime.jobreport.txt"
  #file <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/gatkPerformanceOverTime/Q-24937@gsa1.jobreport.txt"
  outputPDF <- NA
}

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

plotByNSamples <- function(report) {
  plotCmdByX(report, "nSamples", T, logUnit = 10)
}

plotByNT <- function(report, includeFacet = F) {
  plotCmdByX(report, "nt", includeFacet, logUnit = 2)
}

plotCmdByX <- function(report, X, includeFacet = F, logUnit = 10) {
  report$X <- report[[X]]
  p = ggplot(data=report, aes(x=X, y=runtime, group=gatk, color=gatk))
  if ( includeFacet ) 
    p = p + facet_grid(. ~ assessment, scales="free")
  p = p + geom_jitter()
  #p = p + geom_point()
  p = p + geom_smooth()
  if ( logUnit == 10 ) {
    p = p + scale_x_log10(breaks=unique(report$X)) + scale_y_log10()
  } else if ( logUnit == 2 ) {
    p = p + scale_x_continuous(trans = "log2", breaks=unique(report$X)) + scale_y_continuous(trans = "log2")
  }
  p = p + xlab(X) + ylab(paste("Runtime", RUNTIME_UNITS))
  p = p + opts(title=report$analysisName)
  p = p + geom_boxplot(aes(group=interaction(X, gatk)))
  p
}

plotNormalizedByNSamples <- function(report) {
  norm = ddply(report, .(nSamples, assessment), transform, normRuntime = runtime / mean(runtime))
  p = ggplot(data=norm, aes(x=nSamples, y=normRuntime, group=gatk, color=gatk))
  p = p + facet_grid(. ~ assessment, scales="free")
  p = p + geom_boxplot(aes(group=interaction(nSamples, gatk)), outlier.colour="blue")
  #p = p + geom_jitter()
  #p = p + geom_point()
  p = p + geom_smooth()
  p = p + scale_x_log10()# + scale_y_log10()
  p = p + opts(title=paste("Runtime per nSamples relative to nSamples mean value", report$analysisName))
  print(p)
}

plotByGATKVersion <- function(report) {
  #report$runtime <- replicate(length(report$runtime), rnorm(1))
  p = ggplot(data=report, aes(x=gatk, y=runtime, group=gatk, color=gatk))
  p = p + geom_boxplot()
  p = p + geom_jitter()
  #p = p + scale_x_log10()# + scale_y_log10()
  p = p + xlab("GATK version") + ylab(paste("Runtime", RUNTIME_UNITS))
  p = p + opts(title=paste("Runtime", report$analysisName))
  print(p)
}

convertUnits <- function(gatkReportData) {
  convertGroup <- function(g) {
    g$runtime = g$runtime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g$startTime = g$startTime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g$doneTime = g$doneTime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    subset(g, runtime > MIN_RUN_TIME_FOR_SUCCESSFUL_JOB)
  }
  lapply(gatkReportData, convertGroup)
}

# -------------------------------------------------------
# Actually invoke the above plotting functions 
# -------------------------------------------------------

# load the data.
#if ( onCMDLine || LOAD_DATA ) {
  allReports <- gsa.read.gatkreport(file)
  allReports <- convertUnits(allReports)
#}

if ( ! is.na(outputPDF) ) {
  pdf(outputPDF, height=8.5, width=11)
}

# Create reports for per N sample CountLoci and UG
if ( "UnifiedGenotyper" %in% names(allReports) ) {
  for ( report in list(allReports$CountLoci, allReports$UnifiedGenotyper) ) {
    #print(head(report))
    print(plotByNSamples(report))
    plotNormalizedByNSamples(report)
  }
}

# Create reports just doing runtime vs. GATK version
if ( "CombineVariants" %in% names(allReports) ) {
  for ( report in list(allReports$SelectVariants, allReports$CombineVariants,
                       allReports$VariantEval))  {
    #print(head(report))
    plotByGATKVersion(report)
  }
}

getAssessments <- function(report) {
  if ( "assessment" %in% names(report) ) {
    return(levels(report$assessment))
  } else {
    stop(paste("Missing assessment field for report", report[1, "analysisName"]))
  }
}

#
# Plot runtime vs. NT for all of the NT tests
#
ntReports <- c("CombineVariants.nt", "UnifiedGenotyper.nt", "CountLoci.nt", "BaseRecalibrator.nt", "VariantEval.nt")
for ( ntReport in ntReports ) {
  if ( ntReport %in% names(allReports) ) {
    report = allReports[[ntReport]]
    for ( assess in getAssessments(report) ) {
      p = plotByNT(report[report$assessment == assess,])
      p = p + opts(title=paste(ntReport, " performance as a function of nt for ", assess))
      print(p)
    }
  }
}

if ( ! is.na(outputPDF) ) {
  dev.off()
  if (exists("compactPDF")) {
    compactPDF(outputPDF)
  }
}
