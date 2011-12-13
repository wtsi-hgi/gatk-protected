library(gsalib)
library(ggplot2)
#library(gplots)
library(tools)

args <- commandArgs(TRUE)

onCMDLine <- ! is.na(args[1])
LOAD_DATA <- ! exists("allReports")

if ( onCMDLine ) {
  file <- args[1]
  outputPDF <- args[2]
} else {
  file <- "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/gatkPerformanceOverTime/Q-5126@gsa1.jobreport.txt"
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
  p = ggplot(data=report, aes(x=nSamples, y=runtime, group=gatk, color=gatk))
  p = p + facet_grid(. ~ assessment, scales="free")
  #p = p + geom_jitter()
  #p = p + geom_point()
  p = p + geom_smooth()
  p = p + scale_x_log10() + scale_y_log10()
  p = p + opts(title=report$analysisName)
  p = p + geom_boxplot(aes(group=interaction(nSamples, gatk)), outlier.colour="blue")
  print(p)
}

plotNormalizedByNSamples <- function(report) {
  norm = ddply(report, .(nSamples, assessment), transform, normRuntime = runtime / mean(runtime))
  p = ggplot(data=norm, aes(x=nSamples, y=normRuntime, group=gatk, color=gatk))
  p = p + facet_grid(. ~ assessment, scales="free")
  p = p + geom_jitter()
  #p = p + geom_point()
  p = p + geom_smooth()
  p = p + scale_x_log10()# + scale_y_log10()
  p = p + opts(title=paste("Runtime per nSamples relative to nSamples mean value", report$analysisName))
  #p = p + geom_boxplot(aes(group=interaction(nSamples, gatk)), outlier.colour="blue")
  print(p)
}


convertUnits <- function(gatkReportData) {
  convertGroup <- function(g) {
    g$runtime = g$runtime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g$startTime = g$startTime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g$doneTime = g$doneTime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g
  }
  lapply(gatkReportData, convertGroup)
}

# -------------------------------------------------------
# Actually invoke the above plotting functions 
# -------------------------------------------------------

RUNTIME_UNITS = "(hours)"
ORIGINAL_UNITS_TO_RUNTIME_UNITS = 1/1000/60/60

# load the data.
if ( onCMDLine || LOAD_DATA ) {
  allReports <- gsa.read.gatkreport(file)
  allReports <- convertUnits(allReports)
}

if ( ! is.na(outputPDF) ) {
  pdf(outputPDF, height=8.5, width=11)
}

# actually do something
for ( report in list(allReports$CountLoci, allReports$UnifiedGenotyper) ) {
  print(head(report))
  plotByNSamples(report)
  plotNormalizedByNSamples(report)
}


if ( ! is.na(outputPDF) ) {
  dev.off()
  if (exists("compactPDF")) {
    compactPDF(outputPDF)
  }
}
