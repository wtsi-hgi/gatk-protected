library(gsalib)
require("ggplot2")
args = commandArgs(TRUE)
onCMDLine = ! is.na(args[1])

if ( onCMDLine ) {
  inputFileName = args[1]
} else {
  inputFileName = "~/Desktop/broadLocal/GATK/unstable/report.txt"
}

allJobsFromReport <- function(report) {
  names <- c("jobName", "startTime", "analysisName", "doneTime")
  sub <- lapply(report, function(table) table[,names])
  r = do.call("rbind", sub)
  # todo -- sort by time here
  r
}

plotJobsGantt <- function(gatkReport, sortOverall) {
  allJobs = allJobsFromReport(gatkReport)
  if ( sortOverall ) {
    title = "All jobs, by analysis, by start time"
    allJobs = allJobs[order(allJobs$analysisName, allJobs$startTime, decreasing=T), ]
  } else {
    title = "All jobs, sorted by start time"
    allJobs = allJobs[order(allJobs$startTime, decreasing=T), ]
  }
  allJobs$index = 1:nrow(allJobs)
  minTime = min(allJobs$startTime)
  allJobs$relStartTime = allJobs$startTime - minTime
  allJobs$relDoneTime = allJobs$doneTime - minTime
  p <- ggplot(data=allJobs, aes(x=relStartTime, y=index, color=analysisName))
  p <- p + geom_segment(aes(xend=relDoneTime, yend=index), size=3, arrow=arrow(length = unit(0.3, "cm")))
  p <- p + xlab("Start time (relative to first job) (ms)")
  p <- p + opts(title=title)
  #p <- p + scale_x_datetime(format = "%d %H:%M:%S")
  print(p)
}

standardColumns = c("jobName", "startTime", "formattedStartTime", "analysisName", "intermediate", "formattedDoneTime", "doneTime", "runtime")
plotGroup <- function(groupTable) {
  name = unique(groupTable$analysisName)[1]
  groupAnnotations = setdiff(names(groupTable),standardColumns)
  
  groupTable$annotatedName = "x"
  for ( i in 1:nrow(groupTable)) {
    parts = lapply(groupAnnotations, function(x) paste(x, groupTable[i, x], sep="="))
    groupTable[i,]$annotatedName = do.call("paste", c(groupTable[i,c("jobName")], parts, sep="/"))
  }
  
  print(groupTable)
  p <- ggplot(data=groupTable, aes(x=annotatedName, y=runtime))
  p <- p + geom_bar()
  p <- p + xlab("Jobs")
  p <- p + opts(title=paste(name, ": runtime job"))
  print(p)  
}

print("Report")
print(paste("Project          :", inputFileName))

# parseTimes <- function(report) {
#   timeFormat = "%d.%m.%y/%H:%M:%S"
#   one <- function(g) {
#     g$parsedDoneTime <- strptime(g$doneTime, timeFormat)
#     g$parsedStartTime <- strptime(g$startTime, timeFormat)
#     print(g)
#   }
#   lapply(report, one)
# }

gatkReportData <- gsa.read.gatkreport(inputFileName)
#gatkReportData <- parseTimes(gatkReportData)
#print(gatkReportData$BigCombine)
print(summary(gatkReportData))

# plotJobsGantt(gatkReportData, T)
# plotJobsGantt(gatkReportData, F)
# for ( group in gatkReportData ) {
#   plotGroup(group)
# }
plotGroup(gatkReportData$SitesVsGenotypes)