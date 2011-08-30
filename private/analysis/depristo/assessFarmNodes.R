library(gsalib)
require("ggplot2")
require("gplots")

inputFileName = "~/Desktop/broadLocal/GATK/unstable/report.txt"

pdf("farm.pdf")
p <- ggplot(data=gatkReportData$CountingDBSNPRecords, aes(x=runtime))
p <- p + geom_histogram()
print(p)

gatkReportData$CountingDBSNPRecords$relStartTime = gatkReportData$CountingDBSNPRecords$startTime - min(gatkReportData$CountingDBSNPRecords$startTime)
p <- ggplot(data=gatkReportData$CountingDBSNPRecords, aes(x=startTime, y=runtime))
p <- p + geom_point(aes(color=exechosts))
p <- p + geom_smooth(se=T, alpha=0.5)
print(p)

plotBox <- function(minStartTime, maxStartTime) {
  p <- ggplot(data=subset(gatkReportData$CountingDBSNPRecords, startTime > minStartTime & startTime < maxStartTime), aes(x=exechosts, y=runtime, color=exechosts)) 
  p <- p + geom_boxplot(outlier.size=0)
  p <- p + geom_jitter(position=position_jitter(width=0.1))
  p <- p + opts(axis.text.x=theme_text(angle=-45))
  print(p)
}

timeOfSecondRuns = 1.314559e12
plotBox(-Inf, Inf)
plotBox(timeOfSecondRuns, Inf)
plotBox(-Inf, timeOfSecondRuns)
dev.off()

