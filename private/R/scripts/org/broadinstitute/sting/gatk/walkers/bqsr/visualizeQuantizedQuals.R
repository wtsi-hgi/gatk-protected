library(ggplot2)
library(gsalib)
library(gplots)

data <- gsa.read.gatkreport("~/Desktop/broadLocal/GATK/unstable/quantize.quals.report.txt")
qualHist <- data$QualHistogram
intervals <- data$QualQuantizerIntervals

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

MAX_QUAL = 50

p = ggplot(data=intervals, aes(x=merge.order, y=qual, ymin=qStart, ymax=qEnd+0.5, color=root.node, label=penalty))
p = p + geom_crossbar()
p = p + geom_text(aes(y=50))
p = p + ylim(0, MAX_QUAL+1)
p = p + coord_flip() 
p = p + opts(legend.position="top")
intervalsGraph <- p
#print(p)

p = ggplot(data=qualHist, aes(x=qual, y=count))
p = p + geom_bar(stat="identity")
p = p + xlim(0, MAX_QUAL+1)
histGraph <- p
print(histGraph)

distributeGraphRows(list(intervalsGraph, histGraph), c(2,1))

printIntervals <- function(d) {
  d = d[order(d$qStart),]
  print(d)
}
printIntervals(subset(intervals, root.node == "true"))