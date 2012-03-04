library(ggplot2)
library(gsalib)

data <- gsa.read.gatkreport("~/Desktop/broadLocal/GATK/unstable/quantize.quals.report.txt")
data <- data$QualQuantizerIntervals

p = ggplot(data=data, aes(x=merge.order, y=qual, ymin=qStart, ymax=qEnd+0.5, color=root.node, label=penalty))
p = p + geom_linerange(size=2)
p = p + geom_text(aes(y=50))
p = p + ylim(0, 51)
p = p + coord_flip() 
print(p)

printIntervals <- function(d) {
  d = d[order(d$qStart),]
  print(d)
}


printIntervals(subset(data, root.node == "true"))