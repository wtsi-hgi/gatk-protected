library(ggplot2)
library(reshape2)
library(plyr)

d <- read.table("~/Desktop/broadLocal/GATK/unstable/banded.txt", header=T)
print(head(d))

d$cell.improvement <- d$n.cells.evaluated / d$n.cells.overall
d$runtime.improvement <- d$banded.nano / d$logless.nano 
#d$runtime.improvement <- log10(d$banded.nano / d$logless.nano)
d$banded.nano = log10(d$banded.nano)
d$logless.nano = log10(d$logless.nano)
#d <- subset(d, runtime.improvement < 2 & n.cells.overall < 100000)
molten <- melt(d, id.vars=c("index", "ref.length", "read.length", "max.length"))

p <- ggplot(molten, aes(x=max.length, y=value, color=variable))
#p <- p + geom_point() 
p <- p + geom_smooth() 
#p <- p + geom_boxplot(aes(group=max.length))
p <- p + facet_grid(variable ~ ., scales="free")
#print(p)

d$hapSize <- "<50"
d[d$ref.length > 50 & d$ref.length <= 100,]$hapSize <- "51-100"
d[d$ref.length > 100 & d$ref.length <= 150,]$hapSize <- "101-150"
d[d$ref.length > 150,]$hapSize <- ">150"
p <- ggplot(d, aes(x=logless.n.evaluated, y=n.cells.evaluated, color=hapSize))
p <- p + geom_point()
p <- p + facet_grid(single.run ~ .)
p <- p + geom_abline(slope=1, linetype="dashed")
print(p)

p <- ggplot(d, aes(x=logless.nano, y=banded.nano, color=hapSize))
p <- p + geom_point()
p <- p + facet_grid(single.run ~ .)
p <- p + geom_abline(slope=1, linetype="dashed")
p <- p + geom_abline(intercept=-1, slope=1, linetype="dashed", color="red")
print(p)

p <- ggplot(d, aes(x=logless.nano, y=banded.nano, color=factor(round_any(n.haplotypes,50))))
p <- p + geom_point()
p <- p + facet_grid(single.run ~ .)
p <- p + geom_abline(slope=1, linetype="dashed")
p <- p + geom_abline(intercept=-1, slope=1, linetype="dashed", color="red")
print(p)

d$banded.nano.per.cell = log10(d$banded.nano / d$n.cells.evaluated)
d$logless.nano.per.cell = log10(d$logless.nano / d$logless.n.evaluated)
molten <- melt(d, id.vars=c("index", "hapSize", "single.run"), measure.vars=c("banded.nano.per.cell", "logless.nano.per.cell"))
p <- ggplot(molten, aes(x=value, color=hapSize, linetype=single.run))
p <- p + geom_density()
p <- p + facet_grid(variable ~ .)
print(p)

p <- ggplot(d, aes(x=logless.nano.per.cell, y=banded.nano.per.cell, color=interaction(n.haplotypes)))
p <- p + geom_point()
p <- p + facet_grid(single.run ~ .)
p <- p + geom_abline(slope=1, linetype="dashed")
print(p)

#pre.opt <- read.table("~/Desktop/broadLocal/GATK/unstable/banded.preopt.txt", header=T)
# only works if we are running with the baseline comparator
opt = data.frame(pre.opt=d$logless.nano, cur=d$banded.nano, improvement=d$logless.nano-d$banded.nano, hapSize=d$hapSize)
p <- ggplot(opt, aes(x=pre.opt, y=improvement, color=hapSize))
p <- p + geom_point()
p <- p + geom_abline(slope=0, color="red")
#p <- p + geom_abline(slope=1, color="red")
p <- p + facet_grid(. ~ hapSize)
p <- p + geom_smooth(method="lm")
print(p)

p <- ggplot(opt, aes(x=improvement, fill=hapSize))
p <- p + geom_histogram(binwidth=0.01)
p <- p + facet_grid(hapSize ~ ., scales="free")
p <- p + geom_vline(x=mean(opt$improvement), color="red")
p <- p + xlim(-1,1)
print(p)
