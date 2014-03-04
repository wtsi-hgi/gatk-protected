require("ggplot2")

chr20Ns = 3520000
processDataset <- function (filename) {
  x = read.table(filename, header=T)
  x[1,]$Filtered = x[1,]$Filtered - chr20Ns
  x = data.frame(x$Coverage, x$Filtered, x$Filtered/mean(x$Filtered))
  names(x) = c("Coverage", "Bases", "Normalized")
  x
}
p101 = read.table("~/sandbox/101.tbl", header=T)
d101 = read.table("~/sandbox/trio.tbl", header=T)

xd101 = processDataset("~/sandbox/CEUTrio.HiSeq.WGS.b37.NA12878.grp")
xp101 = processDataset("~/sandbox/PCRFree.2x101.Illumina.WGS.b37.NA12878.grp")
davg = sum(as.numeric(xd101$Coverage*xd101$Bases))/sum(xd101$Bases)
pavg = sum(as.numeric(xp101$Coverage*xp101$Bases))/sum(xp101$Bases)

d = merge(p101, d101, by="POS")
names(d) = c("position", "pcrfree", "pcrplus")

p = qplot(pcrplus, pcrfree, data=d, geom="point", alpha=0.1, xlim=c(7,50), ylim=c(7,50)) 
p = p + geom_abline(intercept=0, slope=1, linetype="longdash")
p = p + xlab("PCR+") + ylab("PCR-")
p = p + ggtitle("Coverage over poorly covered targets")
p = p + theme_bw()
p = p + geom_rug(alpha=0.3)
p = p + scale_alpha(guide="none")
p