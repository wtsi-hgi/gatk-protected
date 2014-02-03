require("ggplot2")
require("reshape2")

add_dataset <- function(rlen, filename, m, meancov) {
  d <- cbind(rlen, read.delim(filename))
  names(d) <- c("dataset", "position", "raw_depth", "not_depth")
  d <- cbind(d, (d$not_depth/meancov))
  names(d) <- c("dataset", "position", "raw_depth", "not_depth", "depth")
  merge(m, d, by="position")
}

old = read.delim("~/Dropbox/sandbox/CEUTrio.HiSeq.WGS.b37.NA12878ug.bqsr.vcf.large_indels.tbl")
names(old) <- c("position", "raw_depth", "depth")
master = data.frame(old$position, old$depth/64)
names(master) <- c("position", "old")

x = add_dataset(101, "~/Dropbox/sandbox/PCRFree.2x101.Illumina.WGS.b37.NA12878ug.bqsr.vcf.large_indels.tbl", master, 27)
x = add_dataset(250, "~/Dropbox/sandbox/NA12878-2x250.bwasw.chr20.dedup.clean.ug.bqsr.vcf.large_indels.tbl", x, 41)
x = add_dataset(32, "~/Dropbox/sandbox/NA12878-32.se.dedup.clean.recal.hc.bqsr.vcf.large_indels.tbl", x, 3)

#d = data.frame(x$old, x$depth.x, x$depth.y, x$raw_depth/3.8)
#names(d) <- c("old", "101", "250", "32")
d = data.frame(x$old, x$depth.x, x$depth.y)
names(d) <- c("old", "101", "250")
m = melt(d, id.vars=c("old"))

qplot(old, value, data=m, group=variable, color=variable,geom="jitter", ylab="new coverage", xlab = "old coverage", log="xy" ) + geom_abline(intercept=0, slope=1)