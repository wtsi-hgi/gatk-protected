chr20Ns = 3520000


calcCoverage <- function (x) {
  x[1,]$Filtered = x[1,]$Filtered - chr20Ns
  avg = sum(as.numeric(x$Coverage*x$Filtered)) / sum(x$Filtered)
  p = avg * 0.2
  c(sum(x[1:p,]$Filtered)/ sum(x$Filtered), sum(x[1,]$Filtered)/ sum(x$Filtered))
}

calcCoverage(read.table("~/Dropbox/sandbox/PCRFree.2x101.Illumina.WGS.b37.NA12878.grp", header=T))
calcCoverage(read.table("~/Dropbox/sandbox/PCRFree.2x250.bwasw.grp", header=T))
calcCoverage(read.table("~/Dropbox/sandbox/CEUTrio.HiSeq.WGS.b37.NA12878.grp", header=T))
calcCoverage(read.table("~/Dropbox/sandbox/PCRFree.1x32.grp", header=T))
calcCoverage(read.table("~/Dropbox/sandbox/PCRFree.2x400.grp", header=T))
calcCoverage(read.table("~/Dropbox/sandbox/CEUTrio.HiSeq.WEx.b37.NA12878.grp", header=T))

