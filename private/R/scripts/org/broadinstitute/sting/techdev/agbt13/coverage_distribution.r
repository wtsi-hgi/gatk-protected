require("ggplot2")

chr20size = 63025520

process_dataset <- function(dataType, readLength, filename, nReads) {
  d  <- cbind(paste(dataType, readLength, sep=""), read.table(filename, header=T, quote=""))
  d  <- d[with(d, order(Coverage)),]
  d  <- cbind(d, cumsum(d$Filtered)/sum(d$Filtered), d$Coverage/max(d$Coverage), dbinom(x=d$Coverage, size=nReads, prob=readLength/chr20size), d$Filtered/max(d$Filtered))
  names(d) <- c("Sample", "Coverage", "Count", "Filtered", "CumSum", "NormalizedCumulativeCoverage", "BinomialProbability", "NormalizedCoverage")
  return(d)
}

#pcr400  <- process_dataset("pcrfree", 400, "/Users/carneiro/sandbox/PCRFree.2x400.grp", 18605161)
pcr250  <- process_dataset("pcrfree", 250, "/Users/carneiro/sandbox/PCRFree.2x250.bwasw.grp", 14578670)
pcr101  <- process_dataset("pcrfree", 101, "/Users/carneiro/sandbox/PCRFree.2x101.Illumina.WGS.b37.NA12878.grp", 21917166)
pcr32  <- process_dataset("pcrfree", 32,   "/Users/carneiro/sandbox/PCRFree.1x32.grp", 7462217)
old_wgs <- process_dataset("pcrwgs", 101, "/Users/carneiro/sandbox/CEUTrio.HiSeq.WGS.b37.NA12878.grp", 45938316)
#old_wex <- process_dataset("pcrwex", 76, "/Users/carneiro/sandbox/CEUTrio.HiSeq.WEx.b37.NA12878.grp", 3260490)

#d = rbind(pcr400, pcr250, pcr101, pcr32, old_wgs, old_wex)
d = rbind(pcr250, pcr101, pcr32, old_wgs)

x = rep(pcr250$Coverage, times=pcr250$Count)
bin = rbinom(n=10000, size=chr20size, prob=250/14578670)
y = data.frame(x, bin)

qplot(qqplot(bin, x)


qplot(NormalizedCumulativeCoverage, CumSum, data=d, geom="line", group=Sample, color=Sample, 
      xlim = c(0,0.15),
      main = "cumulative sum of loci coverage in PCR-Free samples (BWA short reads)",
      xlab = "normalized coverage ratio of a locus",
      ylab = "normalized cumulative sum of loci covered")

qplot(BinomialProbability, NormalizedCoverage, data=d, geom="point", group=Sample, color=Sample) + geom_abline(intercept=0, slope=1)
qplot(Coverage, NormalizedCoverage/BinomialProbability, data=d, group=Sample, color=Sample, geom="point")

ggplot(data=subset(d, Coverage < 150), aes(x=Coverage, y=log10(BinomialProbability/NormalizedCoverage), group=Sample, color=Sample)) + geom_point()