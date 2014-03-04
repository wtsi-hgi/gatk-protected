library(gsalib)
library(ggplot2)
library(scales)

d <- gsa.read.gatkreport("~/Downloads/agbt2013/pcr_bqsr.gatkreport.txt")
assess <- d$NA12878Assessment
assess$Name = factor(assess$Name, levels=c("hc_pcrfree_32_bqsr", "ug_pcr_nobqsr", "ug_pcr_bqsr", "ug_pcrfree_nobqsr", "ug_pcrfree_bqsr", "hc_pcr_nobqsr", "hc_pcr_bqsr", "hc_pcrfree_nobqsr", "hc_pcrfree_bqsr", "hc_pcrfree_250_nobqsr", "hc_pcrfree_250_bqsr", "ug_pcrfree_250_bqsr", "hc_bubble_pcrfree_250_nobqsr", "hc_bubble_pcrfree_250_bqsr"), ordered=T)
assess = assess[order(assess$Name),]

assoc <- function(data) {
  m = matrix(data, nrow = 2)
  fisher.test(m)
}

assoc(c(10, 1000, 25, 2000))
assoc(c(10, 1000, 1, 2000))
assoc(c(250, 1000, 400, 2000))

assoc(c(10, 1000, 0, 500))
assoc(c(10, 1000, 1, 2000))
assoc(c(10, 1000, 10, 20000))
assoc(c(10, 1000, 100, 200000))


assessTable <- function(data) {
  data$FP_BURDEN = data$REASONABLE_FILTERS_WOULD_FILTER_FP_SITE + data$FALSE_POSITIVE_SITE_IS_FP
  data$UNKNOWN = data$CALLED_NOT_IN_DB_AT_ALL + data$CALLED_IN_DB_UNKNOWN_STATUS
  data$FN = data$FALSE_NEGATIVE_NOT_CALLED_AT_ALL + data$FALSE_NEGATIVE_CALLED_BUT_FILTERED
  data$TOTAL_TRUE_VARIANTS = data$FN + data$TRUE_POSITIVE
  print(data[c("Name", "VariantType", "TOTAL_TRUE_VARIANTS", "TRUE_POSITIVE", "FP_BURDEN", "UNKNOWN", "FN")])
}

snps = subset(assess, VariantType == "SNPS")
indels = subset(assess, VariantType == "INDELS")

# for pcr vs. pcr-free
assessTable(subset(rbind(snps, indels), Name %in% c("ug_pcr_nobqsr", "ug_pcrfree_nobqsr")))

# for HC smith-waterman vs. bubble
assessTable(subset(rbind(snps, indels), Name %in% c("hc_bubble_pcrfree_250_bqsr", "hc_pcrfree_250_bqsr")))

# for bqsr vs no-bqsr
assessTable(subset(rbind(snps, indels), Name %in% c("ug_pcrfree_nobqsr", "hc_pcrfree_nobqsr", "hc_pcrfree_bqsr", "hc_pcrfree_250_nobqsr", "hc_pcrfree_250_bqsr")))

# for Mauricio
#assessTable(subset(rbind(snps, indels), Name %in% c("ug_pcrfree_bqsr", "hc_pcrfree_32_bqsr", "ug_pcrfree_250_bqsr")))

#x = subset(assess, VariantType == "indel" & (Name == "hc_pcr_nobqsr" | Name == "hc_pcr_bqsr"))
assessTable(gsa.read.gatkreport("~/Downloads/agbt2013/oddKB.gatkreport.txt")$NA12878Assessment)

if ( ! exists("bqsr.pcrfree") || ! exists("bqsr.pcr")) {
  bqsr.pcr <- gsa.read.gatkreport("~/Downloads/agbt2013/bqsr/CEUTrio.HiSeq.WGS.b37.NA12878.grp")$RecalTable2
  bqsr.pcr$PCRType = "PCR"
  bqsr.pcrfree <- gsa.read.gatkreport("~/Downloads/agbt2013/bqsr/PCRFree.2x101.Illumina.WGS.b37.NA12878.grp")$RecalTable2
  bqsr.pcrfree$PCRType = "PCR-free"
  bqsr = rbind(bqsr.pcr, bqsr.pcrfree)
  
  cycle <- subset(bqsr, CovariateName == "Cycle" & QualityScore == 45) 
  cycle$CovariateValue = as.numeric(as.character(cycle$CovariateValue))
  cycle$ErrorRate = 10^(cycle$EmpiricalQuality/-10)
  
  p = ggplot(data=subset(cycle, EventType == "D" & CovariateValue > 0), aes(x=CovariateValue, y=ErrorRate, color=PCRType))
  p = p + theme_bw()
  #p = p + geom_point(size=3)
  p = p + geom_line(size=1.5)
  p = p + xlab("Machine cycle") + ylab("Error rate by cycle")
  p = p + ggtitle("Deletion error rate by cycle")
  p = p + theme(text = element_text(size=18))
  p = p + scale_y_continuous(labels = scientific)
  p = p + opts(legend.position="none")
  print(p)
}

#p = ggplot(data=subset(bqsr, CovariateName == "Context" & QualityScore==45 & EventType == "D"), aes(x=CovariateValue, y=EmpiricalQuality, color=PCRType))
#p = p + geom_point()
#print(p)

# cov <- gsa.read.gatkreport("~/Downloads/datasets/PCRFree.2x101.Illumina.WGS.b37.NA12878.grp")$BaseCoverageDistribution
# cov[1,]$Count = cov[1,]$Count - 3.5e6 
# cov[1,]$Filtered = cov[1,]$Filtered - 3.5e6 
# 
# 
# qqnorm(cov$Filtered)
# qqnorm(rbinom(100000, size=2197166, prob=100/63e6))
# 
# x <- rnorm(1000)
# qqplot(x, qnorm(x))
# 
# size = 2197166
# prob = 100/59e6
# x <- rbinom(100000, size=size, prob=prob)
# qqplot(x, rbinom(length(x), size=size, prob=prob))
# 
# selCov <- subset(cov, Coverage < 80)
# x = sample(selCov$Coverage, 10000, replace=T, prob=selCov$Count/sum(selCov$Count))
# qqplot(x, rbinom(length(x), size=size, prob=prob))


# p = 100 / 63e6
# readCov <- function(name, file) {
#   x = gsa.read.gatkreport(file)$BaseCoverageDistribution
#   x$type = name
#   x$FilteredNoNs <- x$Filtered
#   x$FilteredNoNs[1] = x$FilteredNoNs[1] - chr20NCount
#   x$RelativeCoverage = x$FilteredNoNs / sum(x$FilteredNoNs)
#   #x$MeanCoverage = weighted.mean(x$Coverage, x$RelativeCoverage)
#   x$MeanCoverage = x$Coverage[x$RelativeCoverage == max(x$RelativeCoverage)]
#   x$MeanAt0 = x$Coverage - x$MeanCoverage
#   x$RelativeVariantAdjustedCoverage = x$RelativeCoverage / x$MeanCoverage
#   x
# }
# 
# chr20NCount <- 63025520 - 59505520
# cov = rbind(
#   readCov("PCR-free-2x250", "~/Desktop/broadLocal/GATK/unstable/pcrfreetables/PCRFree.2x250.bwasw.grp"),
#   readCov("PCR-2x101", "~/Desktop/broadLocal/GATK/unstable/pcrfreetables/CEUTrio.HiSeq.WGS.b37.NA12878.grp"),
#   readCov("PCR-free-2x101", "~/Desktop/broadLocal/GATK/unstable/pcrfreetables/PCRFree.2x101.Illumina.WGS.b37.NA12878.grp"),
#   readCov("PCR-free-2x400", "~/Desktop/broadLocal/GATK/unstable/pcrfreetables/PCRFree.2x400.grp")
#   #readCov("PCR-free-1x32", "~/Desktop/broadLocal/GATK/unstable/pcrfreetables/PCRFree.1x32.grp")
# )
# print(head(subset(cov, Coverage == 0 | Coverage == 1)))
# 
# p <- ggplot(cov, aes(x=Coverage,y=RelativeCoverage, color=type))
# p <- p + geom_line()
# p <- p + geom_vline(aes(xintercept=MeanCoverage, color=type), linetype="dashed")
# #p <- p + scale_y_log10()
# p <- p + xlim(0,150)
# print(p)
# 
# p <- ggplot(cov, aes(x=MeanAt0,y=RelativeCoverage, color=type))
# p <- p + geom_line()
# #p <- p + scale_y_log10()
# p <- p + xlim(-100,100)
# print(p)

