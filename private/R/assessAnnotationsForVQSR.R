library(ggplot2)
library(entropy)
#d <- read.table("~/Desktop/broadLocal/GATK/unstable/vqsr.table", header=T)
#d <- read.table("~/Desktop/broadLocal/GATK/unstable/vqsr.indels.table", header=T)

training = subset(d, ! is.na(POSITIVE_TRAIN_SITE) | ! is.na(NEGATIVE_TRAIN_SITE))

training$train = "NONE"
training$train[! is.na(training$POSITIVE_TRAIN_SITE)] = "pos"
training$train[! is.na(training$NEGATIVE_TRAIN_SITE)] = "neg"
training$train = factor(training$train)

p <- ggplot(training, aes(x=LikelihoodRankSumTest, color=train))
p <- p + geom_density()
print(p)

forRFF <- na.omit(training[,c("FS", "ReadPosRankSum", "LikelihoodRankSumTest", "BaseQRankSum", "QD", "MQRankSum", "DP", "train")])
#forRFF <- forRFF
forRFF <- forRFF[sample(dim(forRFF)[1], 10000),]
rrf <- RRF(forRFF[,c("FS", "ReadPosRankSum", "BaseQRankSum", "LikelihoodRankSumTest", "QD", "MQRankSum", "DP")], y=forRFF$train)
print(rrf)
varImpPlot(rrf)

mutualInformation <- function(feature, breaks, label) {
  choppedFeature = cut(feature, breaks)
  #print(table(choppedFeature, label))
  mi.plugin(prop.table(table(choppedFeature, label)), unit="log2")
}

mutualInformation(forRFF$DP, seq(0, 600, 5), forRFF$train)
mutualInformation(forRFF$FS, 0:100, forRFF$train)
mutualInformation(forRFF$QD, 0:40, forRFF$train)
mutualInformation(forRFF$MQRankSum, seq(-15, 15, 0.5), forRFF$train)
mutualInformation(forRFF$ReadPosRankSum, seq(-15, 15, 0.5), forRFF$train)
mutualInformation(forRFF$BaseQRankSum, seq(-15, 15, 0.5), forRFF$train)
mutualInformation(forRFF$LikelihoodRankSumTest, seq(-15, 15, 0.5), forRFF$train)
mutualInformation(rnorm(length(forRFF$train)), seq(-15, 15, 0.5), forRFF$train)
