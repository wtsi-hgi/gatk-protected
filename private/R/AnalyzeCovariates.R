## Assumes the data is organized like this:
#Recalibration  ReadGroup	Covariate	Value	QualityScore	EmpiricalQuality	Errors	Observations
#Original	    B00EG.4	    Cycle	    2	    2.000	        11.825	            9011003	137161904
#Original	    B00EG.4	    Cycle	    4	    4.000	        5.001	            526054	1663850
#Original	    B00EG.4	    Cycle	    5	    5.000	        5.817	            522185	1993294
#Original	    B00EG.4	    Context	    TGCCGT	6.000	        6.725	            937794	4411346
#Original	    B00EG.4	    Context	    ACGTTC	7.000	        7.476	            1028074	5749046
#Original	    B00EG.4	    Cycle	    8	    8.000	        8.006	            883264	5580387

library(ggplot2);
library(gsalib);

args <- commandArgs(TRUE)
before = gsa.read.gatkreport(args[1])
after = gsa.read.gatkreport(args[2])
a = cbind("Original", before$RecalTable2)
b = cbind("Recalibrated", after$RecalTable2)
names(a) = c("Recalibration", names(before$RecalTable2))
names(b) = c("Recalibration", names(after$RecalTable2))
data = rbind(a,b)

data$Accuracy = data$EmpiricalQuality - data$QualityScore # calculate the delta in accuracy up front as a short hand
numRG = length(unique(data$ReadGroup))
blankTheme = opts(
    panel.grid.major = theme_blank(),
    panel.grid.minor = theme_blank(),
    panel.background = theme_blank(),
    axis.ticks = theme_blank()
  )

pdf(args[3],height=7,width=11)
for(cov in levels(data$CovariateName)) { # for each covariate in turn
  d = data[data$CovariateName==cov,] # pull out just the data for this covariate so we can treat the non-numeric values appropriately
  if (cov == "Context") 
    d$Value = as.character(d$CovariateValue)
  if (cov == "Cycle")
    d$Value = as.numeric(levels(d$CovariateValue))[as.integer(d$CovariateValue)] # efficient way to convert factors back to their real values

  # plot Accuracy by the covariate value
  p <- qplot(data=d, x=Value, y=Accuracy, size=Errors, alpha=Observations, color=Recalibration, xlab = paste(cov,"Covariate"), ylab = "Base Quality Score Accuracy (Empirical - Reported)") + geom_abline(intercept=0, slope=0, linetype=2) + scale_color_manual(values=c("maroon1","blue")) + blankTheme
  print(p + opts(title = paste("Base recalibration accuracy", cov, "All read groups", sep=" - ")))
  print(p + facet_wrap(~ ReadGroup) + opts(title = paste("Base recalibration accuracy", cov, "Per read group", sep=" - "), axis.text.x=theme_text(angle=70, hjust=0)))

  # plot Reported quality versus Empirical quality as a function of the covariate value
  orig = subset(d, Recalibration == "Original")
  rmseOriginal = sqrt(sum((orig$QualityScore-orig$EmpiricalQuality)^2 * orig$Observations) / sum(orig$Observations))
  recal = subset(d, Recalibration == "Recalibrated")
  rmseRecalibrated = sqrt(sum((recal$QualityScore-recal$EmpiricalQuality)^2 * recal$Observations) / sum(recal$Observations))   
  p <- ggplot(d, aes(x=QualityScore,y=EmpiricalQuality,size=Errors,alpha=Observations)) +
    geom_abline(intercept=0, slope=1, linetype=2) + 
    xlab("Reported Base Quality") +
    ylab("Empirical Base Quality") +
    blankTheme
  print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) +
    opts(title = paste("Base recalibration curve", cov, "All read groups", sep=" - ")))
    grid.text(paste("Original RMSE",round(rmseOriginal,digits=3),sep=" = "), x=0.1, y=0.89, just="left", gp = gpar(col = "maroon1",fontsize=12))
    grid.text(paste("Recalibrated RMSE",round(rmseRecalibrated,digits=3),sep=" = "), x=0.1, y=0.865, just="left", gp = gpar(col = "blue",fontsize=12))
  print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_wrap(~ ReadGroup) +
    opts(title = paste("Base recalibration curve", cov, "Per read group", sep=" - ")))
  
  # plot mean Reported quality as a function of the covariate value
  p <- ggplot(d, aes(x=Value,y=QualityScore,size=Errors,alpha=Observations)) +
    xlab(paste(cov,"Covariate")) +
    ylab("Mean Base Quality Score") +
    blankTheme
  print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) +
    opts(title = paste("Mean base quality score", cov, "All read groups", sep=" - ")))
  print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_wrap(~ ReadGroup) +
    opts(title = paste("Mean base quality score", cov, "Per read group", sep=" - ")) +
    opts(axis.text.x=theme_text(angle=70, hjust=0)))


  # plot a histogram of the number of observations of the covariate value
  p <- ggplot(d, aes(x=Value)) +
    xlab(paste(cov,"Covariate")) +
    ylab("Number of Observations") +
    blankTheme
  print(p + geom_histogram(aes(fill=Recalibration,weight=Observations),alpha=0.6,binwidth=1,position="identity") + scale_fill_manual(values=c("maroon1","blue")) +      
    scale_y_continuous(formatter="comma") + opts(title = paste("Histogram of observations", cov, "All read groups", sep=" - ")))
  print(p + geom_histogram(aes(fill=Recalibration,weight=Observations),position="identity",alpha=0.6,binwidth=1) + 
    scale_fill_manual(values=c("maroon1","blue")) + facet_wrap(~ ReadGroup) + scale_y_continuous(formatter="comma") +
    opts(title = paste("Histogram of observations", cov, "Per read group", sep=" - ")) +
    opts(axis.text.x=theme_text(angle=70, hjust=0)))
}
