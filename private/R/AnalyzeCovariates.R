library("ggplot2")
data = read.table("NA12892.recal.dat", head=T)
data$Accuracy = data$Qempirical - data$Qreported # calculate the delta in accuracy up front as a short hand
numRG = length(unique(data$ReadGroup))
blankTheme = opts(
    panel.grid.major = theme_blank(),
    panel.grid.minor = theme_blank(),
    panel.background = theme_blank(),
    axis.ticks = theme_blank()
  )

pdf("NA12892.recal.pdf",height=7,width=11)
for(cov in levels(data$Covariate)) { # for each covariate in turn  
  d = data[data$Covariate==cov,] # pull out just the data for this covariate so we can treat the non-numeric values appropriately
  if( cov == "Dinuc" ) {
    d$Value = as.character(d$Value)
  } else {
    d$Value = as.numeric(levels(d$Value))[as.integer(d$Value)] # efficient way to convert factors back to their real values
    if( cov == "QualityScore" ) {
      d = subset(d,Value > 5) # cutting off low qual bases here, but eventually this will be done in the BQSRv2 code
    }
  }
  # plot Accuracy by the covariate value
  p <- ggplot(d, aes(x=Value,y=Accuracy,size=nMismatches)) +
    geom_abline(intercept=0, slope=0, linetype=2) + 
    xlab(paste(cov,"Covariate")) +
    ylab("Base Quality Score Accuracy (Empirical - Reported)") +
    blankTheme
  print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) +
    opts(title = paste("Base recalibration accuracy", cov, "All read groups", sep=" - ")))
  print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_wrap(~ ReadGroup) +
    opts(title = paste("Base recalibration accuracy", cov, "Per read group", sep=" - ")) +
    opts(axis.text.x=theme_text(angle=70, hjust=0)))
  
  # plot Reported quality versus Empirical quality as a function of the covariate value
  orig = subset(d, Recalibration == "Original")
  rmseOriginal = sqrt(sum((orig$Qreported-orig$Qempirical)^2 * orig$nBases) / sum(orig$nBases))
  recal = subset(d, Recalibration == "Recalibrated")
  rmseRecalibrated = sqrt(sum((recal$Qreported-recal$Qempirical)^2 * recal$nBases) / sum(recal$nBases))   
  p <- ggplot(d, aes(x=Qreported,y=Qempirical,size=nMismatches)) +
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
  if( cov != "QualityScore" ) {
    p <- ggplot(d, aes(x=Value,y=Qreported,size=nMismatches)) +
      xlab(paste(cov,"Covariate")) +
      ylab("Mean Base Quality Score") +
      blankTheme
    print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) +
      opts(title = paste("Mean base quality score", cov, "All read groups", sep=" - ")))
    print(p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_wrap(~ ReadGroup) +
      opts(title = paste("Mean base quality score", cov, "Per read group", sep=" - ")) +
      opts(axis.text.x=theme_text(angle=70, hjust=0)))
  }

  # plot a histogram of the number of observations of the covariate value
  p <- ggplot(d, aes(x=Value)) +
    xlab(paste(cov,"Covariate")) +
    ylab("Number of Observations") +
    blankTheme
  print(p + geom_histogram(aes(fill=Recalibration,weight=nBases),alpha=0.6,binwidth=1,position="identity") + scale_fill_manual(values=c("maroon1","blue")) +      
    scale_y_continuous(formatter="comma") + opts(title = paste("Histogram of observations", cov, "All read groups", sep=" - ")))
  print(p + geom_histogram(aes(fill=Recalibration,weight=nBases),position="identity",alpha=0.6,binwidth=1) + 
    scale_fill_manual(values=c("maroon1","blue")) + facet_wrap(~ ReadGroup) + scale_y_continuous(formatter="comma") +
    opts(title = paste("Histogram of observations", cov, "Per read group", sep=" - ")) +
    opts(axis.text.x=theme_text(angle=70, hjust=0)))
}
dev.off()
