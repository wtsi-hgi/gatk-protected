library("ggplot2")

data = read.table("NA12878.20GAV.1.recal.errorModels.noN.dat", head=T,sep="\t")
data <- within(data, ErrorModel <- factor(ErrorModel, levels = rev(levels(ErrorModel))))
data <- within(data, Recalibration <- factor(Recalibration, levels = rev(levels(Recalibration))))

data$Qempirical = -10 * log10((data$nMismatches + 1) / (data$nBases + 1))
data$Accuracy = data$Qempirical - data$Qreported # calculate the delta in accuracy up front as a short hand
numRG = length(unique(data$ReadGroup))
blankTheme = opts(
    panel.grid.major = theme_blank(),
    panel.grid.minor = theme_blank(),
    panel.background = theme_blank(),
    axis.ticks = theme_blank()
  )

# Viewport (layout 2 graphs top to bottom)
distributeGraphRows <- function(graphs, heights = c()) {
  if (length(heights) == 0) {
    heights <- rep.int(1, length(graphs))
  }
  heights <- heights[!is.na(graphs)]
  graphs <- graphs[!is.na(graphs)]
  numGraphs <- length(graphs)
  Layout <- grid.layout(nrow = numGraphs, ncol = 1, heights=heights)
  grid.newpage()
  pushViewport(viewport(layout = Layout))
  subplot <- function(x) viewport(layout.pos.row = x, layout.pos.col = 1)
  for (i in 1:numGraphs) {
    print(graphs[[i]], vp = subplot(i))
  }
}


for(cov in levels(data$Covariate)) { # for each covariate in turn  
  d = data[data$Covariate==cov,] # pull out just the data for this covariate so we can treat the non-numeric values appropriately
  if( cov == "Context" ) {
    d$Value = as.character(d$Value)
  } else {
    d$Value = as.numeric(levels(d$Value))[as.integer(d$Value)] # efficient way to convert factors back to their real values
    if( cov == "QualityScore" ) {
      d = subset(d,Value > 5) # cutting off low qual bases here, but eventually this will be done in the BQSRv2 code
    }
  }
  d=subset(d,nBases>2000) # only show bins which have enough data to acually estimate the quality
  
  if( cov != "QualityScore" ) {    
    p <- ggplot(d, aes(x=Value,y=Accuracy,alpha=log10(nBases))) +
      geom_abline(intercept=0, slope=0, linetype=2) + 
      xlab(paste(cov,"Covariate")) +
      ylab("Quality Score Accuracy") +
      blankTheme
    if(cov == "Cycle") {
      b <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~ErrorModel) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))
      
      p <- ggplot(d, aes(x=Value,y=Qreported,alpha=log10(nBases))) +
      xlab(paste(cov,"Covariate")) +
      ylab("Mean Quality Score") +
      blankTheme
      e <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~ErrorModel) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))
      
      
    } else {
      c <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~ErrorModel) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))
      p <- ggplot(d, aes(x=Value,y=Qreported,alpha=log10(nBases))) +
      xlab(paste(cov,"Covariate")) +
      ylab("Mean Quality Score") +
      blankTheme
      f <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~ErrorModel) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))

    }
  } else {
    p <- ggplot(d, aes(x=Qreported,y=Qempirical,alpha=log10(nBases))) +
      geom_abline(intercept=0, slope=1, linetype=2) + 
      xlab("Reported Quality Score") +
      ylab("Empirical Quality Score") +
      blankTheme
    a <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~ErrorModel)
    
    p <- ggplot(d, aes(x=Value)) +
    xlab(paste(cov,"Covariate")) +
    ylab("Number of Observations") +
    blankTheme
    d <- p + geom_histogram(aes(fill=Recalibration,weight=nBases),alpha=0.6,binwidth=1,position="identity") + scale_fill_manual(values=c("maroon1","blue")) + facet_grid(.~ErrorModel) +     
    scale_y_continuous(formatter="comma")

  }
}

pdf("NA12878.20GAV.1.bqsr.recal.pdf",height=9,width=15)
distributeGraphRows(list(a,b,c), c(1,1,1))
distributeGraphRows(list(d,e,f), c(1,1,1))
dev.off()