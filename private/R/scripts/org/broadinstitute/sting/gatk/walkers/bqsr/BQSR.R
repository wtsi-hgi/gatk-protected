library("ggplot2")

args <- commandArgs(TRUE)
data <- read.csv(args[1])
data <- within(data, EventType <- factor(EventType, levels = rev(levels(EventType))))

numRG = length(unique(data$ReadGroup))
blankTheme = opts(panel.grid.major = theme_blank(), panel.grid.minor = theme_blank(), panel.background = theme_blank(), axis.ticks = theme_blank())

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


for(cov in levels(data$CovariateName)) {    # for each covariate in turn  
  d = data[data$CovariateName==cov,]        # pull out just the data for this covariate so we can treat the non-numeric values appropriately
  if( cov == "Context" ) {
    d = subset(d, grepl(pattern="AAAAA...", CovariateValue))
    d$CovariateValue = as.character(d$CovariateValue)
  } else {
    d$CovariateValue = as.numeric(levels(d$CovariateValue))[as.integer(d$CovariateValue)] # efficient way to convert factors back to their real values
  }
#  d=subset(d,Observations>2000) # only show bins which have enough data to acually estimate the quality
  
  if( cov != "QualityScore" ) {    
    p <- ggplot(d, aes(x=CovariateValue,y=Accuracy,alpha=log10(Observations))) +
      geom_abline(intercept=0, slope=0, linetype=2) + 
      xlab(paste(cov,"Covariate")) +
      ylab("Quality Score Accuracy") +
      blankTheme
    if(cov == "Cycle") {
      b <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~EventType) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))
      
      p <- ggplot(d, aes(x=CovariateValue,y=AverageReportedQuality,alpha=log10(Observations))) +
        xlab(paste(cov,"Covariate")) +
        ylab("Mean Quality Score") +
        blankTheme
      e <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~EventType) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))
      
      
    } else {
      c <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~EventType) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))
      p <- ggplot(d, aes(x=CovariateValue,y=AverageReportedQuality,alpha=log10(Observations))) +
        xlab(paste(cov,"Covariate")) +
        ylab("Mean Quality Score") +
        blankTheme
      f <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~EventType) +
        opts(axis.text.x=theme_text(angle=90, hjust=0))
      
    }
  } else {
    p <- ggplot(d, aes(x=AverageReportedQuality,y=EmpiricalQuality,alpha=log10(Observations))) +
      geom_abline(intercept=0, slope=1, linetype=2) + 
      xlab("Reported Quality Score") +
      ylab("Empirical Quality Score") +
      blankTheme
    a <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("maroon1","blue")) + facet_grid(.~EventType)
    
    p <- ggplot(d, aes(x=CovariateValue)) +
      xlab(paste(cov,"Covariate")) +
      ylab("Number of Observations") +
      blankTheme
    d <- p + geom_histogram(aes(fill=Recalibration,weight=Observations),alpha=0.6,binwidth=1,position="identity") + scale_fill_manual(values=c("maroon1","blue")) + facet_grid(.~EventType) +     
      scale_y_continuous(formatter="comma")
    
  }
}

pdf(args[2],height=9,width=15)
distributeGraphRows(list(a,b,c), c(1,1,1))
distributeGraphRows(list(d,e,f), c(1,1,1))
dev.off()