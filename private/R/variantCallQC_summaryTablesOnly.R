library(gplots)
library(gsalib)
library(ggplot2)
library(reshape)

args <- commandArgs(TRUE)

projectName <- args[1]
bySampleEval <- args[2]
outputPDF <- args[3]    

theme_set(theme_bw())

print("Report")
print(paste("Project          :", projectName))
print(paste("bySampleEval     :", bySampleEval))
print(paste("outputPDF        :", outputPDF))

expandVEReport <- function(d) {
  d$TiTvVariantEvaluator$tiTvRatio <- round(d$TiTvVariantEvaluator$tiTvRatio,2) 
  d$CountVariants$insertionDeletionRatio <- round(d$CountVariants$insertionDeletionRatio,2) 
  d$CountVariants$nIndels <- d$CountVariants$nInsertions + d$CountVariants$nDeletions
  d$CountVariants$nVariants <- d$CountVariants$nIndels + d$CountVariants$nSNPs + d$CountVariants$nMNPs 
  return(d)
}

getSize <- function(d) {
	x <- d$CountVariants$nProcessedLoci[[1]]
	return(x)
}

# Filters out any metric that isnt cumulative for the specified column name
selectCumulativeMetrics <- function(data,column) {
  lapply(data,function(x) { subset(x,x[column]=='all') })
}

# create sites section
summaryTable <- function(reportMetrics) { 
  
  allSmaples <- subset(reportMetrics, Sample == "all") 
  raw <- melt(allSmaples, id.vars=c("Novelty"), measure.vars=c("nProcessedLoci", "nVariants", "nSNPs", "tiTvRatio", "nMNPs", "nIndels", "insertionDeletionRatio"))
   
  table <- cast(raw, Novelty ~ ...,sum)
  
  # doesn't work with textplot
  colnames(table) <- c("Novelty", "Size (bp)", "Variants", "SNPs", "Ti/Tv", "MNPs", "Indels", "Ins/Del")
  return(table)
}

#create the per-sample section
sampleSummaryTable <- function(reportMetrics){ 
  
  metricsBySamples <- subset(reportMetrics, Sample != "all")
  raw <- melt(metricsBySamples, id.vars=c("Novelty", "Sample"), measure.vars=c("nProcessedLoci", "nVariants", "nSNPs", "tiTvRatio", "nMNPs", "nIndels","insertionDeletionRatio"))
  
  table <- cast(raw, Novelty ~ variable, mean)
  table$nMNPs <- round(table$nMNPs, 0)
  table$nSNPs <- round(table$nSNPs, 0)
  table$nIndels <- round(table$nIndels, 0)
  table$tiTvRatio <- round(table$tiTvRatio, 2)
  table$insertionDeletionRatio <- round(table$insertionDeletionRatio, 2)
  table$nVariants <- round(table$nVariants, 0)
  
  colnames(table) <- c("Novelty", "Size (bp)", "Variants", "SNPs", "Ti/Tv", "MNPs", "Indels", "Ins/Del")
  return(table)
}

overallSummaryTable <- function(reportMetrics, size){ 
  sitesSummary <- as.data.frame(summaryTable(reportMetrics))
  sitesSummary$"Metric Type" <- "Sites"
  
  sampleSummary <- as.data.frame(sampleSummaryTable(reportMetrics))
  sampleSummary$"Metric Type" <- "Per-sample"
  
  #create the expected values raw as a seperate table
  nExSNPs <- round(size/1000, 0)
  nExMNPs <- "?"
  nExIndels <- round(size/10000, 0)
  ExTiTvRatio <- "2.1 - 2.3"
  nExVariants <- "?"
  ExInsertionDeletionRatio <- "?"
  expectedLine <- matrix(c("all",size,nExVariants,nExSNPs,ExTiTvRatio,nExMNPs,nExIndels,ExInsertionDeletionRatio,"expected values per sample"),nrow = 1)	
  expextedLine <- as.table(expectedLine)
  colnames(expectedLine) <- c("Novelty", "Size (bp)", "Variants", "SNPs", "Ti/Tv", "MNPs", "Indels", "Ins/Del","Metric Type")
  
  # that last item puts the metric.type second in the list
  table <- rbind(sitesSummary, expectedLine, sampleSummary)[, c(1,9,2,3,4,5,6,7,8)]
  
  # remove columns with all NA/NaN
  table <- table <- table[,colSums(is.na(table))<nrow(table)]
  
  #transpose the table
  return(as.data.frame(t(as.matrix(table))))
}

CompSummaryTP <- function(gatkReport) {
	
	report <- selectCumulativeMetrics(gatkReport,'FunctionalClass')
	validationReport <-report$ValidationReport
	
	omni <- subset(validationReport, CompRod=="omni")
    omni_all <- CompSummaryAllSamplesTPByComp(omni)
    omni_all$"Comp" <- "Omni"
    
    gsIndels <- subset(validationReport, CompRod=="GSindels")
    gsIndels_all <- CompSummaryAllSamplesTPByComp(gsIndels)
    gsIndels_all$"Comp" <- "GS-Indels"
    
    table_allSamples <- rbind(omni_all, gsIndels_all)
    table <- table_allSamples[, c(1,7,6)]
    colnames(table) <- c("Novelty", "Comp", "overlap")
    table <- subset(table, Novelty != "novel")
     
    return(table)
}


CompSummaryFP <- function(gatkReport) {
	
	report <- selectCumulativeMetrics(gatkReport,'FunctionalClass')
	validationReport <- report$ValidationReport
    
    omni_mono <- subset(validationReport, CompRod=="omni_mono")
    omni_mono_all <- CompSummaryAllSamplesFPByComp(omni_mono)
    omni_mono_all$"Comp" <- "FP - omni mono"
    
    fp_MVL <- subset(validationReport, CompRod=="fp_MVL")
    fp_MVL_all <- CompSummaryAllSamplesFPByComp(fp_MVL)
    fp_MVL_all$"Comp" <- "FP MVL"
    
    table_allSamples <- rbind(omni_mono_all, fp_MVL_all)
    table <- table_allSamples[, c(1,7,6)]
    
    colnames(table) <- c("Novelty", "Comp","1-specificity")
    table <- subset(table, Novelty != "novel")

	return(table)
}


CompSummaryAllSamplesTPByComp <- function(compTable) {    
	comp <- subset(compTable, Sample == "all")
    raw <- melt(comp, id.vars=c("Novelty"), measure.vars=c("nComp", "TP", "FN", "sensitivity"))
    table <- cast(raw, Novelty ~ ...)
    table$sensitivity <- round(table$sensitivity, 2)
    table$sensitivityCalc <- apply(table,1,function(raw) {y <- paste(as.character(raw[4]) ,"(", as.character(raw[2]),"/",as.character(raw[3]+raw[2]),")", sep = " ")})
    return(table) 
}

CompSummaryAllSamplesFPByComp <- function(compTable) {    
	comp <- subset(compTable, Sample == "all")
    raw <- melt(comp, id.vars=c("Novelty"), measure.vars=c("nComp", "TP","FN","sensitivity"))
    table <- cast(raw, Novelty ~ ...)
    table$sensitivity <- round(table$sensitivity, 2)
    table$OneMinusSpecificity <- apply(table,1,function(raw) {y <- paste(as.character(raw[4]) ,"(", as.character(raw[2]),"/",as.character(raw[3]+raw[2]),")", sep = " ")})
    return(table) 
}

createReportMetrics <- function(gatkReport) {
  reportMetrics <- expandVEReport(selectCumulativeMetrics(gatkReport,'FunctionalClass'))
  r <- merge(reportMetrics$TiTvVariantEvaluator, reportMetrics$CountVariants)
  r <- merge(r, reportMetrics$CompOverlap)
  r <- subset(r, CompRod=="dbsnp")
  x <- subset(r, Novelty=="all")
  # order the samples by nSNPs -- it's the natural ordering.
  r$Sample <- factor(x$Sample, levels=x$Sample[order(x$nSNPs)])
  return(r) 
}

# -------------------------------------------------------
# Actually invoke the above functions 
# -------------------------------------------------------

# load the data.
gatkReport <- gsa.read.gatkreport(bySampleEval)
reportMetrics <- createReportMetrics(gatkReport)

size <- getSize(selectCumulativeMetrics(gatkReport,'Sample'))
openPDF(outputPDF)
cols <- c("black", "black", "black","red","black", "black", "black")

# Table of overall counts and quality
textplot(overallSummaryTable (reportMetrics, size), show.colnames=F,col.data=matrix(cols, nrow=9, byrow=T, ncol=7))
title(paste("Summary metrics for project", projectName), cex=3)

# comp tables
table1 <- as.data.frame(CompSummaryTP(gatkReport))
table2 <- as.data.frame(CompSummaryFP(gatkReport))
textplot(table1 , show.rownames=F) 
title(paste("TP Comp metrics for project", projectName), cex=3)

textplot(table2 , show.rownames=F) 
title(paste("FP Comp metrics for project", projectName), cex=3)

closePDF(outputPDF)
