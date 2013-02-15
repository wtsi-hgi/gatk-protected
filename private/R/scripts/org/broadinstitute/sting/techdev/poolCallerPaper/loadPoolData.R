outputDir="/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/processedTables"

loadData <- T
if (loadData) {
  allPoolStats <- read.table(sprintf("%s/allPoolStats.table",outputDir),header=T)
  allPoolData <- read.table(sprintf("%s/allPoolData.table",outputDir),header=T)
  allStats <- read.table(sprintf("%s/allStats.table",outputDir),header=T)
  allData <- read.table(sprintf("%s/allData.table",outputDir),header=T)
  
}

print("number of LOF sites we have:") 
print(nrow(subset(allData,inLOF==T)))
print("SNPS:") 
print(nrow(subset(allData,inLOF==T & TYPE=="SNP")))

print("indels:")
print(nrow(subset(allData,inLOF==T& TYPE=="INDEL")))

print("SNPs and indels that passed filters in pool caller:") 
print(nrow(subset(allData,inLOF==T & TYPE=="SNP" & goodSite==T)))
print(nrow(subset(allData,inLOF==T& TYPE=="INDEL"& goodSite==T)))

print("SNPs and Indels that passed filters and validated in Pool caller as polymorphic:")
print(nrow(subset(allData,inLOF==T & TYPE=="SNP"& goodSite==T & AC.LOF > 0)))
print(nrow(subset(allData,inLOF==T & TYPE=="INDEL"& goodSite==T & AC.LOF > 0)))

## do the same for only HC LOF sites:
print("HC sites that passed filters in pool caller:")
print(nrow(subset(allData,inLOF==T & TYPE=="SNP" & goodSite==T &  LOF=="HC")))
print(nrow(subset(allData,inLOF==T & TYPE=="INDEL"& goodSite==T &  LOF=="HC")))
print("high confidence SNPs and Indels that passed filters and validated in Pool caller as polymorphic:")
print(nrow(subset(allData,inLOF==T & TYPE=="SNP" & goodSite==T & AC.LOF > 0 & LOF=="HC")))
print(nrow(subset(allData,inLOF==T & TYPE=="INDEL"& goodSite==T & AC.LOF > 0 & LOF=="HC")))
