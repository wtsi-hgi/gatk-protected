args = commandArgs(TRUE)
onCMDLine = ! is.na(args[1])

if ( onCMDLine ) {
  reference_dataset = '/Users/mhanna/src/StingUnstable/private/R/preqc.database'
  inputTSV = args[1]
  outputPDF = args[2]
} else {
  reference_dataset = '/Users/mhanna/src/StingUnstable/private/R/preqc.database'
  inputTSV = 'GoT2D_exomes_batch_005_per_sample_metrics.tsv'
  outputPDF = 'T2D.pdf'
}

require('ggplot2')

data <- read.table(inputTSV,header=T)

complete <- read.table(reference_dataset,header=T)
#novel <- subset(complete,exon_intervals == "whole_exome_agilent_1.1_refseq_plus_3_boosters"&Novelty=="novel"&FunctionalClass=="all")
novel <- subset(complete,Novelty=="novel"&FunctionalClass=="all")
selected_samples <- novel$sample %in% data$sample
novel_with_highlights <- cbind(novel,selected_samples)

# provide a reordering of the samples based on Last_Sequenced_WR
samples <- unique(rbind(data.frame(sample=paste(as.character(novel$sample),as.character(novel$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=novel$Last_Sequenced_WR_Created_Date),
			data.frame(sample=paste(as.character(data$sample),as.character(data$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=data$Last_Sequenced_WR_Created_Date)))
novel$sample <- factor(paste(novel$sample,novel$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR)])
data$sample <- factor(paste(data$sample,data$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR,samples$sample)])

fingerprint_lods = list()
for(i in 1:nrow(data)) {
  fingerprint_lods[[as.character(data$sample[i])]] <- eval(parse(text=data$FINGERPRINT_LODS[i]))
}

fingerprint_lod_order = order(unlist(lapply(fingerprint_lods,median),use.names=F))

if(onCMDLine) {
    pdf(outputPDF)
}

boxplot(fingerprint_lods[fingerprint_lod_order],las=3,main='Fingerprint LOD Scores By Sample',xlab='Sample',ylab='LOD Score Distribution',cex.axis=0.65)

x_axis_formatting = opts(axis.ticks=theme_blank(),axis.text.x=theme_blank())

ggplot(novel,aes(sample,PCT_SELECTED_BASES)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,PCT_SELECTED_BASES),color='blue') + opts(title='Mean Target Coverage per Sample') + x_axis_formatting
ggplot(novel,aes(sample,MEAN_TARGET_COVERAGE)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,MEAN_TARGET_COVERAGE),color='blue') + opts(title='Mean Target Coverage per Sample') + x_axis_formatting
ggplot(novel,aes(sample,ZERO_CVG_TARGETS_PCT)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,ZERO_CVG_TARGETS_PCT),color='blue') + opts(title='% of Targets with <2x Coverage per Sample')  + x_axis_formatting
ggplot(novel,aes(sample,FOLD_80_BASE_PENALTY)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,FOLD_80_BASE_PENALTY),color='blue') + opts(title='Fold 80 Base Penalty per Sample') + x_axis_formatting
ggplot(novel,aes(sample,PCT_TARGET_BASES_20X)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,PCT_TARGET_BASES_20X),color='blue') + opts(title='% Target Bases Achieving >20x Coverage per Sample') + x_axis_formatting
ggplot(novel,aes(sample,PCT_PF_READS)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,PCT_PF_READS),color='blue') + opts(title='% PF Reads Aligned per Sample') + x_axis_formatting
ggplot(novel,aes(sample,PF_HQ_ERROR_RATE)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,PF_HQ_ERROR_RATE),color='blue') + opts(title='% HQ Bases mismatching the Reference per Sample') + x_axis_formatting
ggplot(novel,aes(sample,MEAN_READ_LENGTH)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,MEAN_READ_LENGTH),color='blue') + opts(title='Mean Read Length per Sample') + x_axis_formatting
ggplot(novel,aes(sample,BAD_CYCLES)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,BAD_CYCLES),color='blue') + opts(title='# Bad Cycles per Sample') + x_axis_formatting
ggplot(novel,aes(sample,STRAND_BALANCE)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,STRAND_BALANCE),color='blue') + opts(title='% PF Reads Aligned to the + Strand per Sample') + x_axis_formatting
ggplot(novel,aes(sample,TOTAL_SNPS)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,TOTAL_SNPS),color='blue') + opts(title='# SNPs called per Sample') + x_axis_formatting
ggplot(novel,aes(sample,PCT_DBSNP)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,PCT_DBSNP),color='blue') + opts(title='% SNPs in dbSNP per Sample') + x_axis_formatting
qplot(PCT_DBSNP,data=data,geom="histogram") + opts(title='% SNPs in dbSNP per Sample') + x_axis_formatting
#ggplot(novel,aes(sample,MEDIAN_INSERT_SIZE_RF)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,MEDIAN_INSERT_SIZE_RF),color='blue') + opts(title='Median Insert Size per Sample (RF pairs)') + x_axis_formatting + scale_y_discrete()
#ggplot(novel,aes(sample,MEDIAN_INSERT_SIZE_FR)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,MEDIAN_INSERT_SIZE_FR),color='blue') + opts(title='Median Insert Size per Sample (FR pairs') + x_axis_formatting + scale_y_discrete()
#ggplot(novel,aes(sample,MEDIAN_INSERT_SIZE_TANDEM)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,MEDIAN_INSERT_SIZE_TANDEM),color='blue') + opts(title='Median Insert Size per Sample (tandem pairs') + x_axis_formatting + scale_y_discrete()
ggplot(novel,aes(sample,PCT_CHIMERAS)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,PCT_CHIMERAS),color='blue') + opts(title='% Chimera Read Pairs per Sample') + x_axis_formatting
ggplot(novel,aes(sample,PCT_ADAPTER)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,PCT_ADAPTER),color='blue') + opts(title='% Unaligned Reads Matching an Adapter Sequence per Sample') + x_axis_formatting
ggplot(novel,aes(sample,NOVEL_SNPS)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,NOVEL_SNPS),color='blue') + opts(title='# Novel SNPs called per Sample') + x_axis_formatting
ggplot(novel,aes(sample,DBSNP_TITV)) + geom_point(alpha=I(1/10)) + geom_point(data=data,aes(sample,DBSNP_TITV),color='blue') + opts(title='TiTv of SNPs in dbSNP per Sample') + x_axis_formatting

if(onCMDLine) {
    dev.off()
}

#qplot(sample,Library_Size_HS,data=novel_with_highlights,color=selected_samples) + opts(title='Hybrid Sequencing Library Size per Sample')
