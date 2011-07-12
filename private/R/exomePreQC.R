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

# Load the complete dataset and add in a column to indicate sample 'goodness'; fill in with false for now, until we find a way to indicate sample quality
complete <- read.table(reference_dataset,header=T)
complete <- cbind(complete,good=F)

#novel <- subset(complete,exon_intervals == "whole_exome_agilent_1.1_refseq_plus_3_boosters"&Novelty=="novel"&FunctionalClass=="all")
novel <- subset(complete,Novelty=="novel"&FunctionalClass=="all")

# provide a reordering of the samples based on Last_Sequenced_WR
samples <- unique(rbind(data.frame(sample=paste(as.character(novel$sample),as.character(novel$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=novel$Last_Sequenced_WR_Created_Date),
			data.frame(sample=paste(as.character(data$sample),as.character(data$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=data$Last_Sequenced_WR_Created_Date)))
novel$sample <- factor(paste(novel$sample,novel$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR)])
data$sample <- factor(paste(data$sample,data$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR,samples$sample)])

samples = c()
fingerprint_lod_values = c()
fingerprint_lod_median = c()
for(i in 1:nrow(data)) {
   fingerprint_lods_for_sample <- eval(parse(text=data$FINGERPRINT_LODS[i]))
   samples <- c(samples,rep(data$sample[i],length(fingerprint_lods_for_sample)))
   fingerprint_lod_values = c(fingerprint_lod_values,fingerprint_lods_for_sample)
   fingerprint_lod_median = c(fingerprint_lod_median,rep(median(fingerprint_lods_for_sample),length(fingerprint_lods_for_sample)))
}
fingerprint_lods = data.frame(sample=samples,median=fingerprint_lod_median,FINGERPRINT_LODS=fingerprint_lod_values)
fingerprint_lods$sample = factor(fingerprint_lods$sample,levels=unique(fingerprint_lods$sample[order(fingerprint_lods$median)]))

if(onCMDLine) {
    pdf(outputPDF)
}

qplot(sample,FINGERPRINT_LODS,data=fingerprint_lods,geom="boxplot",outlier.size=0,main='Fingerprint LOD Scores By Sample') + opts(axis.text.x = theme_text(angle = 90,size=7)) + xlab('Sample') + ylab('LOD Score Distribution')

formatting = opts(axis.ticks=theme_blank(),axis.text.x=theme_blank(),panel.grid.major=theme_blank(),panel.background=theme_blank())

ggplot(novel,aes(sample,PCT_SELECTED_BASES)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_SELECTED_BASES),color='red') + opts(title='Mean Target Coverage per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,MEAN_TARGET_COVERAGE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,MEAN_TARGET_COVERAGE),color='red') + opts(title='Mean Target Coverage per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,ZERO_CVG_TARGETS_PCT)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,ZERO_CVG_TARGETS_PCT),color='red') + opts(title='% of Targets with <2x Coverage per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,FOLD_80_BASE_PENALTY)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,FOLD_80_BASE_PENALTY),color='red') + opts(title='Fold 80 Base Penalty per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,PCT_TARGET_BASES_20X)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_TARGET_BASES_20X),color='red') + opts(title='% Target Bases Achieving >20x Coverage per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

novel_sampling = subset(novel,(PCT_PF_READS<1.0&PCT_PF_READS>0.5)|good==T)
ggplot(novel_sampling,aes(sample,PCT_PF_READS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_PF_READS),color='red') + opts(title='% PF Reads Aligned per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

ggplot(novel,aes(sample,PF_HQ_ERROR_RATE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PF_HQ_ERROR_RATE),color='red') + opts(title='% HQ Bases mismatching the Reference per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,MEAN_READ_LENGTH)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100))+ geom_point(data=data,aes(sample,MEAN_READ_LENGTH),color='red') + opts(title='Mean Read Length per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,BAD_CYCLES)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,BAD_CYCLES),color='red') + opts(title='# Bad Cycles per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
#ggplot(novel,aes(sample,HS_LIBRARY_SIZE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,HS_LIBRARY_SIZE),color='red') + opts(title='TiTv of SNPs in dbSNP per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

novel_sampling = subset(novel,(STRAND_BALANCE>0.4&STRAND_BALANCE<0.6)|good==T)
ggplot(novel_sampling,aes(sample,STRAND_BALANCE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,STRAND_BALANCE),color='red') + opts(title='% PF Reads Aligned to the + Strand per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

novel_sampling = subset(novel,TOTAL_SNPS>15000|good==T)
ggplot(novel_sampling,aes(sample,TOTAL_SNPS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,TOTAL_SNPS),color='red') + opts(title='# SNPs called per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

novel_sampling = subset(novel,(PCT_DBSNP>0.8)|good==T)
pct_dbsnp_scatter <- ggplot(novel_sampling,aes(sample,PCT_DBSNP)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_DBSNP),color='red') + opts(title='% SNPs in dbSNP per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
pct_dbsnp_merged <- rbind(data.frame(sample=novel_sampling$sample,PCT_DBSNP=novel_sampling$PCT_DBSNP,type='reference'),data.frame(sample=data$sample,PCT_DBSNP=data$PCT_DBSNP,type='new'))
pct_dbsnp_density <- qplot(PCT_DBSNP,data=pct_dbsnp_merged,geom="density",fill=type) + opts(title='% SNPs in dbSNP per Sample')

layout <- grid.layout(nrow=2,ncol=1,heights=c(2,1))
grid.newpage()
pushViewport(viewport(layout=layout))
subplot <- function(x) viewport(layout.pos.row=x,layout.pos.col=1)
print(pct_dbsnp_scatter,vp=subplot(1))
print(pct_dbsnp_density,vp=subplot(2))

median_insert_size = rbind(data.frame(sample=novel$sample,MEDIAN_INSERT_SIZE=novel$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='reference'),
		           data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='new'),
			   data.frame(sample=novel$sample,MEDIAN_INSERT_SIZE=novel$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='reference'),
			   data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='new'),
			   data.frame(sample=novel$sample,MEDIAN_INSERT_SIZE=novel$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='reference'),
			   data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='new'))
ggplot(median_insert_size,aes(sample,MEDIAN_INSERT_SIZE,color=data_type)) + geom_point() + geom_rug(aes(x=NULL),alpha=I(1/100)) + facet_grid(insert_type ~ .) + opts(title='Median Insert Size per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

#median_insert_size_ref = rbind(data.frame(sample=novel$sample,MEDIAN_INSERT_SIZE=novel$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='reference'),
#			       data.frame(sample=novel$sample,MEDIAN_INSERT_SIZE=novel$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='reference'),
#			       data.frame(sample=novel$sample,MEDIAN_INSERT_SIZE=novel$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='reference'))
#median_insert_size_new = rbind(data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='new'),
#			       data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='new'),
#			       data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='new'))
#ggplot(median_insert_size_ref,aes(sample,MEDIAN_INSERT_SIZE,color=data_type,alpha=I(1/10))) + geom_point() + geom_point(data=median_insert_size_new,aes(sample,MEDIAN_INSERT_SIZE),color='red') + facet_grid(insert_type ~ .)

ggplot(novel,aes(sample,PCT_CHIMERAS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_CHIMERAS),color='red') + opts(title='% Chimera Read Pairs per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,PCT_ADAPTER)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_ADAPTER),color='red') + opts(title='% Unaligned Reads Matching an Adapter Sequence per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,NOVEL_SNPS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,NOVEL_SNPS),color='red') + opts(title='# Novel SNPs called per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
ggplot(novel,aes(sample,DBSNP_TITV)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,DBSNP_TITV),color='red') + opts(title='TiTv of SNPs in dbSNP per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

if(onCMDLine) {
    dev.off()
}

