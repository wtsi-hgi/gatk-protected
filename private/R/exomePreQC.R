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

trim_to_95_pct <- function(data,complete_column,sampled_column) {
  mean <- mean(complete_column,na.rm=T)
  sd <- sd(complete_column,na.rm=T)
  min <- mean - 2*sd
  max <- mean + 2*sd
  return(subset(data,sampled_column>=min&sampled_column<=max))
}

add_distribution <- function(plot,column) {
  confidence_interval <- t.test(column)$conf.int
  mean <- mean(column,na.rm=T)
  sd <- sd(column,na.rm=T)
  text=c('-2*sigma','-1*sigma','mean','1*sigma','2*sigma')
  value=c(mean-2*sd,mean-1*sd,mean,mean+1*sd,mean+2*sd)
  plot <- plot + geom_hline(yintercept=value,color='blue',alpha=I(1/10))
  plot <- plot + geom_text(aes(x,y,label=label),data=data.frame(x=levels(complete$sample)[1],y=value,label=text),hjust=0,vjust=0,color='blue',size=2)
  return(plot)
}

# Return the number of months in the difference date-reference
number_of_months <- function(date,reference) {
  date_posix <- as.POSIXlt(date)
  reference_posix <- as.POSIXlt(reference)
  return((date_posix$year*12+date_posix$mon)-(reference_posix$year*12+reference_posix$mon))
}

# Load the complete dataset and add in a column to indicate sample 'goodness'; fill in with false for now, until we find a way to indicate sample quality
complete <- read.table(reference_dataset,header=T)
complete <- cbind(complete,good=F)
complete <- data.frame(complete,months_to_current_project=number_of_months(as.Date(complete$Last_Sequenced_WR_Created_Date),min(as.Date(data$Last_Sequenced_WR_Created_Date),na.rm=T)))
complete <- subset(complete,months_to_current_project>=-6)

# provide a reordering of the samples based on Last_Sequenced_WR
samples <- unique(rbind(data.frame(sample=paste(as.character(complete$sample),as.character(complete$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=complete$Last_Sequenced_WR_Created_Date),
			data.frame(sample=paste(as.character(data$sample),as.character(data$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=data$Last_Sequenced_WR_Created_Date)))
complete$sample <- factor(paste(complete$sample,complete$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR)])
data$sample <- factor(paste(data$sample,data$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR,samples$sample)])

#novel <- subset(complete,exon_intervals == "whole_exome_agilent_1.1_refseq_plus_3_boosters"&Novelty=="novel"&FunctionalClass=="all")
novel_sampled <- subset(complete,Novelty=="novel"&FunctionalClass=="all")
novel_sampled <- trim_to_95_pct(novel_sampled,complete$PCT_SELECTED_BASES,novel_sampled$PCT_SELECTED_BASES)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$MEAN_TARGET_COVERAGE,novel_sampled$MEAN_TARGET_COVERAGE)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$ZERO_CVG_TARGETS_PCT,novel_sampled$ZERO_CVG_TARGETS_PCT)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$PCT_TARGET_BASES_20X,novel_sampled$PCT_TARGET_BASES_20X)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$PCT_PF_READS,novel_sampled$PCT_PF_READS)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$PF_HQ_ERROR_RATE,novel_sampled$PF_HQ_ERROR_RATE)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$PF_INDEL_RATE,novel_sampled$PF_INDEL_RATE)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$MEAN_READ_LENGTH,novel_sampled$MEAN_READ_LENGTH)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$BAD_CYCLES,novel_sampled$BAD_CYCLES)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$STRAND_BALANCE,novel_sampled$STRAND_BALANCE)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$PCT_CHIMERAS,novel_sampled$PCT_CHIMERAS)
novel_sampled <- trim_to_95_pct(novel_sampled,complete$PCT_ADAPTER,novel_sampled$PCT_ADAPTER)

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

p <- ggplot(novel_sampled,aes(sample,PCT_SELECTED_BASES)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_SELECTED_BASES),color='red') + opts(title='% Selected Bases per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$PCT_SELECTED_BASES)
p

p <- ggplot(novel_sampled,aes(sample,MEAN_TARGET_COVERAGE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,MEAN_TARGET_COVERAGE),color='red') + opts(title='Mean Target Coverage per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$MEAN_TARGET_COVERAGE)
p

p <- ggplot(novel_sampled,aes(sample,ZERO_CVG_TARGETS_PCT)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,ZERO_CVG_TARGETS_PCT),color='red') + opts(title='% of Targets with <2x Coverage per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$ZERO_CVG_TARGETS_PCT)
p

p <- ggplot(novel_sampled,aes(sample,PF_INDEL_RATE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PF_INDEL_RATE),color='red') + opts(title='Indels per PF Read by Smaple') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$PF_INDEL_RATE)
p

p <- ggplot(novel_sampled,aes(sample,PCT_TARGET_BASES_20X)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_TARGET_BASES_20X),color='red') + opts(title='% Target Bases Achieving >20x Coverage per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$PCT_TARGET_BASES_20X)
p

p <- ggplot(novel_sampled,aes(sample,PCT_PF_READS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_PF_READS),color='red') + opts(title='% PF Reads Aligned per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$PCT_PF_READS)
p

p <- ggplot(novel_sampled,aes(sample,PF_HQ_ERROR_RATE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PF_HQ_ERROR_RATE),color='red') + opts(title='% HQ Bases mismatching the Reference per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$PF_HQ_ERROR_RATE)
p

p <- ggplot(novel_sampled,aes(sample,MEAN_READ_LENGTH)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100))+ geom_point(data=data,aes(sample,MEAN_READ_LENGTH),color='red') + opts(title='Mean Read Length per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$MEAN_READ_LENGTH)
p

p <- ggplot(novel_sampled,aes(sample,BAD_CYCLES)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,BAD_CYCLES),color='red') + opts(title='# Bad Cycles per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$BAD_CYCLES)
p

p <- ggplot(novel_sampled,aes(sample,STRAND_BALANCE)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,STRAND_BALANCE),color='red') + opts(title='% PF Reads Aligned to the + Strand per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$STRAND_BALANCE)
p

p <- ggplot(novel_sampled,aes(sample,TOTAL_SNPS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,TOTAL_SNPS),color='red') + opts(title='# SNPs called per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$TOTAL_SNPS)
p

pct_dbsnp_scatter <- ggplot(novel_sampled,aes(sample,PCT_DBSNP)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_DBSNP),color='red') + opts(title='% SNPs in dbSNP per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
pct_dbsnp_scatter <- add_distribution(p,complete$PCT_DBSNP)

pct_dbsnp_merged <- rbind(data.frame(sample=novel_sampled$sample,PCT_DBSNP=novel_sampled$PCT_DBSNP,type='reference'),data.frame(sample=data$sample,PCT_DBSNP=data$PCT_DBSNP,type='new'))
pct_dbsnp_density <- qplot(PCT_DBSNP,data=pct_dbsnp_merged,geom="density",fill=type) + opts(title='% SNPs in dbSNP per Sample')

layout <- grid.layout(nrow=2,ncol=1,heights=c(2,1))
grid.newpage()
pushViewport(viewport(layout=layout))
subplot <- function(x) viewport(layout.pos.row=x,layout.pos.col=1)
print(pct_dbsnp_scatter,vp=subplot(1))
print(pct_dbsnp_density,vp=subplot(2))

median_insert_size = rbind(data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='reference'),
		           data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='new'),
			   data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='reference'),
			   data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='new'),
			   data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='reference'),
			   data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='new'))
ggplot(median_insert_size,aes(sample,MEDIAN_INSERT_SIZE,color=data_type)) + geom_point() + geom_rug(aes(x=NULL),alpha=I(1/100)) + facet_grid(insert_type ~ .) + opts(title='Median Insert Size per Sample') + formatting + xlab('Sample (ordered by sequencing date)')

#median_insert_size_ref = rbind(data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='reference'),
#			       data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='reference'),
#			       data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='reference'))
#median_insert_size_new = rbind(data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_RF,insert_type='RF',data_type='new'),
#			       data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_FR,insert_type='FR',data_type='new'),
#			       data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM',data_type='new'))
#ggplot(median_insert_size_ref,aes(sample,MEDIAN_INSERT_SIZE,color=data_type,alpha=I(1/10))) + geom_point() + geom_point(data=median_insert_size_new,aes(sample,MEDIAN_INSERT_SIZE),color='red') + facet_grid(insert_type ~ .)

p <- ggplot(novel_sampled,aes(sample,PCT_CHIMERAS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_CHIMERAS),color='red') + opts(title='% Chimera Read Pairs per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$PCT_CHIMERAS)
p

p <- ggplot(novel_sampled,aes(sample,PCT_ADAPTER)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,PCT_ADAPTER),color='red') + opts(title='% Unaligned Reads Matching an Adapter Sequence per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$PCT_ADAPTER)
p

p <- ggplot(novel_sampled,aes(sample,NOVEL_SNPS)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,NOVEL_SNPS),color='red') + opts(title='# Novel SNPs called per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$NOVEL_SNPS)
p

p <- ggplot(novel_sampled,aes(sample,DBSNP_TITV)) + geom_point(alpha=I(1/10)) + geom_rug(aes(x=NULL),alpha=I(1/100)) + geom_point(data=data,aes(sample,DBSNP_TITV),color='red') + opts(title='TiTv of SNPs in dbSNP per Sample') + formatting + xlab('Sample (ordered by sequencing date)')
p <- add_distribution(p,complete$DBSNP_TITV)
p

if(onCMDLine) {
    dev.off()
}

