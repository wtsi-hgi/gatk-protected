#
# Requirements:
#   Access to the reference data in the exomePreQC database, available in /humgen/gsa-scr1/GATK_Data.
#   A tsv file generated for the current project via $STING_HOME/private/python/generate_per_sample_metrics.py.
#
# To run:
#   Rscript exomePreQC.R <input tsv> <output pdf>
#
args = commandArgs(TRUE)
onCMDLine = ! is.na(args[1])

if ( onCMDLine ) {
  reference_dataset = '/humgen/gsa-scr1/GATK_Data/preqc.database'
  inputTSV = args[1]
  outputPDF = args[2]
} else {
  reference_dataset = '/humgen/gsa-scr1/GATK_Data/preqc.database'
  inputTSV = 'GWASeq_batches.metrics'
  outputPDF = 'GWASeq_batches.pdf'
}

require('ggplot2')
require('gplots')

data <- read.table(inputTSV,header=T)

trim_to_95_pct <- function(column) {
  mean <- mean(column,na.rm=T)
  sd <- sd(column,na.rm=T)
  if(sd > 0) {
      min <- mean - 2*sd
      max <- mean + 2*sd
      not_within_bounds <- function(value) {
        return(as.numeric(value<=min|value>=max))
      }
      return(sapply(column,not_within_bounds))
  }
  else {
      # Dataset was completely uniform.  Do not attempt to trim outliers.
      return(rep(0,length(column)))
  }
}

create_base_plot <- function(title,reference_dataset,new_dataset,column_name,include_sigmas=T) {
  # Create the basic plot with reference data
  p <- ggplot(reference_dataset,aes_string(x='sample',y=column_name))
  p <- p + opts(title=title,axis.ticks=theme_blank(),axis.text.x=theme_blank(),panel.grid.major=theme_blank(),panel.background=theme_blank())
  p <- p + xlab('Sample (ordered by sequencing date)')
  p <- p + geom_point(alpha=0.1) + geom_rug(aes(x=NULL),alpha=0.01)

  # Add in the new points and color them red.
  p <- p + geom_text(data=new_dataset,aes_string(x='sample',y=column_name,label='sample'),color='red',size=2)

  if(include_sigmas) {
    # Lines for the mean,+/- one sigma,+/- two sigma
    mean <- mean(complete[,column_name],na.rm=T)
    sd <- sd(complete[,column_name],na.rm=T)
    text=c('-2*sigma','-1*sigma','mean','1*sigma','2*sigma')
    value=c(mean-2*sd,mean-1*sd,mean,mean+1*sd,mean+2*sd)
    p <- p + geom_hline(yintercept=value,color='blue',alpha=0.1)
    p <- p + geom_text(aes(x,y,label=label),data=data.frame(x=levels(samples$sample)[1],y=value,label=text),hjust=0,vjust=0,color='blue',size=2)
  }

  return(p)
}

create_density_plot <- function(reference_dataset,new_dataset,column_name) {
  merged_dataset <- rbind(data.frame(sample=reference_dataset$sample,data=reference_dataset[,column_name],type='reference'),data.frame(sample=new_dataset$sample,data=new_dataset[,column_name],type=new_dataset$INITIATIVE))
  density_plot <- ggplot(merged_dataset,aes_string(x='data',y='..density..',fill='type')) + geom_density(na.rm=T) + ylab(column_name)
  density_plot <- density_plot + scale_fill_manual(values=c(alpha('black',0.1),alpha('red',0.5))) 
  density_plot <- density_plot + xlab(NULL) + ylab(NULL) + opts(legend.text=theme_text(size=6)) + labs(x=NULL,y=NULL)
  return(density_plot)
}

create_stock_plots <- function(title,reference_dataset,new_dataset,column_name) {
  base_plot <- create_base_plot(title,reference_dataset,new_dataset,column_name)
  density_plot <- create_density_plot(reference_dataset,new_dataset,column_name)

  layout <- grid.layout(nrow=2,ncol=1,heights=c(2,1))
  grid.newpage()
  pushViewport(viewport(layout=layout))
  subplot <- function(x) viewport(layout.pos.row=x,layout.pos.col=1)
  print(base_plot,vp=subplot(1))
  print(density_plot,vp=subplot(2))
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
complete <- subset(complete,months_to_current_project>=-12)

novel_sampled <- subset(complete,Novelty=="novel"&FunctionalClass=="all")
comparing_to_all = TRUE
if(any(novel_sampled$BAIT_SET %in% data$BAIT_SET)) {
    novel_sampled <- novel_sampled[novel_sampled$BAIT_SET %in% data$BAIT_SET,]
    comparing_to_all = FALSE
}

violations <- trim_to_95_pct(novel_sampled$PCT_SELECTED_BASES)
violations <- violations + trim_to_95_pct(novel_sampled$MEAN_TARGET_COVERAGE)
violations <- violations + trim_to_95_pct(novel_sampled$ZERO_CVG_TARGETS_PCT)
violations <- violations + trim_to_95_pct(novel_sampled$PCT_TARGET_BASES_20X)
violations <- violations + trim_to_95_pct(novel_sampled$PCT_PF_READS)
violations <- violations + trim_to_95_pct(novel_sampled$PF_HQ_ERROR_RATE)
violations <- violations + trim_to_95_pct(novel_sampled$PF_INDEL_RATE)
violations <- violations + trim_to_95_pct(novel_sampled$MEAN_READ_LENGTH)
violations <- violations + trim_to_95_pct(novel_sampled$BAD_CYCLES)
violations <- violations + trim_to_95_pct(novel_sampled$STRAND_BALANCE)
violations <- violations + trim_to_95_pct(novel_sampled$PCT_CHIMERAS)
violations <- violations + trim_to_95_pct(novel_sampled$PCT_ADAPTER)
novel_sampled <- subset(data.frame(novel_sampled,violations=violations),violations==0)	

# provide a reordering of the samples based on Last_Sequenced_WR
samples <- unique(rbind(data.frame(sample=paste(as.character(novel_sampled$sample),as.character(novel_sampled$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=novel_sampled$Last_Sequenced_WR_Created_Date),
			data.frame(sample=paste(as.character(data$sample),as.character(data$Last_Sequenced_WR_Created_Date)),Last_Sequenced_WR_Created_Date=data$Last_Sequenced_WR_Created_Date)))
novel_sampled$sample <- factor(paste(novel_sampled$sample,novel_sampled$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR)])
data$sample <- factor(paste(data$sample,data$Last_Sequenced_WR_Created_Date),levels=samples$sample[order(samples$Last_Sequenced_WR,samples$sample)])

# Write to PDF as necessary.
if(onCMDLine) {
    pdf(outputPDF)
}

# Specify a project header.
initiative <- as.character(unique(data$INITIATIVE))
num_samples <- length(unique(data$sample))
bait_set <- as.character(unique(data$BAIT_SET))
if(comparing_to_all) {
    bait_set <- paste(bait_set,'(Comparing custom bait set to all prior runs)',sep=' ')
}
total_reads <- sum(as.numeric(data$TOTAL_READS))
total_pf_reads <- sum(as.numeric(data$PF_READS))
total_pf_aligned_reads <- sum(as.numeric(data$PF_READS_ALIGNED))

num_reference_samples <- length(unique(complete$sample))
num_curated_reference_samples <- length(unique(novel_sampled$sample))

summary <- data.frame(keys=c('Initiative:','Number of Samples:','Bait Set:','Total Reads:','PF Reads:','PF Reads Aligned:','','Number of Samples in Reference Database:','Number of Samples in + Curated Database:'),values=c(initiative,num_samples,bait_set,total_reads,total_pf_reads,total_pf_aligned_reads,'',num_reference_samples,num_curated_reference_samples))
textplot(summary,show.rownames=F,show.colnames=F,valign=c("top"))
title('Project Summary Metrics')

fingerprint_samples = c()
fingerprint_lod_values = c()
fingerprint_lod_median = c()
for(i in 1:nrow(data)) {
   fingerprint_lods_for_sample <- eval(parse(text=data$FINGERPRINT_LODS[i]))
   # No fingerprint data for this sample?  Drop in an NA so that the fingerprint_samples database actually has reasonable values.
   if(is.null(fingerprint_lods_for_sample)) {
     fingerprint_lods_for_sample = c(NA)
   }
   fingerprint_samples <- c(fingerprint_samples,rep(as.character(data$sample[i]),length(fingerprint_lods_for_sample)))
   fingerprint_lod_values = c(fingerprint_lod_values,fingerprint_lods_for_sample)
   fingerprint_lod_median = c(fingerprint_lod_median,rep(median(fingerprint_lods_for_sample),length(fingerprint_lods_for_sample)))
}
fingerprint_lods = data.frame(sample=fingerprint_samples,median=fingerprint_lod_median,FINGERPRINT_LODS=fingerprint_lod_values)
fingerprint_lods$sample = factor(fingerprint_lods$sample,levels=unique(fingerprint_lods$sample[order(fingerprint_lods$median)]))

qplot(sample,FINGERPRINT_LODS,data=fingerprint_lods,geom="boxplot",outlier.size=0,main='Fingerprint LOD Scores By Sample') + opts(axis.text.x = theme_text(angle = 90,size=7)) + xlab('Sample') + ylab('LOD Score Distribution')

create_stock_plots('% Selected Bases per Sample',novel_sampled,data,'PCT_SELECTED_BASES')

create_stock_plots('Mean Target Coverage per Sample',novel_sampled,data,'MEAN_TARGET_COVERAGE')

create_stock_plots('% of Targets with <2x Coverage per Sample',novel_sampled,data,'ZERO_CVG_TARGETS_PCT')

create_stock_plots('# of Indels per PF Read by Sample',novel_sampled,data,'PF_INDEL_RATE')

create_stock_plots('% Target Bases Achieving >20x Coverage per Sample',novel_sampled,data,'PCT_TARGET_BASES_20X')

create_stock_plots('% PF Reads Aligned per Sample',novel_sampled,data,'PCT_PF_READS')

create_stock_plots('% HQ Bases mismatching the Reference per Sample',novel_sampled,data,'PF_HQ_ERROR_RATE')

create_base_plot('Mean Read Length per Sample',novel_sampled,data,'MEAN_READ_LENGTH')

create_base_plot('# Bad Cycles per Sample',novel_sampled,data,'BAD_CYCLES')

create_stock_plots('% PF Reads Aligned to the + Strand per Sample',novel_sampled,data,'STRAND_BALANCE')

create_stock_plots('# SNPs called per Sample',novel_sampled,data,'TOTAL_SNPS')

create_stock_plots('% SNPs in dbSNP per Sample',novel_sampled,data,'PCT_DBSNP')

median_insert_size_ref <- rbind(data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_RF,insert_type='RF'),		           
			        data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_FR,insert_type='FR'),
			        data.frame(sample=novel_sampled$sample,MEDIAN_INSERT_SIZE=novel_sampled$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM'))
median_insert_size_new <- rbind(data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_RF,insert_type='RF'),
			        data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_FR,insert_type='FR'),
			        data.frame(sample=data$sample,MEDIAN_INSERT_SIZE=data$MEDIAN_INSERT_SIZE_TANDEM,insert_type='TANDEM'))
create_base_plot('Median Insert Size per Sample',median_insert_size_ref,median_insert_size_new,'MEDIAN_INSERT_SIZE',include_sigmas=F) + facet_grid(insert_type ~ .)

create_stock_plots('% Chimera Read Pairs per Sample',novel_sampled,data,'PCT_CHIMERAS')

create_stock_plots('% Unaligned Reads Matching an Adapter Sequence per Sample',novel_sampled,data,'PCT_ADAPTER')

create_stock_plots('# Novel SNPs called per Sample',novel_sampled,data,'NOVEL_SNPS')

create_stock_plots('TiTv of SNPs in dbSNP per Sample',novel_sampled,data,'DBSNP_TITV')

if(onCMDLine) {
    dev.off()
}

