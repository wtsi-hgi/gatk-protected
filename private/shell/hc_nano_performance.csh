#!/bin/tcsh

#set loc = "20:10,000,000-10,002,000"
#set loc = "20:10,000,000-10,010,000"
#set loc = "20:10,000,000-11,000,000"
set loc = "20:10,000,000-15,000,000"
set GATK = GATK/dist/GenomeAnalysisTK.jar
set performanceBAM = /dev/shm/tmp.bam
set BAM = /humgen/gsa-hpprojects/dev/depristo/oneOffProjects/ceuTrioBestPractices/ceuTrio.bwaaln.bam.list

set root = "java -Xmx4g -jar $GATK -T HaplotypeCaller -I $performanceBAM -R $BUNDLE/b37/human_g1k_v37.fasta -L $loc -o /dev/null"

# load the bam to the local cache
java -Xmx4g -jar $GATK -T PrintReads -I $BAM -R $BUNDLE/b37/human_g1k_v37.fasta -L $loc -o $performanceBAM

# serial calculations
foreach nct (32 24 16 12)
  $root -nct $nct -log hc_nano_${nct}_iteration_$1.log
end

foreach nct (16 12 8 7 6 5 4 3 2 1)
  $root -nct $nct -log hc_nano_${nct}_iteration_$1.log &
end
