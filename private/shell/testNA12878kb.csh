#!/bin/tcsh

set args = ""
if ( $1 == "local" ) then
set args = " -useLocal"
shift
endif

#set loc = 20:10,000,000-10,010,000
#set loc = 20:1-30,000,000
set loc = 20:10,000,000-11,000,000
#set loc = 20:10019093
#set loc = 20:1-10,009,259
#set loc = 20:1-1,000,000
#set loc = 20

set root = "java -Xmx2g -jar dist/GenomeAnalysisTK.jar -R ~/Desktop/broadLocal/localData/human_g1k_v37.fasta -L $loc $args"

if ( $1 == 1 ) then
$root -T ExportReviews -o export_reviews.vcf
endif

set assess = "$root -T AssessNA12878 -badSites falseNegativesAndPositives.vcf -L /Users/depristo/Desktop/broadLocal/localData/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list -isr INTERSECTION"

if ( $1 == 2 ) then
$assess -V /Users/depristo/Desktop/broadLocal/localData/CEUTrio.HiSeq.WGS.b37_decoy.recalibrated.vcf -AssessmentsToExclude CALLED_NOT_IN_DB_AT_ALL
endif

if ( $1 == 3.0 ) then
$assess -V hcSingleExome/hc.vcf 
endif

if ( $1 == 3.1 ) then
$assess -V hcSingleExome/ug.vcf 
endif

