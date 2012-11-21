#!/bin/tcsh
# 
# must be run from updated GATK source directory (looks for private/testdata/na12878db)

set args = ""
#set source = ~/Desktop/broadLocal/localData/na12878.db/
set source = /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase

#if ( ! -e $source ) then
#endif

if ( $1 == "local" ) then
set args = " -useLocal"
shift
endif

if ( $1 == "dev" ) then
set args = " -dbToUse DEV"
shift
endif


#set loc = 20:10,000,000-10,010,000
#set loc = 20:1-30,000,000
set loc = 20:10,000,000-11,000,000
#set loc = 20:1-10,009,259
#set loc = 20:1-1,000,000
#set loc = 20:10019093
#set loc = 20

set root = "java -Xmx2g -jar dist/GenomeAnalysisTK.jar -R $BUNDLE/b37/human_g1k_v37.fasta -L $loc $args"

# import all callsets
#   - must enumerate each call set individually from $source
if ( $1 == 1 ) then
set import = "$root -T ImportCallset"
$import -reset -callSetName Mills_1000G_GS_indels -assumedCallTruth UNKNOWN -howToTreatFilteredSites SKIP -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $source/Mills_and_1000G_gold_standard.indels.b37.na12878.20.vcf
$import -callSetName OMNI2.5Poly -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $source/Omni25_genotypes_2141_samples.b37.na12878.20.vcf
$import -callSetName OMNI2.5Mono -assumedCallTruth FALSE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V $source/Omni25_monomorphic_2141_samples.b37.na12878.20.vcf
$import -callSetName CEUTrio_best_practices -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $source/CEUTrio.HiSeq.WGS.b37.bestPractices.phased.b37.na12878.20.vcf
$import -callSetName 1000GPilot1Liftover -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 SKIP -V $source/ceu_yri_trios.genotypes.b37.na12878.20.vcf
$import -callSetName HapMap3.3 -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $source/genotypes_r27_nr.b37_fwd.na12878.20.vcf
$import -callSetName GoldIndelGenotyped -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $source/NA12878.indel.filtered.truth.na12878.20.vcf
endif

# import reviews from private testdata
foreach review (private/testdata/na12878kb/reviews.vcf)
$root -T ImportReviews -V $review
end

# update review
$root -T UpdateConsensus

