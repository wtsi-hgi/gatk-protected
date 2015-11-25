#!/bin/tcsh
# 
# must be run from updated GATK source directory (looks for private/testdata/na12878db)

set args = ""
#set source = ~/Desktop/broadLocal/localData/na12878.db/
set source = /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase
set BUNDLE = /humgen/gsa-hpprojects/GATK/bundle/current/b37/
set data = /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/data/
set jarPath = target/GenomeAnalysisTK-internal.jar
#if ( ! -e $source ) then
#endif

if ( $1 == "local" ) then
set args = " -useLocal"
shift
endif

if ( $1 == "production" ) then
set args = "$args -dbToUse PRODUCTION"
shift
endif


#set loc = 20:10,000,000-10,010,000
#set loc = 20:1-30,000,000
#set loc = 20:10,000,000-11,000,000
#set loc = 20:1-10,009,259
#set loc = 20:1-1,000,000
#set loc = 20:10019093
#set loc = 20
#set root = "echo java -Xmx2g -jar dist/GenomeAnalysisTK.jar -R $BUNDLE/human_g1k_v37.fasta $args"

set gencode = /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/gencode.v12_broad.agilent_merged.50bpPadding.interval_list
set root = "java -Xmx2g -jar $jarPath -R $BUNDLE/human_g1k_v37.fasta $args"

# import all callsets
#   - must enumerate each call set individually from $source
if ( $1 != 1 ) then
set import = "$root -T ImportCallset"
$import -reset -callSetName Mills_1000G_GS_indels -assumedCallTruth UNKNOWN -howToTreatFilteredSites SKIP -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $BUNDLE/Mills_and_1000G_gold_standard.indels.b37.vcf
$import -callSetName OMNI2.5Poly -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_genotypes_2141_samples.b37.vcf
$import -callSetName OMNI2.5Mono -assumedCallTruth FALSE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_monomorphic_2141_samples.b37.vcf
$import -untrusted -callSetName CEUTrio_best_practices -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $BUNDLE/CEUTrio.HiSeq.WGS.b37.bestPractices.phased.b37.vcf
$import -callSetName 1000GPilot2Liftover -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 SKIP -V $data/ceu_yri_trios.genotypes.b37.vcf
$import -callSetName HapMap3.3 -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf
$import -callSetName GoldIndelGenotyped -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $data/NA12878.indel.filtered.truth.vcf
$import -callSetName 1000G_250sites_indelValidation_POLY -assumedCallTruth UNKNOWN -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V $source/1000G.250sites.indelValidation.polymorphic.alleles.na12878.kbInterval.vcf
$import -callSetName 1000G_250sites_indelValidation_MONO -assumedCallTruth FALSE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 FALSE_POSITIVE -V $source/1000G.250sites.indelValidation.monomorphic.alleles.na12878.kbInterval.vcf
$import -callSetName 1000G_exomeChip -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Illumina_ExomeChip/1000G.exomechip.20121009.snps_only.genotypes.vcf
$import -callSetName 1000G_variousValidations -assumedCallTruth UNKNOWN -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/1000G_validation_experiments/1000G.validationExperiments.polymorphic.b37.vcf
$import -callSetName 1000G_snpChip -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/1kg_chip_jan2010/1000Genome.chip.b37.filtered.vcf
$import -callSetName AffyAxiom -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Affymetrix_Axiom/Affymetrix_Axiom_DB_2010_v4_b37.vcf
$import -callSetName 1000G_exomeIndels -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V /humgen/gsa-hpprojects/dev/validationExperiments/exomeIndels/HaplotypeCallerGGA/exomeIndels.HC.GGA.MiSeq.chr20.raw.vcf
$import -callSetName MiSeqLargeIndels -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 FALSE_POSITIVE -V /humgen/gsa-hpprojects/dev/validationExperiments/largeIndelValidation/largeIndels.HC.GGA.MiSeq.chr20.raw.vcf
$import -callSetName largeScaleValidationPools_POLY -assumedCallTruth UNKNOWN -howToTreatFilteredSites SKIP -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/largeScaleValidationSites_pools_run20121030EASFix.POLYMORPHIC.sites.vcf
$import -callSetName largeScaleValidationPools_MONO -assumedCallTruth FALSE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 FALSE_POSITIVE -V /humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/largeScaleValidationSites_pools_run20121030EASFix.MONOMORPHIC.sites.vcf
$import -callSetName NIST_GenomesInABottle -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/NA12878Collection/NIST/v2.17_09112013/NISTIntegratedCalls_13datasets_130719_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.17_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf
$import -callSetName AffyExomePlus_POLY -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $source/Axiom_Exome_Plus.genotypes.NA12878.poly.vcf
$import -callSetName AffyExomePlus_MONO -assumedCallTruth FALSE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 FALSE_POSITIVE -V $source/Axiom_Exome_Plus.genotypes.all_populations.monomorphic.biallelic.vcf
$import -callSetName BroadExomeLOF_POLY -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 MARK_AS_NON_POLYMORPHIC -V $source/ALL.wex.broad_lof_exome_chip_beadstudio.20130703.snps_and_indels.chip.genotypes.NA12878.poly.vcf
$import -callSetName BroadExomeLOF_MONO -assumedCallTruth FALSE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 FALSE_POSITIVE -V $source/ALL.wex.broad_lof_exome_chip_beadstudio.20130703.snps_and_indels.chip.genotypes.mono.vcf
#$import -callSetName largeScaleValidationNA12878 -assumedCallTruth TRUE_POSITIVE -howToTreatFilteredSites FALSE_POSITIVE -howToTreatAC0 SKIP -V /humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/HaplotypeCallerAnalysis/largeScaleValidation.HC.NA12878.chr20.recalibrated.vcf
endif

# import reviews from private testdata
set reviewArchive = `ls -td /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/reviewsBackup/**.vcf | head -n 1`
foreach review (private/testdata/na12878kb/reviews.vcf private/testdata/na12878kb/omnipoly.reviews.vcf private/testdata/na12878kb/reviews.ebanks.vcf private/testdata/na12878kb/reviews.rpoplin.vcf private/testdata/na12878kb/reviews.justinzook.vcf $reviewArchive)
$root -T ImportReviews -V $review
end

# update review
$root -T UpdateConsensus

