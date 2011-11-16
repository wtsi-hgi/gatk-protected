package org.broadinstitute.sting.queue.qscripts.reducedreads

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.exceptions.UserException

class CompareCalls extends QScript {
  @Argument(shortName = "rb", doc = "Reduced BAM / BAM list", required=true) val reducedBAM: File = null
  @Argument(shortName = "fb", doc = "Full BAM / BAM list", required=true) val fullBAM: File = null
  @Argument(shortName = "g", doc = "Goldstandard callset from the Full BAMs", required=false) val truthCallset: File = null
  @Argument(shortName = "r", doc = "Reference sequence", required=false) val referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "i", doc = "Intervals file", required=false) val intervalsFile: File = new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")
  @Argument(shortName = "s", doc = "scatterCount", required=false) val scatterCount: Int = 400


  val hapmap = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val dbsnp = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf"
  val omni = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"
  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
  val projectConsensus_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf"
  val badSites_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"


  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.intervals :+= intervalsFile
    this.memoryLimit = 4
  }

  def getVCFName (bam: File) : File = {
    if (bam.endsWith(".bam")) {
      swapExt(bam, ".bam", ".vcf")
    } else if (bam.endsWith(".list")) {
      swapExt(bam, ".list", ".vcf")
    } else {
      throw new UserException("input file must be a BAM (ending with .bam) or a list of bam files (ending with .list)")
    }
  }

  def script = {

    call(reducedBAM)
    call(fullBAM)
  }

  def call(bam: File) {
    val callsVCF = getVCFName(bam)
    val outVCF = swapExt(callsVCF, ".vcf", ".recalibrated.vcf")
    val tranchesFile = swapExt(callsVCF, ".vcf", ".tranches")
    val recalFile = swapExt(callsVCF, ".vcf", ".recal_file")

    val ug = new UnifiedGenotyper() with UNIVERSAL_GATK_ARGS
    ug.input_file :+= bam
    ug.dbsnp = dbsnp
    ug.out = callsVCF
    ug.scatterCount = scatterCount
    ug.nt = 2
    ug.stand_call_conf = 30.0
    ug.stand_emit_conf = 30.0

    val vqsr = new VariantRecalibrator() with UNIVERSAL_GATK_ARGS
    vqsr.nt = 2
    vqsr.input :+= callsVCF
    vqsr.resource :+= new TaggedFile( hapmap, "training=true,truth=true,prior=15.0" )
    vqsr.resource :+= new TaggedFile( omni, "training=true,truth=true,prior=12.0" )
    vqsr.resource :+= new TaggedFile( training_1000G, "training=true,prior=10.0" )
    vqsr.resource :+= new TaggedFile( dbsnp, "known=true,prior=2.0" )
    vqsr.resource :+= new TaggedFile( projectConsensus_1000G, "prior=8.0" )
    vqsr.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "InbreedingCoeff")
    vqsr.resource :+= new TaggedFile( badSites_1000G, "bad=true,prior=2.0")
    vqsr.mG = 6
    vqsr.tranches_file = tranchesFile
    vqsr.recal_file = recalFile
    vqsr.allPoly = true
    vqsr.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    vqsr.rscriptFile = swapExt(callsVCF, ".vcf", ".vqsr.r")

    val cut = new ApplyRecalibration() with UNIVERSAL_GATK_ARGS
    cut.memoryLimit = 6
    cut.input :+= callsVCF
    cut.tranchesFile = tranchesFile
    cut.recalFile = recalFile
    cut.ts_filter_level = 98.0
    cut.out = outVCF

    add(ug, vqsr, cut)
  }

}
