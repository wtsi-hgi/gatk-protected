package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils


class BQSR_with_different_known_sets extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  val devDir = "/humgen/gsa-hpprojects/dev/ebanks/oneOffProjects/phaseIIpipeline/"

  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val lowQualBam = new File("/humgen/1kg/DCC/ftp/data/NA12878/alignment/NA12878.chrom1.ILLUMINA.bwa.CEU.high_coverage.20100311.bam")
  val highQualBam = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.bam")

  val dbSNP_129 = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf")
  val dbSNP_132 = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf")
  val v2b = new File("/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf")

  val omni = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"
  val hapmap = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf"

  val queueLogDir = ".qlog/"

  var nContigs:Int = 0

  def script = {
    nContigs = QScriptUtils.getNumberOfContigs(highQualBam)

    doWork("null", null, null);
    doWork("129", dbSNP_129, null);
    doWork("132", dbSNP_132, null);
    doWork("1000G", null, v2b);
    doWork("134", dbSNP_132, v2b);
  }

  def doWork(prefix: String, dbSNP: File, KG: File) {
    val high = new File(devDir + prefix + ".highQual.recal")
    add(new cov(highQualBam, high, dbSNP, KG))
    add(new recal(highQualBam, high, new File(devDir + prefix + ".highQual.bam")))

    val low = new File(devDir + prefix + ".lowQual.recal")
    add(new cov(lowQualBam, low, dbSNP, KG))
    add(new recal(lowQualBam, low, new File(devDir + prefix + ".lowQual.bam")))
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    jarFile = gatkJarFile;
    reference_sequence = b37
  }

  class cov (inBam: File, outRecalFile: File, dbSNP: File, KG: File) extends CountCovariates {
    this.memoryLimit = 3
    if (dbSNP != null) {
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    }
    if (KG != null) {
      this.rodBind :+= RodBind("mask", "VCF", KG)
    }
    if (dbSNP == null && KG == null) {
      this.run_without_dbsnp_potentially_ruining_quality = true
    }
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
    this.useOriginalQualities = true
    this.reference_sequence = b37
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
  }

  class recal (inBam: File, inRecalFile: File, outBam: File) extends TableRecalibration {
    this.memoryLimit = 3
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.useOriginalQualities = true
    this.out = outBam
    this.reference_sequence = b37
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
  }

}
