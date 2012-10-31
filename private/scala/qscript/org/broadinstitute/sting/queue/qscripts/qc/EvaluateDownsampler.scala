package org.broadinstitute.sting.queue.qscripts.qc

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import java.lang.Math
import org.broadinstitute.sting.utils.PathUtils
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.utils.baq.BAQ

class EvaluateDownsampler extends QScript {
  val BAM = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam"
  val OMNI = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_omni2.5.b37.sites.vcf"
  val MILLS = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf"
  val INTERVAL = "1"

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta"
    this.memoryLimit = 4
    this.intervalsString = List(INTERVAL)
  }

  def script() {
    val orig = new Call("orig")
    val expt = new Call("expt")
    expt.enable_experimental_downsampling = true

    val CV = new CombineVariants with UNIVERSAL_GATK_ARGS
    CV.V +:= new TaggedFile(orig.out, "orig")
    CV.V +:= new TaggedFile(expt.out, "expt")
    CV.V +:= new TaggedFile(OMNI, "omni")
    CV.V +:= new TaggedFile(MILLS, "mills")
    CV.out = new File("combined.vcf")

    val VT = new VariantsToTable with UNIVERSAL_GATK_ARGS
    VT.V = List(CV.out)
    VT.F = List("CHROM", "POS", "set", "TYPE")
    VT.out = "combined.table"

    add(orig, expt, CV, VT)
  }

  class Call(val prefix: String) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file :+= BAM
    this.stand_call_conf = 50.0
    this.glm = GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.dcov = 60
    this.G = List()
    this.nosl = true
    this.o = new File(prefix + ".vcf")
    this.log = new File(prefix + ".log")
    this.nt = 8
  }
}
