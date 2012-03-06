package org.broadinstitute.sting.queue.qscripts.bqsr

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class EvaluateQuantizedQuals extends QScript {
  @Argument(shortName = "resources", doc="resources", required=false)
  val resourcesDir: String = "/humgen/gsa-hpprojects/GATK/bundle/current/b37"

  @Argument(shortName = "test", doc="test", required=false)
  val TEST: Boolean = false;

  val INPUT_BAM_FILENAME = "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam"
  val INPUT_BAM_QUAL_HIST = "NA12878.HiSeq.b37.chr20.10_11mb.qual_dist.gatkreport.txt"
  val INPUT_BAM_ALLELES = "NA12878.omni.vcf"
  val b37_FILENAME = "human_g1k_v37.fasta"

  def makeResource(x: String): File = new File("%s/%s".format(resourcesDir, x))

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = makeResource(b37_FILENAME);
    this.memoryLimit = 4
    this.intervalsString :+= (if ( TEST ) "20:10,000,000-11,000,000" else "20" )
  }

  def script = {
    def makeOneEval(nLevels: Int, maybeUG: Option[UnifiedGenotyper] = None): UnifiedGenotyper = {
      val inputBAM = makeResource(INPUT_BAM_FILENAME)
      val outputRoot = "levels.%s".format(nLevels)
      val QQ = new QuantizeQuals with UNIVERSAL_GATK_ARGS
      QQ.input_file :+= inputBAM
      QQ.o = new File(outputRoot + ".bam")
      QQ.simplifyBAM = true
      QQ.nLevels = nLevels
      QQ.report = new File(outputRoot + ".quantize.report.txt")
      QQ.qualityHistogram = INPUT_BAM_QUAL_HIST
      add(QQ)

      val CGL = new CalibrateGenotypeLikelihoods with UNIVERSAL_GATK_ARGS
      CGL.input_file :+= QQ.o
      CGL.alleles = new File(INPUT_BAM_ALLELES)
      CGL.out = new File(outputRoot + ".gcl")
      add(CGL)

      val UG = new UnifiedGenotyper with UNIVERSAL_GATK_ARGS
      UG.input_file :+= QQ.o
      UG.dbsnp = makeResource("dbsnp_132.b37.vcf")
      UG.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
      UG.out = new File(outputRoot + ".ug.vcf")
      add(UG)

      if ( ! maybeUG.isEmpty ) {
        val VE = new VariantEval with UNIVERSAL_GATK_ARGS
        VE.eval :+= new TaggedFile(UG.out, "eval")
        VE.out = new File(outputRoot + ".variant_eval.gatkreport.txt")
        VE.comp :+= new TaggedFile(maybeUG.get.out, "fullVCF")
        VE.dbsnp = new TaggedFile(maybeUG.get.out, "dbSNP")
        add(VE)
      }

      UG
    }

    val fullUG = makeOneEval(64)

    for ( nLevels <- if (TEST) List(4) else List(1, 2, 4, 8, 16, 32) ) {
      makeOneEval(nLevels, Some(fullUG))
    }
  }
}
