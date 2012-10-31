package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import java.lang.Math
import org.broadinstitute.sting.utils.PathUtils
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.queue.extensions.gatk.BaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.HaplotypeCaller
import org.broadinstitute.sting.queue.extensions.gatk.DelocalizedBaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.gatk.filters.SingleReadGroupFilter
import org.broadinstitute.sting.queue.extensions.gatk.SingleReadGroup

class HaplotypeCallerHMMPerformance extends QScript {

  @Argument(shortName = "it", doc="iterations", required=false)
  val iterations: Int = 4

  @Argument(shortName = "maxAltAllele", doc="maxAltAllele", required=false)
  val maxAlt: Int = 3

  @Argument(shortName = "mem", doc="memory limit", required=false)
  val mem: Int = 4

  @Argument(shortName = "I", doc="input bam list", required=false)
  val INPUT_BAM_LIST: String = "/humgen/gsa-hpprojects/dev/rpoplin/hmm_timing/AFR191.chr20.recal.bam.list"

  @Input(shortName="L", doc="interval file", required=true)
  var interval: List[File] = Nil

  @Input(shortName = "XL", doc="Exclude intervals list", required=false)
  var excludeIntervals: List[File] = Nil

  val b37 = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta"

  val BEFORE_JAR = new File("/humgen/gsa-hpprojects/dev/rpoplin/hmm_timing/before/dist/GenomeAnalysisTK.jar")
  val AFTER_JAR = new File("/humgen/gsa-hpprojects/dev/rpoplin/hmm_timing/after/dist/GenomeAnalysisTK.jar")
  val GATKs: Map[String, File] = Map(
    "before" -> BEFORE_JAR,
    "after" -> AFTER_JAR)

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = new File(b37)
    this.memoryLimit = mem
  }

  def script() {
    for ( iteration <- 0 until iterations ) {
      for ( (gatkName, gatkJar) <- GATKs ) {
          enqueueHC(iteration, gatkName, gatkJar)
      }
    }
  }

  def enqueueHC( iteration: Int, gatkName: String, gatkJar: File ) {
    trait VersionOverrides extends CommandLineGATK {
      this.configureJobReport(Map( "iteration" -> iteration, "jar" -> gatkName ))
    }

    val hc = new HaplotypeCaller with UNIVERSAL_GATK_ARGS with VersionOverrides
    hc.jarFile = gatkJar
    hc.intervals = interval
    hc.excludeIntervals = excludeIntervals
    hc.input_file :+= new File(INPUT_BAM_LIST)
    hc.o = new File("/dev/null")
    hc.analysisName = "HaplotypeCaller"
    hc.stand_call_conf = 10.0
    hc.stand_emit_conf = 10.0
    hc.minPruning = 2
    hc.max_alternate_alleles = maxAlt
    add(hc)

    val ug = new UnifiedGenotyper with UNIVERSAL_GATK_ARGS with VersionOverrides
    ug.jarFile = gatkJar
    ug.intervals = interval
    ug.excludeIntervals = excludeIntervals
    ug.input_file :+= new File(INPUT_BAM_LIST)
    ug.o = new File("/dev/null")
    ug.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    ug.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE
    ug.analysisName = "UnifiedGenotyper"
    ug.stand_call_conf = 8.0
    ug.stand_emit_conf = 8.0
    ug.max_alternate_alleles = maxAlt
    add(ug)
  }
}