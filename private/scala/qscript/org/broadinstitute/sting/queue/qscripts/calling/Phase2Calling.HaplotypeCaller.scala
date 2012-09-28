package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class Phase2CallingHaplotypeCaller extends QScript {
  qscript =>

  @Argument(shortName="out", doc="output file", required=true)
  var out: String = "."
  @Argument(shortName="ref", doc="ref file", required=true)
  var ref: String = "."
  @Argument(shortName="bam", doc="bam file", required=true)
  var bam: String = "."
  @Argument(shortName="original", doc="bam file", required=true)
  var original: String = "."
  @Argument(shortName="mem", doc="mem", required=true)
  var mem: Double = 4
  @Argument(shortName="interval", doc="interval file", required=true)
  var interval: List[File] = Nil
  @Input(doc="Exclude intervals list", shortName = "XL", required=false)
  var excludeIntervals: List[File] = Nil

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
  }

  def script = {
    val hc = new HaplotypeCaller with UNIVERSAL_GATK_ARGS
    hc.reference_sequence = new File(ref)
    hc.intervals = interval
    hc.excludeIntervals = excludeIntervals
    hc.scatterCount = 790
    hc.input_file :+= new File(bam)
    hc.o = new File(out + ".hc.vcf")
    hc.analysisName = "HaplotypeCaller"
    hc.stand_call_conf = 8.0
    hc.stand_emit_conf = 8.0
    hc.minPruning = 1
    hc.javaGCThreads = 4
    hc.memoryLimit = Some(mem)
    add(hc)

    val ug = new UnifiedGenotyper with UNIVERSAL_GATK_ARGS
    ug.reference_sequence = new File(ref)
    ug.intervals = interval
    ug.excludeIntervals = excludeIntervals
    ug.scatterCount = 65
    ug.input_file :+= new File(original)
    ug.o = new File(out + ".ug.vcf")
    ug.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    ug.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    ug.analysisName = "UnifiedGenotyper"
    ug.stand_call_conf = 8.0
    ug.stand_emit_conf = 8.0
    ug.javaGCThreads = 4
    ug.memoryLimit = Some(mem)
    add(ug)
  }
}
