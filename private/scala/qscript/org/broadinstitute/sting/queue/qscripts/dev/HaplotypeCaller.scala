package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class HaplotypeCallerScript extends QScript {
  qscript =>

  @Argument(shortName="out", doc="output file", required=true)
  var out: String = "."
  @Argument(shortName="ref", doc="ref file", required=true)
  var ref: String = "."
  @Argument(shortName="bam", doc="bam file", required=true)
  var bam: String = "."
  @Argument(shortName="interval", doc="interval file", required=true)
  var interval: String = "."
  @Argument(shortName="sc", doc="scatter count", required=false)
  var scatterCount: Int = 100
  @Argument(shortName="dr", doc="downsampling", required=false)
  var downsampling: Int = 1000
  @Argument(shortName="lowpass", doc="lowpass", required=false)
  var lowpass: Boolean = false

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    memoryLimit = 4;
  }

  def script = {
    val hc = new HaplotypeCaller with UNIVERSAL_GATK_ARGS
    hc.reference_sequence = new File(ref)
    hc.intervalsString ++= List(interval)
    hc.scatterCount = scatterCount
    hc.input_file :+= new File(bam)
    hc.o = new File(out + ".hc.vcf")
    hc.dr = downsampling
    hc.analysisName = "HaplotypeCaller"
    if(lowpass) {
      hc.stand_call_conf = 4.0
      hc.stand_emit_conf = 4.0
    }
    add(hc)

    val ug = new UnifiedGenotyper with UNIVERSAL_GATK_ARGS
    ug.reference_sequence = new File(ref)
    ug.intervalsString ++= List(interval)
    ug.scatterCount = scatterCount / 2
    ug.input_file :+= new File(bam)
    ug.o = new File(out + ".ug.vcf")
    ug.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    ug.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    ug.analysisName = "UnifiedGenotyper"
    if(lowpass) {
      ug.stand_call_conf = 4.0
      ug.stand_emit_conf = 4.0
    }
    add(ug)
  }
}
