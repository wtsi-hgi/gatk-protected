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
  @Argument(shortName="recalFile", doc="recal file", required=true)
  var recalFile: String = "."
  @Argument(shortName="sc", doc="scatter count", required=false)
  var scatterCount: Int = 135
  @Argument(shortName="dr", doc="downsampling", required=false)
  var downsampling: Int = 1000

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    memoryLimit = 2;
  }

  def script = {
    val hc = new HaplotypeCaller with UNIVERSAL_GATK_ARGS
    hc.reference_sequence = new File(ref)
    hc.intervalsString ++= List(interval)
    hc.scatterCount = scatterCount
    hc.input_file :+= new File(bam)
    hc.recalFile = new File(recalFile)
    hc.o = new File(out + ".hc.vcf")
    hc.dr = downsampling
    hc.analysisName = "HaplotypeCaller"
    add(hc)

    val ug = new UnifiedGenotyper with UNIVERSAL_GATK_ARGS
    ug.reference_sequence = new File(ref)
    ug.intervalsString ++= List(interval)
    ug.scatterCount = scatterCount
    ug.input_file :+= new File(bam)
    ug.o = new File(out + ".ug.vcf")
    ug.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    ug.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    ug.multiallelic = true
    ug.analysisName = "UnifiedGenotyper"
    add(ug)
  }
}
