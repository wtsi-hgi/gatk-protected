package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import org.broadinstitute.sting.utils.baq.BAQ

class LeftAlignVariantsPipeline extends QScript {

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "vcf", doc="The file to left align", required=true)
  val vcfToLA : File = _

  @Argument(shortName= "o", doc="the output file",required=true)
  val outFile : File = _

  @Argument(shortName = "sc", doc = "the scatter count", required = false)
  val scatterCount : Int = 10

  val lav : LeftAlignVariants = new LeftAlignVariants
  lav.reference_sequence = referenceFile
  lav.variant = vcfToLA
  lav.out = outFile
  lav.scatterCount = scatterCount

}