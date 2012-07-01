/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 6/28/12
 * Time: 4:13 PM
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class VCFTester extends QScript {
  // Create an alias 'qscript' to be able to access variables
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="VCF file for locations.", shortName="V")
  var vcfFile: File = _

  @Input(doc="Samples file.", shortName="sf", required=false)
  var samplesFile: File = _

  @Input(doc="Intervals file.", shortName="L", required=false)
  var intervalsFile: File = _

  @Input(doc="Number of clients.", shortName="c")
  var numClients: Int = _

  def script() {
    val selectVariants = new SelectVariants

    selectVariants.reference_sequence = referenceFile
    selectVariants.variant = vcfFile

    if (samplesFile != null ) {
      selectVariants.sf :+= samplesFile
      selectVariants.out = swapExt(qscript.samplesFile, "samples", "%d.vcf".format(numClients))
    }
    else {
      selectVariants.out = "no_samples.vcf"
    }

    selectVariants.memoryLimit = 4
    selectVariants.scatterCount = numClients

    if (intervalsFile != null) {
      selectVariants.intervals :+= intervalsFile
    }

    add(selectVariants)
  }
}


