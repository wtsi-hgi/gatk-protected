/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 5/3/12
 * Time: 4:23 PM
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * Tests MongoDB by running SelectVariantsFromMongo in parallel
 * This edition runs the full processing on each node
 */
class MultipleMongoDBTester extends QScript {
  // Create an alias 'qscript' to be able to access variables
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="VCF file for locations.", shortName="V")
  var vcfFile: File = _

  @Input(doc="Samples file.", shortName="sf")
  var samplesFile: File = _

  @Input(doc="Number of clients.", shortName="c")
  var numClients: Int = _

  def script() {
    for (i <- 1 to numClients) {
      val selectVariants = new SelectVariantsFromMongo

      selectVariants.reference_sequence = referenceFile
      selectVariants.variant = vcfFile
      selectVariants.sf :+= samplesFile
      selectVariants.memoryLimit = 2

      selectVariants.out = swapExt(qscript.samplesFile, "samples", "%d_client%d.vcf".format(numClients, i))

      add(selectVariants)
    }
  }
}
