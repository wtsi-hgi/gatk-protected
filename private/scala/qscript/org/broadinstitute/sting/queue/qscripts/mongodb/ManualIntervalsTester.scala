/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 5/3/12
 * Time: 4:23 PM
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.interval.IntervalSetRule

/**
 * Tests MongoDB by running SelectVariantsFromMongo in parallel
 * This edition uses manually-created intervals files for scattering
 */
class ManualIntervalsTester extends QScript {
  // Create an alias 'qscript' to be able to access variables
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="VCF file for locations.", shortName="V")
  var vcfFile: File = _

  @Input(doc="Samples file.", shortName="sf", required=false)
  var samplesFile: File = _

  @Input(doc="Number of clients.", shortName="c")
  var numClients: Int = _

  def script() {
    for (i <- 1 to numClients) {
      val selectVariants = new SelectVariantsFromMongo

      selectVariants.reference_sequence = referenceFile
      selectVariants.variant = vcfFile

      if (samplesFile != null ) {
        selectVariants.sf :+= samplesFile
        selectVariants.out = swapExt(qscript.samplesFile, "samples", "%d.vcf".format(numClients))
      }
      else {
        selectVariants.sites_only = true
        selectVariants.out = "no_samples.vcf"
      }

      selectVariants.memoryLimit = 4

      //selectVariants.intervals :+= "../chr20_%d.interval_list".format(i)
      selectVariants.intervals :+= "/home/unix/thibault/scratch/intervals/LYPLAL1_%d_of_%d.interval_list".format(i, numClients)
      selectVariants.intervals :+= "/home/unix/thibault/scratch/intervals/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"
      selectVariants.isr = IntervalSetRule.INTERSECTION

      add(selectVariants)
    }
  }
}
