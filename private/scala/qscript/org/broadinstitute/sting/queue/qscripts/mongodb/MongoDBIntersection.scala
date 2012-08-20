/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 6/27/12
 * Time: 3:00 PM
 * To change this template use File | Settings | File Templates.
 */
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.interval.IntervalSetRule

/**
 * Tests MongoDB by running SelectVariantsFromMongo in parallel
 * This edition uses intervals and scatter/gather to split up the processing
 */
class MongoDBIntersection extends QScript {
  // Create an alias 'qscript' to be able to access variables
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="VCF file for locations.", shortName="V")
  var vcfFile: File = _

  @Input(doc="Samples file.", shortName="sf", required=false)
  var samplesFile: File = _

  @Input(doc="Intervals file.", shortName="L", required=true)
  var intervalsFile: File = _

  @Input(doc="Number of clients.", shortName="c")
  var numClients: Int = _

  def script() {
    val selectVariants = new SelectVariantsFromMongo

    selectVariants.unsafe = ValidationExclusion.TYPE.LENIENT_VCF_PROCESSING
    selectVariants.reference_sequence = referenceFile
    selectVariants.variant = vcfFile

    val intersection = true

    if (samplesFile != null ) {
      selectVariants.sf :+= samplesFile
      if (intersection) {
        selectVariants.out = swapExt(qscript.samplesFile, "samples", "%d_intersection.vcf".format(numClients))
      }
      else {
        selectVariants.out = swapExt(qscript.samplesFile, "samples", "%d.vcf".format(numClients))
      }
    }
    else {
      selectVariants.no_samples = true
      selectVariants.out = "no_samples.vcf"
    }

    selectVariants.intervals :+= intervalsFile

    if (intersection) {
      selectVariants.intervals :+= "/home/unix/thibault/scratch/intervals/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"
      selectVariants.isr = IntervalSetRule.INTERSECTION
    }

    selectVariants.memoryLimit = 4
    selectVariants.scatterCount = numClients

    add(selectVariants)
  }
}
