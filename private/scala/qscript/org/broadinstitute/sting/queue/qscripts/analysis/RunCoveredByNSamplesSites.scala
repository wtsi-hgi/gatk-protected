package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.CoveredByNSamplesSites
import org.broadinstitute.sting.commandline.Argument


class RunCoveredByNSamplesSites extends QScript{
  qscript =>

  @Argument(doc="chrs", shortName="C", fullName = "chr")
  var chrs: Seq[String] = Nil

  @Argument(doc="mode (snps or indels)", shortName="mode", fullName = "mode")
  var mode: String = _

  val outputDir = "/humgen/gsa-firehose/ReduceReads_v2.2fixed_JointCalling/MacArthur_RRv2_2_JointCalling/chrs/"

  def script() {

    for (chr <- chrs) {
     val coveredByNSamplesSites = new CoveredByNSamplesSites
      coveredByNSamplesSites.variant = outputDir + chr + "/MacArthur_RRv2_2_JointCalling.chr"+chr+"."+mode+".recalibrated.bcf"
      coveredByNSamplesSites.out = outputDir + chr + "/" + chr +"."+ mode +".covered.intervals"
      coveredByNSamplesSites.reference_sequence = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
      coveredByNSamplesSites.scatterCount = 100
      add(coveredByNSamplesSites)
    }

  }
}