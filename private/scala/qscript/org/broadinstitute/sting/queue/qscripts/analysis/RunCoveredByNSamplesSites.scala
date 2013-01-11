package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.CoveredByNSamplesSites


class RunCoveredByNSamplesSites extends QScript{

  @Argument(doc="chrs", shortName="C", fullName = "chr")
  var chrs: Seq[String] = Nil

  val outputDir = "/humgen/gsa-firehose/ReduceReads_v2.2fixed_JointCalling/MacArthur_RRv2_2_JointCalling/chrs/"

  def script() {

    for (chr <- chrs) {
      val coveredByNSamplesSites = new CoveredByNSamplesSites
      coveredByNSamplesSites.variant = outputDir + chr + "/MacArthur_RRv2_2_JointCalling.chr"+chr+".snps.recalibrated.bcf"
      coveredByNSamplesSites.out = outputDir + chr + "/" + chr +".intervals"
      add(coveredByNSamplesSites)
    }

  }
}