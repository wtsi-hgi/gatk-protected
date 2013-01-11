/**
 * Created with IntelliJ IDEA.
 * User: ami
 * Date: 10/11/12
 * Time: 9:42 AM
 * To change this template use File | Settings | File Templates.
 */

package org.broadinstitute.sting.queue.qscripts.pipeline

import net.sf.picard.io.IoUtil
import org.broadinstitute.sting.pipeline.{PicardSample, PicardAggregationUtils}
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function._
import org.broadinstitute.sting.queue.QScript

class ReduceReadsScript extends QScript {
  @Input(doc="Tab separated squid projects and samples.", shortName="tsv", exclusiveOf="bamList", required=false)
  var projectSampleTsv: File = _

  @Input(doc="BAM list files.", shortName="I", required=false)
  var externalBamList: File = _

  @Argument(doc="Subdirectory to store the reduced bams. By default set to 'reduced'.", shortName="bamDir", required=false)
  var bamDir = "reducedBAMs/"

  @Argument(doc="Reduce reads memory limit.", shortName="rrMem", required=false)
  var reduceReadsMemoryLimit = 4

  def script() {

    val reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
    val intervals = "/humgen/gsa-firehose2/ReduceReads_v2_JointCalling/gencode.v12_broad.agilent_merged.interval_list"

    var bams: Seq[Tuple2[File, File]] = Nil

    if (externalBamList != null) {
      for (originalBam: File <- io.Source.fromFile(externalBamList).getLines().toSeq.map(new File(_))) {
        val reducedBam: File = new File(new File(bamDir, "external"), swapExt(originalBam, ".bam", ".reduced.bam").getName)
        bams :+= Tuple2(originalBam, reducedBam)
      }
    }

    if (projectSampleTsv != null) {
      val picardSamples = PicardAggregationUtils.parseSamples(projectSampleTsv, false)

      for (picardSample: PicardSample <- picardSamples) {
        try {
          val picardIntervals = PicardAggregationUtils.readAnalysisIntervals(Seq(picardSample))
          require(reference == picardIntervals.getReference, "Unexpected reference: " + picardIntervals.getReference)

          val originalBam: File = PicardAggregationUtils.getSampleBam(picardSample.getProject, picardSample.getSample, picardSample.getVersion)
          val reducedBam: File = new File(bamDir, "%1$s/%2$s/v%3$d/%2$s.reduced.bam".format(
            IoUtil.makeFileNameSafe(picardSample.getProject),
            IoUtil.makeFileNameSafe(picardSample.getSample),
            picardSample.getVersion))

          //Use the hardcoded intervals instead of the Picard specified intervals
          bams :+= Tuple2(originalBam, reducedBam)
        } catch {
          case e =>
            println("Skipping: " + picardSample + ": " + e.getMessage());
        }
      }
    }

    for ((originalBam, reducedBam) <- bams) {
      val reduce = new ReduceReads with BadMate with RetryMemoryLimit
      reduce.memoryLimit = reduceReadsMemoryLimit
      reduce.reference_sequence = reference
      reduce.input_file = Seq(originalBam)
      reduce.intervals = Seq(new File(intervals))
      reduce.interval_padding = 50
      reduce.out = reducedBam
      add(reduce)
    }
  }
}