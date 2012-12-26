/**
 * Created with IntelliJ IDEA.
 * User: ami
 * Date: 10/11/12
 * Time: 9:42 AM
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function._
import org.broadinstitute.sting.queue.QScript

class ReduceReadsScript extends QScript {
  //@Input(doc="Tab separated squid projects and samples.", shortName="tsv", exclusiveOf="bamList", required=false)
  //var projectSampleTsv: File = _

  @Input(doc="BAM list files.", shortName="I", required=true)
  var externalBamList: File = _

  @Argument(doc="Subdirectory to store the reduced bams. By default set to 'reduced'.", shortName="bamDir", required=false)
  var bamDir = "reducedBAMs/"

  @Argument(doc="Reduce reads memory limit.", shortName="rrMem", required=false)
  var reduceReadsMemoryLimit = 4

  def script() {

    val reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
    val intervals = "/humgen/gsa-firehose2/ReduceReads_v2_JointCalling/gencode.v12_broad.agilent_merged.interval_list"

    var bams: Seq[Tuple3[File, File, File]] = Nil

    if (externalBamList != null) {
      for (originalBam: File <- io.Source.fromFile(externalBamList).getLines().toSeq.map(new File(_))) {
        val reducedBam: File = new File(new File(bamDir, "external"), swapExt(originalBam, ".bam", ".reduced.bam").getName)
        bams :+= Tuple3(originalBam, reducedBam, new File(intervals))
      }
    }

    //    if (projectSampleTsv != null) {
    //      val picardSamples = PicardAggregationUtils.parseSamples(projectSampleTsv, false)
    //
    //      for (picardSample <- picardSamples) {
    //        try {
    //          val picardIntervals = PicardAggregationUtils.readAnalysisIntervals(Seq(picardSample))
    //          require(reference == picardIntervals.getReference, "Unexpected reference: " + picardIntervals.getReference)
    //
    //          val originalBam: File = PicardAggregationUtils.getSampleBam(picardSample.getProject, picardSample.getSample, picardSample.getVersion)
    //          val reducedBam: File = new File(bamDir, "%1$s/%2$s/v%3$d/%2$s.reduced.bam".format(
    //            IoUtil.makeFileNameSafe(picardSample.getProject),
    //            IoUtil.makeFileNameSafe(picardSample.getSample),
    //            picardSample.getVersion))

    // Use the hardcoded intervals instead of the Picard specified intervals
    //          bams :+= Tuple3(originalBam, reducedBam, new File(intervals))
    //	} catch {
    //          case e =>
    //            println("Skipping: " + picardSample + ": " + e.getMessage());
    //        }
    //      }
    //    }

    trait ChromosomeIntervals extends CommandLineGATK {
      this.intervalsString :+= "1"
      this.interval_set_rule = org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
    }

    for ((originalBam, reducedBam, bamIntervals) <- bams) {
      val reduce = new ReduceReads with BadMate with RetryMemoryLimit with ChromosomeIntervals
      reduce.memoryLimit = reduceReadsMemoryLimit
      reduce.reference_sequence = reference
      reduce.input_file = Seq(originalBam)
      reduce.intervals = Seq(bamIntervals)
      reduce.interval_padding = 50
      reduce.out = reducedBam
      add(reduce)
    }
  }
}