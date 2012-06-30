/**
 * Quick use of GCContentByInterval walker using an unmerged list of intervals
 *
 * Takes a list of intervals and calculates the GC Content of each interval
 * separately so the GATK engine doesn't merge the intervals.
 *
 * @author carneiro
 * @since 2/2/12
 */


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._


class IntervalCoverage extends QScript {
  @Argument(shortName = "i", required = true, doc = "Intervals file") var intervalsFile: List[File] = _
  @Argument(shortName = "b", required = true, doc = "List of BAM files") var bamList: List[File] = _
  @Argument(shortName = "r", required = false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "d", required = false, doc = "dbSNP sequence") var dbSNPFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  @Argument(shortName = "g", required = false, doc = "RGs to combine") var rgList: List[String] = List()
  @Argument(shortName = "sc", required = false, doc = "Scatter count") var ScatterCount: Int = 50

  def script {

    for (interval <- intervalsFile) {
      val output = swapExt(interval, ".list", ".tbl")

      val walker = new CoverageByRG()
      walker.reference_sequence = referenceFile
      walker.intervalsString = List(interval)
      walker.out = output
      walker.D = dbSNPFile
      walker.input_file = bamList
      walker.groupRGs = rgList
      walker.memoryLimit = 4
      walker.scatterCount = ScatterCount
      add(walker)

    }
  }

}
