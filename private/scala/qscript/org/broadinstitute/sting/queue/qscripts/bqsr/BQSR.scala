import org.broadinstitute.sting.queue.extensions.gatk.{TaggedFile, BaseQualityScoreRecalibrator}
import org.broadinstitute.sting.queue.QScript

class BQSR extends QScript {

  @Argument(shortName = "b", required=true,  doc = "List of BAM files")  var bamList: List[File] = _
  @Argument(shortName = "o", required=true,  doc = "output file name")   var output: File = _
  @Argument(shortName = "i", required=false, doc = "Intervals file")     var intervalsFile: List[File] = Nil
  @Argument(shortName = "r", required=false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "s", required=false, doc = "scatter/gather")     var scatterCount: Int = 200
  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="d", required=true) var dbSNP: Seq[File] = Seq()

  def script {
    val walker = new BaseQualityScoreRecalibrator();
    walker.reference_sequence = referenceFile
    walker.intervalsString = intervalsFile
    walker.out = output
    walker.knownSites ++= dbSNP
    walker.input_file = bamList
    walker.memoryLimit = 8
    walker.scatterCount = scatterCount
    add(walker)
  }
}
