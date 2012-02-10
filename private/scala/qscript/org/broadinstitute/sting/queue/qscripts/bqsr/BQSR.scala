import org.broadinstitute.sting.queue.extensions.gatk.{TaggedFile, BaseQualityScoreRecalibration}
import org.broadinstitute.sting.queue.QScript

class BQSR extends QScript {

  @Argument(shortName = "b", required=true,  doc = "List of BAM files")  var bamList: List[File] = _
  @Argument(shortName = "o", required=true,  doc = "output file name")   var output: File = _
  @Argument(shortName = "i", required=false, doc = "Intervals file")     var intervalsFile: List[File] = Nil
  @Argument(shortName = "r", required=false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "d", required=false, doc = "dbSNP sequence")     var dbSNPFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  @Argument(shortName = "s", required=false, doc = "scatter/gather")     var scatterCount: Int = 100

  def script {
    val walker = new BaseQualityScoreRecalibration();
    walker.reference_sequence = referenceFile
    walker.intervalsString = intervalsFile
    walker.out = output
    walker.knownSites :+= new TaggedFile(dbSNPFile, "dbSNP")
    walker.input_file = bamList
    walker.memoryLimit = 8
    walker.scatterCount = 50
    add(walker)
  }
}
