import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.QScript

class Recal extends QScript {

  @Input(shortName = "recal", required=true, doc = "recalibration file")  var recalFile: File = _
  @Input(shortName = "b", required=true,  doc = "List of BAM files")      var bamList: List[File] = _
  @Input(shortName = "i", required=false, doc = "Intervals file")         var intervalsFile: List[File] = Nil
  @Input(shortName = "r", required=false, doc = "Reference sequence")     var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Output(shortName = "o", required=true,  doc = "output file name")      var output: File = _
  @Argument(shortName = "s", required=false, doc = "scatter/gather")      var scatterCount: Int = 200
  @Argument(shortName = "m", required=false, doc = "memory limit")        var memLimit: Int = 4
  @Argument(shortName = "qq", required=false, doc = "quantization lvls")  var nLevels: Int = -1


  def script {
    val walker = new PrintReads();
    walker.reference_sequence = referenceFile
    walker.intervalsString = intervalsFile
    walker.out = output
    walker.input_file = bamList
    walker.memoryLimit = memLimit
    walker.qq = nLevels
    walker.scatterCount = scatterCount
    add(walker)
  }
}
