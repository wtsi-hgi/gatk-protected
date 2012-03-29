import org.broadinstitute.sting.gatk.walkers.bqsr.RecalDataManager
import org.broadinstitute.sting.queue.extensions.gatk.BaseQualityScoreRecalibrator
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils

class BQSR extends QScript {

  @Argument(shortName = "b", required=true,  doc = "List of BAM files")  var bamList: File = _
  @Argument(shortName = "i", required=false, doc = "Intervals file")     var intervalsFile: List[File] = Nil
  @Argument(shortName = "r", required=false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "s", required=false, doc = "scatter/gather")     var scatterCount: Int = 200
  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="d", required=true) var dbSNP: Seq[File] = Seq()

  def script {
    val bams = QScriptUtils.createSeqFromFile(bamList);
    for (bam <- bams) {  
      val walker = new BaseQualityScoreRecalibrator();
      walker.reference_sequence = referenceFile
      walker.intervalsString = intervalsFile
      walker.out = swapExt(bam, ".bam", ".grp")
      walker.knownSites ++= dbSNP
      walker.input_file :+= bam
      walker.memoryLimit = 4
      walker.solid_nocall_strategy = RecalDataManager.SOLID_NOCALL_STRATEGY.PURGE_READ
      walker.scatterCount = scatterCount
      add(walker)
    }
  }
}
