import org.broadinstitute.sting.gatk.walkers.bqsr.RecalDataManager
import org.broadinstitute.sting.queue.extensions.gatk.BaseQualityScoreRecalibrator
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils

class BQSR extends QScript {

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="d", required=true) var dbSNP: Seq[File] = Seq()
  @Input(shortName = "recal", required=false, doc = "recalibration report") var recalFile: File = null

  @Argument(shortName = "b",  required=true,  doc = "List of BAM files")        var bamList: File = _
  @Argument(shortName = "i",  required=false, doc = "Intervals file")           var intervalsFile: List[File] = Nil
  @Argument(shortName = "mcs",required=false, doc = "mismatches context size")  var mismatchesContextSize = 3
  @Argument(shortName = "ics",required=false, doc = "insertions context size")  var insertionsContextSize = 3
  @Argument(shortName = "dcs",required=false, doc = "deletions context size")   var deletionsContextSize = 3
  @Argument(shortName = "r",  required=false, doc = "Reference sequence")       var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "m",  required=false, doc = "memory limit")             var memLimit: Int = 4
  @Argument(shortName = "qq", required=false, doc = "quantization lvls")        var nLevels: Int = 0
  @Argument(shortName = "k",  required=false, doc = "keep intermediate files")  var keepIntermediates: Boolean = true
  @Argument(shortName = "s",  required=false, doc = "scatter/gather")           var scatterCount: Int = 50
  @Argument(shortName = "np",  required=false, doc = "don't generate plots")    var noPlots: Boolean = false




  def script() {
    val bams = QScriptUtils.createSeqFromFile(bamList)
    val recalibrate: Boolean = recalFile == null
    for (bam <- bams) {
      if (recalibrate) {
        recalFile = swapExt(bam, ".bam", ".original.grp")
        val bqsr = new BaseQualityScoreRecalibrator();
        bqsr.reference_sequence = referenceFile
        bqsr.intervalsString = intervalsFile
        bqsr.out = recalFile
        bqsr.knownSites ++= dbSNP
        bqsr.input_file :+= bam
        bqsr.memoryLimit = memLimit
        bqsr.qq = nLevels
        bqsr.mcs = mismatchesContextSize
        bqsr.ics = insertionsContextSize
        bqsr.dcs = deletionsContextSize
        bqsr.no_plots = true
        bqsr.solid_nocall_strategy = RecalDataManager.SOLID_NOCALL_STRATEGY.PURGE_READ
        bqsr.scatterCount = scatterCount
        add(bqsr)
      }

      if (!noPlots) {
        val recal = new BaseQualityScoreRecalibrator();
        recal.reference_sequence = referenceFile
        recal.intervalsString = intervalsFile
        recal.out = swapExt(bam, ".bam", ".recal.grp")
        recal.knownSites ++= dbSNP
        recal.input_file :+= bam
        recal.memoryLimit = memLimit
        recal.qq = nLevels
        recal.mcs = mismatchesContextSize
        recal.ics = insertionsContextSize
        recal.dcs = deletionsContextSize
        recal.solid_nocall_strategy = RecalDataManager.SOLID_NOCALL_STRATEGY.PURGE_READ
        recal.scatterCount = scatterCount
        recal.keep_intermediate_files = keepIntermediates
        recal.BQSR = recalFile
        add(recal)
      }
    }
  }
}
