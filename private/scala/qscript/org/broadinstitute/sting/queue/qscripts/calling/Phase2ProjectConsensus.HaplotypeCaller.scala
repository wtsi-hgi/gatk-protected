import org.broadinstitute.sting.queue.extensions.gatk.BaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.HaplotypeCaller
import org.broadinstitute.sting.queue.extensions.gatk.DelocalizedBaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.gatk.filters.SingleReadGroupFilter
import org.broadinstitute.sting.queue.extensions.gatk.SingleReadGroup

class Phase2ProjectConsensusHaplotypeCaller extends QScript {

  @Argument(shortName = "i",  required=false, doc = "Intervals file")
  var intervalsFile: List[File] = Nil
  @Argument(shortName="alleles", doc="alleles file", required=true)
  var allelesFiles: List[File] = Nil
  @Argument(shortName="I", doc="bam list file", required=true)
  var inputBamList: File = new File("a.bam")


  def script() {

    for( allelesFile: File <- allelesFiles ) {
      	add(HC(allelesFile))
    }
  }

case class HC( allelesFile: File ) extends HaplotypeCaller {
  this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  this.intervalsString = intervalsFile
  this.out = swapExt(allelesFile, ".vcf", ".gga.vcf")
  this.memoryLimit = 5
  this.scatterCount = 25
  this.javaGCThreads = 4
  this.alleles = allelesFile
  this.minPruning = 1
  this.out_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
  this.gt_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
  this.GENOTYPE_GIVEN_ALL_ALLELES_COMBINATORIAL = true
  this.stand_emit_conf = 0.0
  this.stand_call_conf = 0.0
  this.input_file :+= inputBamList
  this.maxAltAlleles = 10
}

}
