package org.broadinstitute.sting.queue.qscripts.dev

import org.broadinstitute.sting.queue.extensions.gatk.BaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.HaplotypeCaller
import org.broadinstitute.sting.queue.extensions.gatk.DelocalizedBaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.gatk.filters.SingleReadGroupFilter
import org.broadinstitute.sting.queue.extensions.gatk.SingleReadGroup

class delocalizedBQSR_comparison_experiment extends QScript {

  @Argument(shortName = "i",  required=false, doc = "Intervals file")              var intervalsFile: List[File] = Nil
  @Argument(shortName = "ai",  required=false, doc = "aIntervals file")              var aintervalsFile: List[File] = Nil
  @Argument(shortName="alleles", doc="alleles file", required=true)
  var allelesFile: String = "."
  @Argument(shortName="allelesName", doc="alleles name", required=false)
  var allelesName: String = "indels"

  def script() {

    for( readGroup: String <- List("20FUK.1","20FUK.6","20FUK.7","20GAV.3","20GAV.5","20GAV.8") ) {
      for( contextSize: Int <- List(7)) {
        for( delocalized: Int <- List(0, 1) ) {

          if( delocalized == 0 ) {
            add(MappingBQSR(readGroup, contextSize, delocalized))
          } else {
            add(DelocalizedBQSR(readGroup, contextSize, delocalized))
          }
          add(HC(readGroup, contextSize, delocalized))

        }
      }
    }
  }

case class MappingBQSR( readGroup: String, contextSize: Int, delocalized: Int ) extends BaseRecalibrator with SingleReadGroup {
this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
this.intervalsString = intervalsFile
this.out = new File("/humgen/gsa-hpprojects/dev/rpoplin/bqsr/delocalized/NA12878."+readGroup+".context"+contextSize.toString()+".delocalized"+delocalized.toString()+".recal.grp")
this.input_file :+= new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam")
this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf"))
this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_135.b37.vcf"))
this.knownSites ++= List(new File("/humgen/gsa-hpprojects/dev/carneiro/bqsr/data/projectConsensus.snps.vcf"))
this.memoryLimit = 7
this.qq = 0
this.mcs = contextSize
this.ics = contextSize
this.goodRG = readGroup
this.scatterCount = 35
this.javaGCThreads = 4
}

case class DelocalizedBQSR( readGroup: String, contextSize: Int, delocalized: Int ) extends DelocalizedBaseRecalibrator with SingleReadGroup {
this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
this.intervalsString = intervalsFile
this.out = new File("/humgen/gsa-hpprojects/dev/rpoplin/bqsr/delocalized/NA12878."+readGroup+".context"+contextSize.toString()+".delocalized"+delocalized.toString()+".recal.grp")
this.input_file :+= new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam")
this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf"))
this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_135.b37.vcf"))
this.knownSites ++= List(new File("/humgen/gsa-hpprojects/dev/carneiro/bqsr/data/projectConsensus.snps.vcf"))
this.memoryLimit = 7
this.qq = 0
this.mcs = contextSize
this.ics = contextSize
this.goodRG = readGroup
this.scatterCount = 35
this.javaGCThreads = 4
}

case class HC( readGroup: String, contextSize: Int, delocalized: Int ) extends HaplotypeCaller {
this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
this.intervalsString = aintervalsFile
this.BQSR = new File("/humgen/gsa-hpprojects/dev/rpoplin/bqsr/delocalized/NA12878."+readGroup+".context"+contextSize.toString()+".delocalized"+delocalized.toString()+".recal.grp")
this.out = new File("/humgen/gsa-hpprojects/dev/rpoplin/bqsr/delocalized/calls/NA12878."+readGroup+".context"+contextSize.toString()+".delocalized"+delocalized.toString()+".recal.hc.vcf")
this.input_file :+= new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam")
this.memoryLimit = 4
this.keepRG = readGroup
this.qq = 0
this.scatterCount = 50
this.javaGCThreads = 4
this.alleles = new File(allelesFile)
this.minPruning = 1
this.out_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
this.gt_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
this.stand_emit_conf = 0.0
this.stand_call_conf = 0.0


}

}
