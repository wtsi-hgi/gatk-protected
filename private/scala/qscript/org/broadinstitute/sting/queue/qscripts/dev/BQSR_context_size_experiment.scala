package org.broadinstitute.sting.queue.qscripts.dev

import org.broadinstitute.sting.gatk.walkers.bqsr.RecalDataManager
import org.broadinstitute.sting.queue.extensions.gatk.BaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.gatk.filters.SingleReadGroupFilter
import org.broadinstitute.sting.queue.extensions.gatk.SingleReadGroup

class BQSR_context_size_experiment extends QScript {

  @Argument(shortName = "i",  required=false, doc = "Intervals file")              var intervalsFile: List[File] = Nil

  def script() {

    for( readGroup: String <- List("20FUK.1","20FUK.6","20GAV.3","20GAV.8") ) {
      for( contextSize: Int <- List(3,4,5,6,7)) {
        add(BQSR(readGroup, contextSize))
      }
    }
  }

  case class BQSR( readGroup: String, contextSize: Int ) extends BaseRecalibrator with SingleReadGroup {
    this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
    this.intervalsString = intervalsFile
    this.out = new File("/humgen/gsa-hpprojects/dev/rpoplin/bqsr/NA12878."+readGroup+".context"+contextSize.toString()+".recal.grp")
    this.input_file :+= new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam")
    this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf"))
    this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf"))
    this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_135.b37.vcf"))
    this.memoryLimit = 8
    this.qq = 0
    this.mcs = contextSize
    this.ics = contextSize
    //this.rf = Seq("SingleReadGroup")
    this.goodRG = readGroup
    this.no_plots = true
    this.keep_intermediate_files = false
    this.scatterCount = 25
    this.javaGCThreads = 4
  }


}
