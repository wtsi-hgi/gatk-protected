package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import java.io.File


class RunJunctionGenotyper extends QScript {

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "r", doc="refSeq", required=false)
  var refSeqFile: File = new File("/humgen/gsa-hpprojects/GATK/data/refGene_b37.sorted.txt")

  @Argument(shortName = "bam", doc="BAM", required=true)
  val bam: File = null;

  @Argument(shortName = "o", doc = "output",required=true)
  val output: File = null;

  @Argument(shortName= "L", doc="interval list",required=false)
  val intervals: File = new File("/humgen/gsa-hphome1/chartl/projects/oneoffs/targets/all_genes.b37.interval_list")

   trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    memoryLimit = 4;
    this.reference_sequence = referenceFile
    this.input_file :+= bam
  }

  def bai(bam: File) = new File(bam + ".bai")

  // 1.) Unified Genotyper Base
  class GenotyperBase (t: File) extends IntronLossGenotyperV2 with UNIVERSAL_GATK_ARGS {
    this.intervals :+= t
    this.scatterCount = 20
    this.refSeq = new TaggedFile(refSeqFile,"REFSEQ")
  }

  class InsertSizeDistributionBase extends InsertSizeDistribution with UNIVERSAL_GATK_ARGS {
    this.intervalsString :+= "20"
  }

  def script = {
    val insertSize : InsertSizeDistribution = new InsertSizeDistributionBase()
    insertSize.out = swapExt(bam,".bam",".isd.table")
    add(insertSize)
    val genotype : IntronLossGenotyperV2 = new GenotyperBase(intervals)
    genotype.insertHistogram = insertSize.out
    genotype.out = output
    add(genotype)
  }
}
