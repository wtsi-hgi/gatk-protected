package org.broadinstitute.sting.queue.qscripts.HaplotypeCalling

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.utils.interval.IntervalSetRule

class SingleExomeCalling extends QScript {
  @Argument(shortName = "test", doc="results", required=false)
  val test: Boolean = false

  val NA12878_COLLECTION = "/humgen/gsa-hpprojects/NA12878Collection/"
  val NA12878_EXOME_BAM = NA12878_COLLECTION + "bams/CEUTrio.HiSeq.WEx.b37.NA12878.clean.dedup.recal.bam"
  val NA12878_BEST_CALLS = NA12878_COLLECTION + "callsets/CEUtrio_BestPractices/CEUTrio.HiSeq.WGS.b37.snps_and_indels.recalibrated.filtered.phased.CURRENT.vcf"
  val EXOME_INTERVALS = "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"

  val b37Reference = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = b37Reference
    this.intervals = List(EXOME_INTERVALS)
    if ( test ) {
      this.intervalsString = List("20")
      //this.intervalsString = List("20:1-1,110,000")
    } else {
      this.intervalsString = List("20")
    }

    this.interval_set_rule = IntervalSetRule.INTERSECTION
    this.memoryLimit = 4
  }

  def script() {
    val hc = new HaplotypeCaller with UNIVERSAL_GATK_ARGS
    hc.input_file :+= NA12878_EXOME_BAM
    hc.out = "hc.vcf"
    hc.debug = true
    add(hc)

    val ug = new UnifiedGenotyper with UNIVERSAL_GATK_ARGS
    ug.input_file :+= NA12878_EXOME_BAM
    ug.out = "ug.vcf"
    ug.glm = GenotypeLikelihoodsCalculationModel.Model.BOTH
    add(ug)

    val best = new SelectVariants with UNIVERSAL_GATK_ARGS
    best.V = NA12878_BEST_CALLS
    best.sample_name = List("NA12878")
    best.env = true
    best.out = if ( test ) "test.best.vcf" else "best.vcf"
    add(best)

    val cv = new CombineVariants with UNIVERSAL_GATK_ARGS
    cv.V :+= new TaggedFile(hc.out, "hc")
    cv.V :+= new TaggedFile(ug.out, "ug")
    cv.V :+= new TaggedFile(best.out, "best")
    cv.out = if ( test ) "test.combined.vcf" else "combined.vcf"
    add(cv)

    val vtt = new VariantsToTable with UNIVERSAL_GATK_ARGS
    vtt.V :+= cv.out
    vtt.F = List("CHROM", "POS", "TYPE", "set")
    vtt.out = swapExt(cv.out, "vcf", "table")
    add(vtt)

    val sv = new SelectVariants with UNIVERSAL_GATK_ARGS
    sv.V = cv.out
    sv.select_expressions = List("set != \"Intersection\"")
    sv.out = if ( test ) "test.discordant.vcf" else "discordant.vcf"
    add(sv)
  }
}
