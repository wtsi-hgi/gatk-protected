package org.broadinstitute.sting.queue.qscripts.HaplotypeCalling

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.utils.interval.IntervalSetRule

class SingleExomeCalling extends QScript {
  @Argument(shortName = "test", doc="results", required=false)
  val test: Boolean = false

  @Argument(shortName = "debugHC", doc="results", required=false)
  val debugHC: Boolean = false

  val NA12878_COLLECTION = "/humgen/gsa-hpprojects/NA12878Collection/"
  val NA12878_EXOME_BAM = NA12878_COLLECTION + "bams/CEUTrio.HiSeq.WEx.b37.NA12878.clean.dedup.recal.bam"
  val NA12878_BEST_CALLS = NA12878_COLLECTION + "callsets/CEUtrio_BestPractices/CEUTrio.HiSeq.WGS.b37.snps_and_indels.recalibrated.filtered.phased.CURRENT.vcf"
  val EXOME_INTERVALS = "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"

  val b37Reference = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta")
  val DBSNP = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_135.b37.excluding_sites_after_129.vcf")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = b37Reference
    this.intervals = List(EXOME_INTERVALS)
    if ( test ) {
      //this.intervalsString = List("20")
      this.intervalsString = List("20:10,000,000-11,010,000")
    } else {
      this.intervalsString = List("20")
    }

    this.interval_set_rule = IntervalSetRule.INTERSECTION
    this.memoryLimit = 4
  }

  def script() {
    val bqsr = new BaseRecalibrator with UNIVERSAL_GATK_ARGS
    bqsr.input_file :+= NA12878_EXOME_BAM
    bqsr.knownSites :+= DBSNP
    bqsr.out = "bqsr.gatkreport.txt"
    bqsr.num_cpu_threads_per_data_thread = 4
    add(bqsr)

    val hc = new HaplotypeCaller with UNIVERSAL_GATK_ARGS
    hc.input_file :+= NA12878_EXOME_BAM
    hc.out = "hc.vcf"
    hc.debug = debugHC
    add(hc)

    val hcBQSR = new HaplotypeCaller with UNIVERSAL_GATK_ARGS
    hcBQSR.input_file :+= NA12878_EXOME_BAM
    hcBQSR.BQSR = bqsr.out
    hcBQSR.out = "hc.bqsr.vcf"
    hcBQSR.debug = debugHC
    add(hcBQSR)

    val hcAR = new HaplotypeCaller with UNIVERSAL_GATK_ARGS
    hcAR.input_file :+= NA12878_EXOME_BAM
    hcAR.out = "hcAR.vcf"
    hcAR.AR = hcAR.intervals
    hcAR.debug = debugHC
    add(hcAR)

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
    cv.V :+= new TaggedFile(hcBQSR.out, "hcBQSR")
    cv.V :+= new TaggedFile(hcAR.out, "hcAR")
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
