package org.broadinstitute.sting.queue.qscripts.reducedreads

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.exceptions.UserException

/**
 * Makes calls on the full bam and the reduced bam. Selects only the high quality and high
 * confidence sites from the full bam callset and checks if they are all present in the reduced bam callset.
 *
 * @author Mauricio Carneiro
 * @since 9/19/12
 */


class ReduceReadsQC extends QScript {

  @Argument(shortName = "R",       fullName = "reference", doc="ref", required=false) var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  @Argument(shortName = "rbam",    fullName = "reduced_bam", doc="If you provide a reduced bam, reducing will be skipped.", required=false) var reducedBam: File = null;
  @Argument(shortName = "fbam",    fullName = "full_bam", doc="If you provide a reduced bam, reducing will be skipped.", required=false) var fullBam: File = null;
  @Argument(shortName = "i",                              doc = "Intervals file", required=false)        var intervalsFile: File = _
  @Argument(shortName = "s",                              doc = "scatterCount", required=false)          var scatterCount: Int = 50

  /**
   * Builds the CommandLineFunctions that will be used to run this script and adds them to this.functions directly or using the add() utility method.
   */
  def script() {
    val fullCallSet = call(fullBam)
    val reducedCallSet = call(reducedBam)

    val sv = new SelectVariants() with UNIVERSAL_GATK_ARGS
    sv.variant = fullCallSet
    sv.discordance = reducedCallSet
    sv.out = new File("missing_in_reduced.vcf")

    add(sv)
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = referenceFile
    this.intervals :+= intervalsFile
    this.memoryLimit = 4
  }


  def call(bam: File): File = {
    val callsVCF: File = getVCFName(bam)
    val filteredVCF: File = swapExt(callsVCF, ".vcf", ".filtered.vcf")

    val ug = new UnifiedGenotyper() with UNIVERSAL_GATK_ARGS
    ug.input_file :+= bam
    ug.out = callsVCF
    ug.scatterCount = scatterCount
    ug.nt = 2
    ug.stand_call_conf = 30.0
    ug.stand_emit_conf = 30.0

    val filter = new SelectVariants() with UNIVERSAL_GATK_ARGS
    filter.variant = callsVCF
    filter.out = filteredVCF
    filter.select :+= "QD > 10 && DP > 50 && FS < 10"

    add(ug, filter)

    return filteredVCF
  }

  def getVCFName (bam: File) : File = {
    if (bam.endsWith(".bam")) {
      swapExt(bam, ".bam", ".vcf")
    } else if (bam.endsWith(".list")) {
      swapExt(bam, ".list", ".vcf")
    } else {
      throw new UserException("input file must be a BAM (ending with .bam) or a list of bam files (ending with .list)")
    }
  }

}