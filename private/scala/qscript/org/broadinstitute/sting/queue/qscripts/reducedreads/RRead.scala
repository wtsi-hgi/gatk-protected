package org.broadinstitute.sting.queue.qscripts.reducedreads


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import io.Source._
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.gatk.walkers.compression.reducereads.ReduceReadsWalker.DownsampleStrategy

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 9/28/11
 * Time: 4:57 PM
 */


class RRead extends QScript {
  @Argument(shortName = "ref",      required = true,  fullName = "reference", doc = "Reference sequence") protected val reference: File = null
  @Argument(shortName = "bam",      required = false, fullName = "bam_file", doc = "") protected val bam: File = null
  @Argument(shortName = "ls",       required = false, fullName = "list", doc ="") protected val list: File = null
  @Argument(shortName = "int",      required = false, fullName = "intervals", doc = "") protected val intervals: File = null
  @Argument(shortName = "sg",       required = false, fullName = "scatterCount", doc ="") protected val scatterCount = 50
  @Argument(shortName = "csmm",     required = false, fullName = "context_size_mismatches", doc = "") protected var contextSizeMismatches: Option[Int] = None
  @Argument(shortName = "csindels", required = false, fullName = "context_size_indels", doc = "") protected var contextSizeIndels: Option[Int] = None
  @Argument(shortName = "minmap",   required = false, fullName = "minimum_mapping_quality", doc ="") protected var minMappingQuality: Option[Int] = None
  @Argument(shortName = "mintail",  required = false, fullName = "minimum_tail_qualities", doc ="") protected var minTailQuality: Option[Byte] = None
  @Argument(shortName = "minvar",   required = false, fullName = "minimum_alt_proportion_to_trigger_variant", doc = "") protected var minAltProportionToTriggerVariant: Option[Double] = None
  @Argument(shortName = "mindel",   required = false, fullName = "minimum_del_proportion_to_trigger_variant", doc = "") protected var minIndelProportionToTriggerVariant: Option[Double] = None
  @Argument(shortName = "minqual",  required = false, fullName = "minimum_base_quality_to_consider", doc = "") protected var minBaseQual: Option[Int] = None
  @Argument(shortName = "ds",       required = false, fullName = "downsample_coverage", doc = "") protected var downsampleCoverage: Option[Int] = None
  @Argument(shortName = "dm",       required = false, fullName = "downsample_strategy", doc = "") protected var downsampleStrategy: DownsampleStrategy = null
  @Argument(shortName = "dl",       required = false, fullName = "", doc = "") protected var debugLevel: Option[Int] = None
  @Argument(shortName = "dr",       required = false, fullName = "", doc = "") protected var debugRead: String = ""
  @Argument(shortName = "e",        required = false, fullName = "expand_intervals", doc = "Expand each target in input intervals by the specified number of bases. By default set to 50 bases.") protected var expandIntervals: Int = 0

    trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = reference;
    this.memoryLimit = 4
  }

  def script = {

    var bamList: List[File] = List()

    if (bam != null)
      bamList :+= bam
    else if (list != null)
      for (file: String <- fromFile(list).getLines())
        bamList :+= new File(file)
    else
      throw new UserException("You have to provide either a BAM or a LIST of bams.")

    val flankIntervals: File = new File("expanded_interval.list")
    var intervalList = List(intervals)

    if (expandIntervals > 0) {
      val writeFlanks = new WriteFlankingIntervalsFunction
      writeFlanks.reference = reference
      writeFlanks.inputIntervals = intervals
      writeFlanks.flankSize = expandIntervals
      writeFlanks.outputIntervals = flankIntervals
      writeFlanks.jobOutputFile = writeFlanks.outputIntervals + ".out"
      add(writeFlanks)
      intervalList :+= flankIntervals
    }


    for (file <- bamList) {
      val reducedBAM = swapExt(file.getName, ".bam", ".reduced.bam")

      // reduce
      val rr = new ReduceReads() with UNIVERSAL_GATK_ARGS
      rr.input_file :+= file
      rr.out = reducedBAM
      rr.scatterCount = scatterCount
      if (intervals != null) rr.intervals = intervalList
      if (contextSizeMismatches != None) rr.context_size_mismatches = Some(contextSizeMismatches)
      if (contextSizeIndels != None) rr.context_size_indels = Some(contextSizeIndels)
      if (minMappingQuality != None) rr.minimum_mapping_quality = Some(minMappingQuality)
      if (minTailQuality != None) rr.minimum_tail_qualities = Some(minTailQuality)
      if (minAltProportionToTriggerVariant != None) rr.minimum_alt_proportion_to_trigger_variant = Some(minAltProportionToTriggerVariant)
      if (minIndelProportionToTriggerVariant != None) rr.minimum_del_proportion_to_trigger_variant = Some(minIndelProportionToTriggerVariant)
      if (minBaseQual != None) rr.minimum_base_quality_to_consider = Some(minBaseQual)
      if (debugLevel != None) rr.debuglevel = Some(debugLevel)
      if (!debugRead.isEmpty) rr.debugread = debugRead
      if (downsampleCoverage != None) rr.downsample_coverage = Some(downsampleCoverage)
      if (downsampleStrategy != null) rr.downsample_strategy = downsampleStrategy

      add(rr)
    }
  }
}
