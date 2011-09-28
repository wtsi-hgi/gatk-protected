package org.broadinstitute.sting.queue.qscripts.reducedreads


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import net.sf.samtools.SAMUtils


class ReduceReads extends QScript {
  @Argument(shortName = "reference", doc = "Reference sequence", required=true)
  var reference: File = null

  @Argument(shortName = "bam", doc = "", required=true)
  val bam: File = null

  @Argument(shortName = "intervals", doc = "", required=true)
  val intervals: File = null

  @Argument(shortName = "scatterCount", doc ="", required = false)
  val scatterCount = 50

  @Argument(fullName = "context_size", shortName = "cs", doc = "", required = false) protected var contextSize: Option[Int] = None
  @Argument(fullName = "minimum_mapping_quality", shortName = "minmap", doc = "", required = false) protected var minMappingQuality: Option[Int] = None
  @Argument(fullName = "minimum_tail_qualities", shortName = "mintail", doc = "", required = false) protected var minTailQuality: Option[Byte] = None
  @Argument(fullName = "minimum_alt_proportion_to_trigger_variant", shortName = "minvar", doc = "", required = false) protected var minAltProportionToTriggerVariant: Option[Double] = None
  @Argument(fullName = "minimum_del_proportion_to_trigger_variant", shortName = "mindel", doc = "", required = false) protected var minIndelProportionToTriggerVariant: Option[Double] = None
  @Argument(fullName = "minimum_base_quality_to_consider", shortName = "minqual", doc = "", required = false) protected var minBaseQual: Option[Int] = None
  @Argument(fullName = "maximum_consensus_base_qual", shortName = "maxqual", doc = "", required = false) protected var maxQualCount: Option[Byte] = None
  @Argument(fullName = "", shortName = "dl", doc = "", required = false) protected var debugLevel: Option[Int] = None
  @Argument(fullName = "", shortName = "dr", doc = "", required = false) protected var debugRead: String = ""
  @Argument(fullName = "downsample_coverage", shortName = "ds", doc = "", required = false) protected var downsampleCoverage: Option[Int] = None

    trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = reference;
    this.memoryLimit = 4
  }

  def script = {

    val reducedBAM = swapExt(bam, ".bam", ".reduced.bam")

    // reduce
    val rr = new org.broadinstitute.sting.queue.extensions.gatk.ReduceReads() with UNIVERSAL_GATK_ARGS
    rr.input_file :+= bam
    rr.out = reducedBAM
    rr.intervals :+= intervals
    rr.scatterCount = scatterCount
    if (contextSize != None) rr.context_size = Some(contextSize)
    if (minMappingQuality != None) rr.minimum_mapping_quality = Some(minMappingQuality)
    if (minTailQuality != None) rr.minimum_tail_qualities = Some(minTailQuality)
    if (minAltProportionToTriggerVariant != None) rr.minimum_alt_proportion_to_trigger_variant = Some(minAltProportionToTriggerVariant)
    if (minIndelProportionToTriggerVariant != None) rr.minimum_del_proportion_to_trigger_variant = Some(minIndelProportionToTriggerVariant)
    if (minBaseQual != None) rr.minimum_base_quality_to_consider = Some(minBaseQual)
    if (maxQualCount != None) rr.maximum_consensus_base_qual = Some(maxQualCount)
    if (debugLevel != None) rr.debuglevel = Some(debugLevel)
    if (!debugRead.isEmpty) rr.debugread = debugRead
    if (downsampleCoverage != None) rr.downsample_coverage = Some(downsampleCoverage)

    add(rr)
  }
}
