package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.baq.BAQ

class ReducedBAMEvaluation extends QScript {
  @Argument(shortName = "R",       fullName = "reference", doc="ref", required=false) var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  @Argument(shortName = "rbam",    fullName = "reduced_bam", doc="If you provide a reduced bam, reducing will be skipped.", required=false) var reducedBam: File = null;
  @Argument(shortName = "fbam",    fullName = "full_bam", doc="If you provide a reduced bam, reducing will be skipped.", required=false) var fullBam: File = null;
  @Argument(shortName = "ri",      fullName = "reduce_interval", doc="Interval to reduce at", required=false) var reduceInterval: String = null;
  @Argument(shortName = "cs",      fullName = "context_size", doc = "", required = false) protected var contextSize: Int = _
  @Argument(shortName = "minmap",  fullName = "minimum_mapping_quality", doc = "", required = false) protected var minMappingQuality: Int = _
  @Argument(shortName = "mintail", fullName = "minimum_tail_qualities", doc = "", required = false) protected var minTailQuality: Byte = _
  @Argument(shortName = "minvar",  fullName = "minimum_alt_proportion_to_trigger_variant", doc = "", required = false) protected var minAltProportionToTriggerVariant: Double = _
  @Argument(shortName = "mindel",  fullName = "minimum_del_proportion_to_trigger_variant", doc = "", required = false) protected var minIndelProportionToTriggerVariant: Double = _
  @Argument(shortName = "minqual", fullName = "minimum_base_quality_to_consider", doc = "", required = false) protected var minBaseQual: Int = _
  @Argument(shortName = "maxqual", fullName = "maximum_consensus_base_qual", doc = "", required = false) protected var maxQualCount: Byte = _
  @Argument(shortName = "ds",      fullName = "downsample_coverage", doc = "", required = false) protected var downsampleCoverage: Int = _

  val dbSNP_b37: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"



  def script = {



    val sliceBAM =  swapExt(bam,".bam",".printreads.bam")
    val sliceVCF = swapExt(sliceBAM,".bam",".filtered.vcf")
    add(SliceBAM(bam, sliceBAM))
    callAndEvaluateBAM(sliceBAM, sliceVCF)

    for {
      a <- cs
      b <- adav
      c <- mbrc
    } yield {
    val header = ".CS-" + toStringLength(cs, a) + ".mrav-" + toStringLength(adav, b) + ".MBRC-" + toStringLength(mbrc, c)
    val reduceBAM = swapExt(bam, ".bam", header + ".reduced.bam")
    val reduceVCF = swapExt(reduceBAM,".bam",".filtered.vcf")
    val combineVCF = swapExt(reduceVCF, ".bam", ".filtered.combined.vcf")

    // Generate the new BAMs
    add(ReduceBAM(bam, reduceBAM, a, b, c ))

    // Call SNPs, filter and get variant eval report
    callAndEvaluateBAM(reduceBAM, reduceVCF)

    // Combine the callsets
    add(combine(reduceVCF, sliceVCF, combineVCF))

    // Get a report of the combined callsets
    val eval = new Eval(combineVCF) // evaluate the combined VCF
    eval.select = List("'set==\"Intersection\"'", "'set==\"fullBAM\"'", "'set==\"reducedBAM\"'", "'set==\"filterInreducedBAM-fullBAM\"'", "'set==\"reducedBAM-filterInfullBAM\"'")
    eval.selectName = List("Intersection", "fullBAM", "reducedBAM", "filterInreducedBAM-fullBAM", "reducedBAM-filterInfullBAM")
    add(eval)
    }

  }
  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
  }


  case class ReduceBAM(bam: File, outVCF: File, a: Int, b: Int, c: Int) extends ReduceReads with UNIVERSAL_GATK_ARGS with CoFoJa {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.o = outVCF
    this.CS = a   // Best safe value 7
    this.ADAV = b     //Best 30
    this.MBRC = c
    this.baq = BAQ.CalculationMode.OFF
    this.intervalsString = List(reduceInterval);
  }

  case class SliceBAM(bam: File, outVCF: File) extends PrintReads with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.o = outVCF
    this.intervalsString = List(reduceInterval);
  }

  case class combine (reducedBAMVCF: File, fullBAMVCF: File, outVCF: File) extends CombineVariants with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("fullBAM", "VCF", fullBAMVCF)
    this.rodBind :+= RodBind("reducedBAM", "VCF", reducedBAMVCF)
    this.rod_priority_list = "reducedBAM,fullBAM"
    this.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    this.out = outVCF
  }

  def callAndEvaluateBAM(inBAM: File, outVCF: File) {
    val rawVCF   = swapExt(inBAM, ".bam", ".vcf")
    val recalVCF = swapExt(inBAM, ".bam", ".filtered.vcf")
    val recal    = swapExt(inBAM, ".bam", ".recal")
    val tranches = swapExt(inBAM, ".bam", ".tranches")

    add(Call(inBAM, rawVCF),
        HardFilter(rawVCF, recalVCF),
//        VQSR(rawVCF, recal, tranches),
//        applyVQSR(rawVCF, recal, tranches, recalVCF),
        Eval(recalVCF)
//        DiffableTable(rawVCF),   // for convenient diffing
//        DiffableTable(recalVCF)  // for convenient diffing
    )
  }

  case class Eval(@Input vcf: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("eval", "VCF", vcf)
    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.doNotUseAllStandardStratifications = true
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants")
    this.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "JexlExpression")
    this.out = swapExt(vcf,".vcf",".eval")
    this.intervalsString = List(reduceInterval);
  }

  case class Call(inBAM: File, outVCF: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file = List(inBAM)
    this.stand_call_conf = 50.0
    this.stand_emit_conf = 50.0
    this.dcov = DCOV;
    this.baq = BAQ.CalculationMode.OFF
    this.o = outVCF
    this.intervalsString = List(reduceInterval);

    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)

    if ( minimalVCF )
      this.group = List("none")

    if ( reduceInterval != null ) {
      this.intervalsString = List(reduceInterval)
    }
  }

  case class HardFilter (inVCF: File, outVCF: File) extends VariantFiltration with UNIVERSAL_GATK_ARGS {
    this.variantVCF = inVCF
    this.filterName = List("SNP_SB", "SNP_QD", "SNP_HRun")
    this.filterExpression = List("\"SB>=0.10\"", "\"QD<5.0\"", "\"HRun>=4\"")
    this.clusterWindowSize = 10
    this.clusterSize = 3
    this.out = outVCF
  }

  case class VQSR(inVCF: File, outRecal: File, outTranches: File) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.intervalsString = List(reduceInterval)
    this.rodBind :+= RodBind("input", "VCF", inVCF)
    this.rodBind :+= RodBind("hapmap", "VCF", hapmap_b37, "known=false,training=true,truth=true,prior=15.0")
    this.rodBind :+= RodBind("omni", "VCF", omni_b37, "known=false,training=true,truth=false,prior=12.0")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP_b37, "known=true,training=false,truth=false,prior=4.0")
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "DP")
    this.tranches_file = outTranches
    this.recal_file = outRecal
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.analysisName = outTranches + "_VQSR"
    this.jobName =  outTranches + ".VQSR"
    this.mG = 5
    this.std = 14
    this.percentBad = 0.04
    this.nt = 5
}

  // 4.) Apply the recalibration table to the appropriate tranches
  case class applyVQSR (inVCF: File, inRecal: File, inTranches: File, outVCF: File) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.intervalsString = List(reduceInterval)
    this.rodBind :+= RodBind("input", "VCF", inVCF )
    this.tranches_file = inTranches
    this.recal_file = inRecal
    this.ts_filter_level = 99.0
    this.out = outVCF
    this.analysisName = outVCF + "_AVQSR"
    this.jobName =  outVCF + ".applyVQSR"
  }


  case class DiffableTable(@Input vcf: File) extends CommandLineFunction {
    @Output var out: File = swapExt(vcf,".vcf",".table")
    def commandLine = "cut -f 1,2,4,5,7 %s > %s".format(vcf, out)
  }
}

